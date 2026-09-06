//! A world in PostgreSQL.
//!
//! # Why this is not a blob
//!
//! `FileStore` already puts a whole world in one place. Doing the same into a
//! `bytea` column would be a very expensive file: it would throw away the only
//! things a database is for — reading one node without reading the rest,
//! writing one node without rewriting the rest, and asking questions about
//! rows.
//!
//! So the node is the row. Its identity is its [`PathKey`], which the tree
//! already guarantees is stable across a node being destroyed and rebuilt,
//! which is exactly what a primary key has to be. The columns are the fields
//! something will eventually query on — parent, tier, depth, aliveness, the two
//! timestamps — and the rest of the node travels as `bytea` in the same wire
//! format a world file uses. That split is deliberate: indexing what you filter
//! on and blobbing what you only ever read whole is the shape that survives
//! being ported.
//!
//! # Writes are proportional to events
//!
//! [`WorldStore::save_since`] writes only the nodes whose dynamics were
//! re-derived or that something happened to. Everything else was *coasted*, and
//! a coasted node's position at any instant is a closed-form function of state
//! the database already holds — so the reader reconstructs it with
//! [`Snapshot::settle`] rather than the writer storing it again.
//!
//! This is the property that decides whether the design carries a large world,
//! and it is worth being precise about what it does and does not buy. It makes
//! write volume track *what happened* rather than *how much world there is*. It
//! does nothing about a single instance being a single instance — one Postgres
//! is one write bottleneck, deliberately accepted for now.
//!
//! # Swapping this out
//!
//! Everything Postgres-specific is in this file, behind the `postgres` feature,
//! and the crate has no dependencies without it. The contract another backend
//! has to meet is [`WorldStore`]: four methods, of which one has a working
//! default. The schema below is ordinary enough — a keyed table of rows with a
//! blob payload — that a document store or SpacetimeDB would express it almost
//! unchanged.

use crate::ids::{NodeIdx, PathKey};
use crate::persist::{dirty_nodes, Flushed, Snapshot, WorldStore, WorldView};
use crate::wire::{Reader, WireError, Writer};
use postgres::{Client, NoTls};

type Result<T> = std::result::Result<T, WireError>;

/// `postgres::Error`'s own `Display` is famously terse — a failed insert
/// renders as the four words "db error", and the column it actually objected to
/// is in the source chain underneath. Walking it turns an unreadable failure
/// into an actionable one, which mattered the first time a schema change met an
/// existing database.
fn db(e: postgres::Error) -> WireError {
    let mut msg = e.to_string();
    let mut src: Option<&(dyn std::error::Error + 'static)> = std::error::Error::source(&e);
    while let Some(inner) = src {
        msg.push_str(": ");
        msg.push_str(&inner.to_string());
        src = inner.source();
    }
    WireError::Io(msg)
}

/// A `u128` key as sixteen big-endian bytes, so it sorts in the database the
/// same way it compares in memory.
fn key_bytes(k: PathKey) -> [u8; 16] {
    k.0.to_be_bytes()
}

fn key_from(b: &[u8]) -> Result<PathKey> {
    if b.len() != 16 {
        return Err(WireError::Truncated { what: "path key", need: 16, have: b.len() });
    }
    let mut a = [0u8; 16];
    a.copy_from_slice(b);
    Ok(PathKey(u128::from_be_bytes(a)))
}

/// Bumped whenever the table layout changes in a way an existing database
/// would not satisfy.
///
/// The tables are created with `IF NOT EXISTS`, which is right for a fresh
/// database and silently wrong for one built against an older layout: the new
/// column is never added, and every insert afterwards fails with an error that
/// names a column rather than the real problem. So the layout carries its own
/// version, checked on connect, and a mismatch is refused the way the file
/// reader refuses an old format — by saying so, rather than by misbehaving.
pub const SCHEMA_VERSION: i32 = 2;

pub const SCHEMA: &str = r#"
CREATE TABLE IF NOT EXISTS schema_version (
    id      int PRIMARY KEY DEFAULT 1 CHECK (id = 1),
    version int NOT NULL
);

CREATE TABLE IF NOT EXISTS world (
    id             int PRIMARY KEY DEFAULT 1 CHECK (id = 1),
    format         int              NOT NULL,
    world_seed     bytea            NOT NULL,
    root           int              NOT NULL,
    instant        double precision NOT NULL,
    pace           double precision NOT NULL,
    time_rate      double precision NOT NULL,
    time_throttle  double precision NOT NULL,
    paced_to       bigint           NOT NULL,
    pace_fixed     boolean          NOT NULL,
    labour_rate    double precision NOT NULL,
    rejected       bigint           NOT NULL,
    tree_stats     bytea            NOT NULL
);

CREATE TABLE IF NOT EXISTS node (
    key             bytea            PRIMARY KEY,
    idx             int              NOT NULL,
    parent          int              NOT NULL,
    depth           int              NOT NULL,
    tier            smallint         NOT NULL,
    alive           boolean          NOT NULL,
    pinned          boolean          NOT NULL,
    epoch           int              NOT NULL,
    instant         double precision NOT NULL,
    last_solved     double precision NOT NULL,
    last_disturbed  double precision NOT NULL,
    payload         bytea            NOT NULL
);
CREATE UNIQUE INDEX IF NOT EXISTS node_idx     ON node(idx);
CREATE INDEX        IF NOT EXISTS node_parent  ON node(parent);
CREATE INDEX        IF NOT EXISTS node_tier    ON node(tier);
CREATE INDEX        IF NOT EXISTS node_solved  ON node(last_solved);

CREATE TABLE IF NOT EXISTS pinned (
    key    bytea PRIMARY KEY,
    bodies bytea NOT NULL
);

CREATE TABLE IF NOT EXISTS fact (
    key      bytea            NOT NULL,
    quantity smallint         NOT NULL,
    value    double precision NOT NULL,
    at       double precision NOT NULL,
    seq      bigint           NOT NULL,
    PRIMARY KEY (key, quantity)
);

CREATE TABLE IF NOT EXISTS ledger_counters (
    id       int PRIMARY KEY DEFAULT 1 CHECK (id = 1),
    sequence bigint NOT NULL,
    queries  bigint NOT NULL,
    commits  bigint NOT NULL
);

CREATE TABLE IF NOT EXISTS environment (
    key  bytea PRIMARY KEY,
    data bytea NOT NULL
);

CREATE TABLE IF NOT EXISTS in_flight (
    seq             bigint           PRIMARY KEY,
    arrives         double precision NOT NULL,
    target          int              NOT NULL,
    kind            smallint         NOT NULL,
    energy          double precision NOT NULL,
    px              double precision NOT NULL,
    py              double precision NOT NULL,
    pz              double precision NOT NULL,
    source_distance double precision NOT NULL
);
CREATE INDEX IF NOT EXISTS in_flight_arrives ON in_flight(arrives);

CREATE TABLE IF NOT EXISTS mailbox_counters (
    id             int PRIMARY KEY DEFAULT 1 CHECK (id = 1),
    delivered      bigint NOT NULL,
    in_flight_peak bigint NOT NULL
);

CREATE TABLE IF NOT EXISTS audit (
    seq          bigint           PRIMARY KEY,
    key          bytea            NOT NULL,
    property     smallint         NOT NULL,
    delta_energy double precision NOT NULL,
    at           double precision NOT NULL
);
"#;

pub struct PostgresStore {
    client: Client,
}

impl PostgresStore {
    /// Connect and make sure the schema is there.
    ///
    /// `url` is an ordinary libpq connection string, e.g.
    /// `host=127.0.0.1 port=5432 user=phys dbname=phys`.
    /// Connect, creating the tables if they are absent and refusing a database
    /// built against a layout this build cannot write.
    pub fn connect(url: &str) -> Result<PostgresStore> {
        let mut client = Client::connect(url, NoTls).map_err(db)?;
        // Is there anything here already? Ask before running the schema, since
        // `CREATE TABLE IF NOT EXISTS` would make an old database look new.
        let existing: Option<i32> = client
            .query_opt(
                "SELECT version FROM schema_version WHERE id = 1",
                &[],
            )
            .ok()
            .flatten()
            .map(|r| r.get(0));
        let has_tables: bool = client
            .query_one(
                "SELECT to_regclass('public.world') IS NOT NULL",
                &[],
            )
            .map_err(db)?
            .get(0);

        match existing {
            Some(v) if v != SCHEMA_VERSION => {
                return Err(WireError::UnsupportedVersion {
                    found: v as u16,
                    supported: SCHEMA_VERSION as u16,
                });
            }
            // Tables from before the version stamp existed, or from a build
            // that wrote a different layout. Either way this build cannot use
            // them, and saying so beats a failed insert naming one column.
            None if has_tables => {
                return Err(WireError::UnsupportedVersion {
                    found: 0,
                    supported: SCHEMA_VERSION as u16,
                });
            }
            _ => {}
        }

        client.batch_execute(SCHEMA).map_err(db)?;
        client
            .execute(
                "INSERT INTO schema_version (id, version) VALUES (1, $1)
                 ON CONFLICT (id) DO UPDATE SET version = EXCLUDED.version",
                &[&SCHEMA_VERSION],
            )
            .map_err(db)?;
        Ok(PostgresStore { client })
    }

    /// Drop every table and build them again at the current layout.
    ///
    /// The way out of the refusal above, and destructive by design: there is no
    /// migration path from an unknown layout, and pretending otherwise would be
    /// worse than saying what this does.
    pub fn reset(url: &str) -> Result<PostgresStore> {
        let mut client = Client::connect(url, NoTls).map_err(db)?;
        client
            .batch_execute(
                "DROP TABLE IF EXISTS node, pinned, fact, environment, audit, world,
                 ledger_counters, in_flight, mailbox_counters, schema_version CASCADE;",
            )
            .map_err(db)?;
        drop(client);
        PostgresStore::connect(url)
    }

    /// Throw the world away. Used by tests and by "start again".
    pub fn clear(&mut self) -> Result<()> {
        self.client
            .batch_execute(
                "TRUNCATE node, pinned, fact, environment, audit, world, ledger_counters,
                  in_flight, mailbox_counters;",
            )
            .map_err(db)
    }

    pub fn client(&mut self) -> &mut Client {
        &mut self.client
    }

    /// How many node rows exist.
    pub fn node_count(&mut self) -> Result<i64> {
        let row = self.client.query_one("SELECT count(*) FROM node", &[]).map_err(db)?;
        Ok(row.get(0))
    }

    fn write_world_row(&mut self, v: &WorldView<'_>) -> Result<()> {
        let mut stats = Writer::new();
        crate::persist::put_tree_stats(&mut stats, &v.tree.stats);
        self.client
            .execute(
                "INSERT INTO world (id, format, world_seed, root, instant, pace, time_rate,
                                    time_throttle, paced_to, pace_fixed, labour_rate, rejected, tree_stats)
                 VALUES (1, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12)
                 ON CONFLICT (id) DO UPDATE SET
                    format = EXCLUDED.format, world_seed = EXCLUDED.world_seed,
                    root = EXCLUDED.root, instant = EXCLUDED.instant, pace = EXCLUDED.pace,
                    time_rate = EXCLUDED.time_rate, time_throttle = EXCLUDED.time_throttle,
                    paced_to = EXCLUDED.paced_to, pace_fixed = EXCLUDED.pace_fixed,
                    labour_rate = EXCLUDED.labour_rate,
                    rejected = EXCLUDED.rejected, tree_stats = EXCLUDED.tree_stats",
                &[
                    &(crate::wire::FORMAT_VERSION as i32),
                    &v.tree.world_seed.to_be_bytes().to_vec(),
                    &(v.tree.root.0 as i32),
                    &v.time,
                    &v.pace,
                    &v.time_rate,
                    &v.time_throttle,
                    &(v.paced_to.0 as i64),
                    &(v.pace_mode == crate::engine::PaceMode::Fixed),
                    &v.labour_rate,
                    &(v.rejected_transactions as i64),
                    &stats.finish(),
                ],
            )
            .map_err(db)?;
        Ok(())
    }

    fn write_nodes(&mut self, v: &WorldView<'_>, which: &[usize]) -> Result<()> {
        let mut tx = self.client.transaction().map_err(db)?;
        let stmt = tx
            .prepare(
                "INSERT INTO node (key, idx, parent, depth, tier, alive, pinned, epoch,
                                   instant, last_solved, last_disturbed, payload)
                 VALUES ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12)
                 ON CONFLICT (key) DO UPDATE SET
                    idx = EXCLUDED.idx, parent = EXCLUDED.parent, depth = EXCLUDED.depth,
                    tier = EXCLUDED.tier, alive = EXCLUDED.alive, pinned = EXCLUDED.pinned,
                    epoch = EXCLUDED.epoch, instant = EXCLUDED.instant,
                    last_solved = EXCLUDED.last_solved,
                    last_disturbed = EXCLUDED.last_disturbed, payload = EXCLUDED.payload",
            )
            .map_err(db)?;
        for &i in which {
            let n = &v.tree.nodes[i];
            let mut w = Writer::new();
            crate::persist::put_node_payload(&mut w, n);
            tx.execute(
                &stmt,
                &[
                    &key_bytes(n.key).to_vec(),
                    &(i as i32),
                    &(n.parent.0 as i32),
                    &(n.depth as i32),
                    &(n.tier.index() as i16),
                    &n.alive,
                    &n.pinned,
                    &(n.epoch as i32),
                    &n.time,
                    &n.last_solved,
                    &n.last_disturbed,
                    &w.finish(),
                ],
            )
            .map_err(db)?;
        }
        tx.commit().map_err(db)?;
        Ok(())
    }

    fn write_side_tables(&mut self, v: &WorldView<'_>) -> Result<()> {
        let mut tx = self.client.transaction().map_err(db)?;

        // Small and rarely changing, so rewritten whole. The node table is the
        // one that has to be incremental.
        tx.execute("DELETE FROM pinned", &[]).map_err(db)?;
        let mut keys: Vec<&PathKey> = v.tree.persisted.keys().collect();
        keys.sort_by_key(|k| k.0);
        for k in keys {
            let mut w = Writer::new();
            crate::persist::put_bodies_pub(&mut w, &v.tree.persisted[k]);
            tx.execute(
                "INSERT INTO pinned (key, bodies) VALUES ($1, $2)",
                &[&key_bytes(*k).to_vec(), &w.finish()],
            )
            .map_err(db)?;
        }

        tx.execute("DELETE FROM fact", &[]).map_err(db)?;
        let mut facts: Vec<_> = v.ledger.entries().collect();
        facts.sort_by_key(|(k, q, _)| (k.0, *q as u8));
        for (k, q, f) in facts {
            tx.execute(
                "INSERT INTO fact (key, quantity, value, at, seq) VALUES ($1,$2,$3,$4,$5)",
                &[
                    &key_bytes(k).to_vec(),
                    &(q as i16),
                    &f.value,
                    &f.time,
                    &(f.sequence as i64),
                ],
            )
            .map_err(db)?;
        }
        tx.execute(
            "INSERT INTO ledger_counters (id, sequence, queries, commits) VALUES (1,$1,$2,$3)
             ON CONFLICT (id) DO UPDATE SET sequence = EXCLUDED.sequence,
                queries = EXCLUDED.queries, commits = EXCLUDED.commits",
            &[
                &(v.ledger.sequence as i64),
                &(v.ledger.queries as i64),
                &(v.ledger.commits as i64),
            ],
        )
        .map_err(db)?;

        tx.execute("DELETE FROM environment", &[]).map_err(db)?;
        let mut envs: Vec<_> = v.environments.iter().collect();
        envs.sort_by_key(|(k, _)| k.0);
        for (k, e) in envs {
            let mut w = Writer::new();
            crate::persist::put_environment_pub(&mut w, e);
            tx.execute(
                "INSERT INTO environment (key, data) VALUES ($1,$2)",
                &[&key_bytes(*k).to_vec(), &w.finish()],
            )
            .map_err(db)?;
        }

        tx.execute("DELETE FROM audit", &[]).map_err(db)?;
        for (i, a) in v.audit.iter().enumerate() {
            tx.execute(
                "INSERT INTO audit (seq, key, property, delta_energy, at) VALUES ($1,$2,$3,$4,$5)",
                &[
                    &(i as i64),
                    &key_bytes(a.key).to_vec(),
                    &(a.property as i16),
                    &a.delta_energy,
                    &a.time,
                ],
            )
            .map_err(db)?;
        }

        // In-flight influences. Small, and durable: an impulse in the light
        // delay between the act and its landing is an action somebody took.
        tx.execute("DELETE FROM in_flight", &[]).map_err(db)?;
        let mut flying: Vec<_> = v.mailbox.in_flight().collect();
        flying.sort_by(|a, b| {
            a.arrives
                .partial_cmp(&b.arrives)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.target.0.cmp(&b.target.0))
        });
        for (i, inf) in flying.iter().enumerate() {
            tx.execute(
                "INSERT INTO in_flight (seq, arrives, target, kind, energy, px, py, pz,
                                        source_distance)
                 VALUES ($1,$2,$3,$4,$5,$6,$7,$8,$9)",
                &[
                    &(i as i64),
                    &inf.arrives,
                    &(inf.target.0 as i32),
                    &(crate::persist::influence_kind_tag(inf.kind) as i16),
                    &inf.energy,
                    &inf.momentum.x,
                    &inf.momentum.y,
                    &inf.momentum.z,
                    &inf.source_distance,
                ],
            )
            .map_err(db)?;
        }
        tx.execute(
            "INSERT INTO mailbox_counters (id, delivered, in_flight_peak) VALUES (1,$1,$2)
             ON CONFLICT (id) DO UPDATE SET delivered = EXCLUDED.delivered,
                in_flight_peak = EXCLUDED.in_flight_peak",
            &[&(v.mailbox.delivered as i64), &(v.mailbox.in_flight_peak as i64)],
        )
        .map_err(db)?;

        tx.commit().map_err(db)?;
        Ok(())
    }
}

impl WorldStore for PostgresStore {
    fn save(&mut self, view: WorldView<'_>) -> Result<()> {
        let all: Vec<usize> = (0..view.tree.nodes.len()).collect();
        self.write_world_row(&view)?;
        self.write_nodes(&view, &all)?;
        self.write_side_tables(&view)
    }

    fn save_since(&mut self, view: WorldView<'_>, since: f64) -> Result<Flushed> {
        let which = dirty_nodes(view.tree, since);
        let total = view.tree.nodes.len();
        self.write_world_row(&view)?;
        self.write_nodes(&view, &which)?;
        self.write_side_tables(&view)?;
        Ok(Flushed {
            nodes_written: which.len(),
            nodes_skipped: total - which.len(),
            incremental: true,
        })
    }

    fn load(&mut self) -> Result<Snapshot> {
        load_impl(self)
    }

    fn size(&self) -> usize {
        0
    }
}

fn load_impl(store: &mut PostgresStore) -> Result<Snapshot> {
    let world = store
        .client
        .query_opt("SELECT * FROM world WHERE id = 1", &[])
        .map_err(db)?
        .ok_or_else(|| WireError::Io("no world row".into()))?;

    let format: i32 = world.get("format");
    if format != crate::wire::FORMAT_VERSION as i32 {
        return Err(WireError::UnsupportedVersion {
            found: format as u16,
            supported: crate::wire::FORMAT_VERSION,
        });
    }
    let seed_bytes: Vec<u8> = world.get("world_seed");
    if seed_bytes.len() != 8 {
        return Err(WireError::Truncated { what: "world seed", need: 8, have: seed_bytes.len() });
    }
    let mut sb = [0u8; 8];
    sb.copy_from_slice(&seed_bytes);
    let world_seed = u64::from_be_bytes(sb);
    let root = NodeIdx(world.get::<_, i32>("root") as u32);
    let time: f64 = world.get("instant");
    let stats_blob: Vec<u8> = world.get("tree_stats");
    let mut sr = Reader::new(&stats_blob);
    let stats = crate::persist::get_tree_stats(&mut sr)?;

    // Ordered by arena slot, so `NodeIdx` means the same thing it did.
    let rows = store
        .client
        .query("SELECT idx, payload FROM node ORDER BY idx", &[])
        .map_err(db)?;
    let mut nodes = Vec::with_capacity(rows.len());
    for (expected, row) in rows.iter().enumerate() {
        let idx: i32 = row.get("idx");
        if idx as usize != expected {
            return Err(WireError::Io(format!(
                "node arena has a hole: expected slot {expected}, found {idx}"
            )));
        }
        let payload: Vec<u8> = row.get("payload");
        let mut r = Reader::new(&payload);
        let n = crate::persist::get_node_payload(&mut r)?;
        r.finish()?;
        nodes.push(n);
    }

    let mut persisted = std::collections::HashMap::new();
    for row in store.client.query("SELECT key, bodies FROM pinned", &[]).map_err(db)? {
        let k = key_from(&row.get::<_, Vec<u8>>("key"))?;
        let blob: Vec<u8> = row.get("bodies");
        let mut r = Reader::new(&blob);
        persisted.insert(k, crate::persist::get_bodies_pub(&mut r)?);
    }

    let mut facts = Vec::new();
    for row in store
        .client
        .query("SELECT key, quantity, value, at, seq FROM fact", &[])
        .map_err(db)?
    {
        let k = key_from(&row.get::<_, Vec<u8>>("key"))?;
        let q = crate::persist::quantity_from(row.get::<_, i16>("quantity") as u8)?;
        facts.push((
            k,
            q,
            crate::observe::Fact {
                value: row.get("value"),
                time: row.get("at"),
                quantity: q,
                sequence: row.get::<_, i64>("seq") as u64,
            },
        ));
    }
    let counters = store
        .client
        .query_opt("SELECT sequence, queries, commits FROM ledger_counters WHERE id = 1", &[])
        .map_err(db)?;
    let (seq, queries, commits) = match counters {
        Some(c) => (
            c.get::<_, i64>("sequence") as u64,
            c.get::<_, i64>("queries") as u64,
            c.get::<_, i64>("commits") as u64,
        ),
        None => (0, 0, 0),
    };
    let ledger = crate::observe::Ledger::restore(facts, seq, queries, commits);

    let mut environments = std::collections::HashMap::new();
    for row in store.client.query("SELECT key, data FROM environment", &[]).map_err(db)? {
        let k = key_from(&row.get::<_, Vec<u8>>("key"))?;
        let blob: Vec<u8> = row.get("data");
        let mut r = Reader::new(&blob);
        environments.insert(k, crate::persist::get_environment_pub(&mut r)?);
    }

    let mut audit = Vec::new();
    for row in store
        .client
        .query("SELECT key, property, delta_energy, at FROM audit ORDER BY seq", &[])
        .map_err(db)?
    {
        audit.push(crate::observe::AuthorEvent {
            key: key_from(&row.get::<_, Vec<u8>>("key"))?,
            property: crate::persist::property_from(row.get::<_, i16>("property") as u8)?,
            delta_energy: row.get("delta_energy"),
            time: row.get("at"),
        });
    }

    let mut in_flight = Vec::new();
    for row in store
        .client
        .query(
            "SELECT arrives, target, kind, energy, px, py, pz, source_distance
             FROM in_flight ORDER BY seq",
            &[],
        )
        .map_err(db)?
    {
        in_flight.push(crate::causal::Influence {
            arrives: row.get("arrives"),
            target: NodeIdx(row.get::<_, i32>("target") as u32),
            kind: crate::persist::influence_kind_from(row.get::<_, i16>("kind") as u8)?,
            energy: row.get("energy"),
            momentum: crate::math::v3(row.get("px"), row.get("py"), row.get("pz")),
            source_distance: row.get("source_distance"),
        });
    }
    let (delivered, peak) = match store
        .client
        .query_opt("SELECT delivered, in_flight_peak FROM mailbox_counters WHERE id = 1", &[])
        .map_err(db)?
    {
        Some(c) => (
            c.get::<_, i64>("delivered") as u64,
            c.get::<_, i64>("in_flight_peak") as usize,
        ),
        None => (0, 0),
    };

    let mut snapshot = Snapshot {
        mailbox_view: crate::causal::Mailbox::restore(in_flight.clone(), delivered, peak),
        in_flight,
        delivered,
        in_flight_peak: peak,
        tree: crate::tree::Tree::restore(nodes, root, world_seed, persisted, stats),
        ledger,
        time,
        pace: world.get("pace"),
        time_rate: world.get("time_rate"),
        time_throttle: world.get("time_throttle"),
        paced_to: NodeIdx(world.get::<_, i64>("paced_to") as u32),
            pace_mode: if world.get::<_, bool>("pace_fixed") {
                crate::engine::PaceMode::Fixed
            } else {
                crate::engine::PaceMode::Follow
            },
        labour_rate: world.get("labour_rate"),
        rejected_transactions: world.get::<_, i64>("rejected") as u64,
        environments,
        audit,
    };
    // Nodes skipped by an incremental write are at whatever instant they were
    // last written at. Carrying them forward is exact.
    snapshot.settle();
    Ok(snapshot)
}
