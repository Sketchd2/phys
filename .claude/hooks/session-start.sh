#!/bin/bash
#
# Bring up the PostgreSQL the feature-gated store tests need.
#
# `tests/postgres.rs` skips when `PHYS_PG` is unset and *fails loudly* when it
# is set and does not work — that asymmetry is deliberate, because the tests
# once reported six passes while running nothing at all. So this script's
# contract is: export `PHYS_PG` only after proving a real connection, and
# otherwise leave it unset and let the tests skip.
#
# Idempotent. Safe to run against a container that already has the cluster up.

set -euo pipefail

# Web sessions only. A developer running locally has their own PostgreSQL and
# would not thank us for starting a second one.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

PGDATA=/var/lib/postgresql/physdata
PGPORT=5433
PGUSER=phys
PGDB=phys
PGLOG=/var/lib/postgresql/phys.log

# The Debian package ships a cluster of its own at 5432. It is down and should
# stay down, but sharing a port with something that might come up is how a
# working setup turns into an intermittent one, so ours sits next door.

BIN=""
for d in /usr/lib/postgresql/*/bin; do
  [ -x "$d/pg_ctl" ] && BIN="$d"
done
if [ -z "$BIN" ]; then
  echo "no PostgreSQL server installed; the store tests will skip" >&2
  exit 0
fi

as_postgres() { su postgres -s /bin/bash -c "$1"; }

if [ ! -d "$PGDATA" ]; then
  echo "creating cluster at $PGDATA"
  install -d -o postgres -g postgres "$PGDATA"
  as_postgres "$BIN/initdb -D $PGDATA -A trust -U postgres" >/dev/null
fi

# Put the port in the cluster's own configuration rather than passing it on the
# command line. It used to live only in a flag, so a restart came back on a
# different port and every connection string in the session was wrong.
if ! grep -qE "^port = $PGPORT" "$PGDATA/postgresql.conf" 2>/dev/null; then
  sed -i '/^port = /d' "$PGDATA/postgresql.conf"
  echo "port = $PGPORT" >> "$PGDATA/postgresql.conf"
  echo "listen_addresses = '127.0.0.1'" >> "$PGDATA/postgresql.conf"
  chown postgres:postgres "$PGDATA/postgresql.conf"
fi

if ! "$BIN/pg_isready" -h 127.0.0.1 -p "$PGPORT" -q 2>/dev/null; then
  echo "starting PostgreSQL on port $PGPORT"
  touch "$PGLOG" && chown postgres:postgres "$PGLOG"
  as_postgres "$BIN/pg_ctl -D $PGDATA -l $PGLOG -w -t 30 start" >/dev/null || {
    echo "PostgreSQL would not start; see $PGLOG. The store tests will skip." >&2
    exit 0
  }
fi

# Which superuser this cluster actually has. `initdb -U postgres` above gives
# one called postgres, but a cluster created by hand earlier may have been
# given any name, so ask rather than assume — the first attempt at this hook
# died on a cluster whose superuser was `phys`.
SUPER=""
for candidate in postgres "$PGUSER" root; do
  if psql -h 127.0.0.1 -p "$PGPORT" -U "$candidate" -d postgres -tAc "SELECT 1" >/dev/null 2>&1; then
    SUPER="$candidate"
    break
  fi
done
if [ -z "$SUPER" ]; then
  echo "PostgreSQL is up on $PGPORT but no superuser we know of can connect;" >&2
  echo "leaving PHYS_PG unset so the store tests skip rather than fail." >&2
  exit 0
fi

psql_root() { psql -h 127.0.0.1 -p "$PGPORT" -U "$SUPER" -d postgres -tAc "$1"; }

if [ "$(psql_root "SELECT 1 FROM pg_roles WHERE rolname='$PGUSER'" || true)" != "1" ]; then
  psql_root "CREATE ROLE $PGUSER LOGIN SUPERUSER" >/dev/null
fi
if [ "$(psql_root "SELECT 1 FROM pg_database WHERE datname='$PGDB'" || true)" != "1" ]; then
  psql_root "CREATE DATABASE $PGDB OWNER $PGUSER" >/dev/null
fi

# Only now, with a connection actually proven, is it safe to point the tests at
# it. A `PHYS_PG` that does not work is worse than none: it fails the suite.
CONN="host=127.0.0.1 port=$PGPORT user=$PGUSER dbname=$PGDB"
if psql "$CONN" -tAc "SELECT 1" >/dev/null 2>&1; then
  if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
    echo "export PHYS_PG='$CONN'" >> "$CLAUDE_ENV_FILE"
  fi
  echo "PostgreSQL ready: $CONN"
else
  echo "PostgreSQL is up but '$PGDB' is not reachable; leaving PHYS_PG unset" >&2
fi

# Warm the registry so the first `--features postgres` build is not also the
# first download. Non-fatal: an offline container should still start.
if command -v cargo >/dev/null 2>&1; then
  cargo fetch --quiet 2>/dev/null || true
fi
