# Running the pair work on many machines

`PLAY.md` E8's data is hundreds of independent calculations (a pair, a trimer,
a cluster cut from a liquid), each taking minutes to hours on one machine. They
share nothing, so the work goes to whatever machines are free through a queue:
`src/queue.rs`, with `phys-queue` (server) and `phys-worker` /
`phys-worker-gpu` (workers). Why it is shaped as it is is in that file's doc
comment; this is how to use it.

## What each machine needs

- The same commit, built for it: `cargo build --release --bin phys-worker`
  (and `-p phys-gpu --bin phys-worker-gpu` where there is a GPU). A Pi 4 can
  build it itself (slowly) or take an `aarch64-unknown-linux-gnu` cross build.
  Neither Linux nor ARM has been run yet: the first thing to do on a new
  kind of machine is `cargo test --test queue` and one timed pair.
- In its working directory, the same `grow-<name>.state` as the server (they
  are in the repo) and, for a plan that uses one, the snapshot file. A snapshot
  that is not the server's is refused by fingerprint.
- A route to the server's port (default 7878). There is **no encryption and
  only a shared token** (`PHYS_QUEUE_TOKEN`, default `open`): trusted networks
  only.

## How machines are classified

Each worker says what it is on every request: its class, its cores and the
memory it has free (read from its own operating system: `MemAvailable` on
Linux, `GlobalMemoryStatusEx` on Windows; `--mem GB` offers less, or says how
much where it cannot be read). The server learns the rest by watching:

- **Speed** is measured, never declared: seconds per unit of task weight, per
  machine, from the results it has returned. A machine nobody has measured is
  treated as middling.
- **Memory a task needs** is what the plan says (`mem=GB`) or, failing that,
  the most any finished task of its kind has used, which every worker reports
  with its result (peak working set, `VmHWM` on Linux). Until a task of a kind
  has finished nothing is known and any machine may try it, so give the first
  tasks of a new kind a `mem=` before they go to the Pis.

From those: a task goes only to a machine it fits in; a machine that has only
too-big tasks left is told so and stops (exit 4) rather than asking again;
among tasks that fit, the heaviest go to the fastest machines and the lightest
to the slowest (equal weights go in order, so a plan of pairs is untouched);
and a lease is as long as that machine's own speed says the task needs (four
times its measured seconds for the task's weight, never less than the base
lease). `--status` prints each machine with its cores, memory, count done and
measured speed, and the most memory seen for each kind.

A task's weight is the work in it relative to a pair (`weight=` in the plan,
default 1). It is a judgement until measured: a hexamer is not 20 pairs
because someone said so, so check it against the speeds `--status` reports
for it and correct the plan.

## Over Tailscale

The protocol has no encryption and only a shared token (see `src/queue.rs`), so
it belongs on a private network. A tailnet is one: WireGuard encrypts the
traffic and only enrolled machines can reach the server.

- **Bind the server to the tailnet address**, not to everything:
  `phys-queue plan.txt --listen 100.x.y.z:7878` (the machine's Tailscale IP,
  from `tailscale ip -4`). The default, `0.0.0.0`, also answers on the home LAN.
- **Windows firewall:** allow inbound TCP 7878 for the Tailscale interface (or
  for the `100.64.0.0/10` range), and nothing wider.
- **Workers name the server by its MagicDNS name or its tailnet IP:**
  `scripts/worker-loop.sh home-pc:7878 pi-017`.
- **Keep the token** (`PHYS_QUEUE_TOKEN`, the same on both ends) as a second
  layer, and if the tailnet has ACLs let the Pis reach port 7878 on the server
  and nothing else of it.
- **Data caps:** the queue's own traffic is under 1 MB a day for a Pi, and
  WireGuard adds a few dozen bytes a packet. Tailscale's own keepalives and path
  discovery cost something too, which has **not been measured**: read one Pi's
  counters over a day before enrolling the fleet. A Pi that cannot make a direct
  connection goes through a relay, which carries the same bytes by a longer route.
- **Installing Tailscale on each Pi is a separate download** from anything here;
  do not count it in the worker's figure.

Per Pi, one-time: the ARM binary (cross-compiled once and copied; building on
each Pi would pull the toolchain and every crate, hundreds of MB), `grow-water.state`
(2.5 KB) and, for a plan that names one, the snapshot (25 KB). A result is one
268-byte line, about 1.5-2 KB on the wire with its connection (an estimate;
the message sizes are measured). A worker with nothing to take asks again every
60 s, about 1 MB a day if it sits idle.

## Memory

A water pair is two molecules with counterpoise, solved with density fitting.
What it takes is measured now: a worker reports its peak, and `--status` shows
the largest per kind (the long-lived pair processes on the desktop sat at
1.8-2.5 GB, with the leak in it). The 8 GB Pis are the ones to try first and a
4 GB one may fit pairs. If a worker is killed for memory the lease runs out
and another machine takes the task, so nothing is lost but the time.
`PHYS_SPILL_DIR` spills the tables to disk, which on an SD card is both slow and
hard on the card; leave it unset on a Pi.

## The plan

`plan.txt`, one range per line (`#` comments):

```text
pair water 300 600 element cpu pairs-water-cpu.txt
snap water bulk-water-298.snap 0 40 element cpu pairs-water-liq-cpu.txt
```

`pair` is indices `from..to` of the random draw, `snap` the same of a
snapshot's pairs; `element` or `grown` is the basis; `cpu` or `gpu` the class of
machine that may do it, and the file the lines go to. **Class matters to the
data**: the GPU's final non-local energy is single precision, the CPU's double,
and the queue never lets one machine's results stand in for the other's. Which
files a fit may combine is the owner's call. An index already in the file is
left out, so a plan is safe to run again; indices in a `cpu` file that start
above the `gpu` file's are new pairs, while the same indices would be the same
geometries computed the other way (useful for the precision test).

## Running

```sh
phys-queue plan.txt                 # on the server; --listen, --lease MINUTES, --log
phys-queue --status server:7878     # who has done what, who holds what
scripts/worker-loop.sh server:7878 pi-017      # on each worker; keeps it going
```

## Checking the data

`phys-recheck file molecule [--sample N] [--seed S] [--indices a,b,c]`
(`phys-recheck-gpu` for a GPU file) reads a pairs file back. It always checks
the file against itself, with no solve: every line whole, no index twice or
missing, each molecule rigid, each atom's type, the separation written the
distance it says. With `--sample` or `--indices` it then solves the chosen pairs
again from the geometry recorded and prints the energies beside the recorded
ones; `--tol` (default 0.01 kcal/mol) only sets the exit status. The geometry
has 8 decimals, so a repeat differs by that rounding as well as by anything
that changed.

A worker does ten tasks per process and then restarts (a long-lived process
slows down; the leak has not been found), waits a minute if the server cannot be
reached and gives up after an hour; it does not restart after exit 4 (nothing that is left fits it). A task that fails on three machines is left
alone and reported on the server's console. A result that cannot be delivered
is kept in `unsent-<name>.txt` on the worker. The server's `queue.log` has one
line a result: time, machine, class, file, index, seconds, which is how to tell
what a Pi is worth.

## Clusters (the many-body test)

`phys-cluster water bulk-water-298.snap --size N --centre I --law law.txt`
(`phys-cluster-gpu` for the GPU's final non-local energy, which is how the pair
data was made and so what a cluster must be compared with) cuts `N` molecules
about molecule `I` of a snapshot, computes the cluster's counterpoise
interaction energy and that of each of its pairs, and says how far the cluster
is from the sum of its pairs and from a law's sum. `--dry` prices a cluster
without solving it. Measured on the desktop (RTX 2060, 16 GB, spill to H:):

```text
molecules  functions  table, at most   first field   later fields
2          428        1.4 GB           72 s          47 s   (table reused)
3          642        4.6 GB           195 s         131-145 s
6          1284       37 GB            4648 s        not yet measured
```

The six-molecule field is disk-bound (13 iterations of 314 s, three passes over
the 30 GB spilled table at 0.28 GB/s); with the table in memory it is not, so
**the run wants a machine with 48 GB or more free**, and the share of free RAM
the table may use is half by default, which on a 64 GB machine still spills:
raise it with `PHYS_TABLE_RAM_FRACTION` (below).

**Interrupted?** Run the same command again. Every pair and every field is on its
line in `cluster-<name>-<N>-c<I>.txt` the moment it is finished, so a run
carries on from the file and loses at most the field it was in. Tested by killing
a trimer after two of its four fields: the restart read the pairs and both
fields from the file, solved the last two, and gave the same answer. A file
begun on a different cluster (another snapshot under the same name) is refused,
and the spill a killed run left behind is deleted at the start (only files
`phys-fit-<pid>-<n>.bin` whose process is known not to exist). It cannot resume
inside a field: an SCF is one process.

### The share of memory the table may keep

`PHYS_TABLE_RAM_FRACTION` (0.05 to 0.9, default 0.5 of the memory free when the
table is built) is how much of the three-centre table stays in RAM; the rest goes
to the spill disk or is rebuilt. It changes speed and nothing else: a pair's
field with 5% of free RAM (spilled) and with the default gave the same energy
bit for bit. For the six-molecule cluster (a table of at most 37 GB) set it so
the whole table fits, e.g. `PHYS_TABLE_RAM_FRACTION=0.8` on a machine with
64 GB free and nothing else running; the solve also holds the grid's basis
values, the densities and the non-local kernel, so do not leave it nothing.
An out-of-range value is ignored with a message. Check what a machine really
has free before trusting a fraction of it.

### Across machines

After the table is built, a cluster's fields and pairs are independent, so
they can go to different machines and be put together:

```sh
# the six-molecule cluster: fields 0 (the cluster) to 6 (molecule 5), pairs I,J
phys-cluster water bulk-water-298.snap --size 6 --centre 0 --field 3 --out field3.txt
phys-cluster water bulk-water-298.snap --size 6 --centre 0 --pair 1,4 --out pair-1-4.txt
# then, on any one machine with all the files copied to it
phys-cluster water bulk-water-298.snap --size 6 --centre 0 --law law.txt     --merge field0.txt field1.txt ... pair-0-1.txt ...
```

`--field K` and `--pair I,J` do one piece, write its line to `--out FILE` (by
default the cluster's own) and stop. `--merge` adds to the main file every
result the others hold that it lacks (never replacing one it has, and refusing
a file begun on a different cluster), and then the run goes on as usual, taking
what is in the file and solving anything still missing. Each machine builds its
own table, so each needs the memory for one: **seven fields on seven machines
is seven 37 GB tables**, not one. A field is the long part (a pair of molecules
takes 1-3 min; a hexamer field took 77 min on this desktop with a disk-bound
table); the 15 pairs are small enough for any machine, the rack servers or the
Pis. Results from machines with different GPUs differ in the last digits; the
GPU's single precision is the same *class* of calculation as the pair data, and
the file does not record which machine made which line, so keep the file names.

Checked on a cluster of two molecules, whose cluster energy is its one pair's:
two fields run separately (`--field 0`, `--field 1`), merged, the rest solved,
gave a cluster energy equal to the pair's, and the non-additive part +0.0000.

## What it does not do

It does not split one calculation across machines (the fit table would have to
cross the network each iteration). Memory is matched on what a machine says is
free *now* and what a task has been seen to use: a task of a kind nothing has
finished yet is tried anywhere unless the plan gives `mem=`, and a machine that
is killed for memory anyway is seen only as a lease running out, so the task
is retried elsewhere. The speeds are for the machines as they were while they
were measured: a Pi that starts throttling is not noticed until its results
come back slower. Classes are only `cpu` and `gpu`.
