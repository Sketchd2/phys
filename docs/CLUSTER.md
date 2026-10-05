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

## What it does not do

It does not split one calculation across machines (the fit table would have to
cross the network each iteration). Memory is matched on what a machine says is
free *now* and what a task has been seen to use: a task of a kind nothing has
finished yet is tried anywhere unless the plan gives `mem=`, and a machine that
is killed for memory anyway is seen only as a lease running out, so the task
is retried elsewhere. The speeds are for the machines as they were while they
were measured: a Pi that starts throttling is not noticed until its results
come back slower. Classes are only `cpu` and `gpu`.
