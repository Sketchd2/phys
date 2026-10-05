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

## Memory

A water pair is two molecules with counterpoise, solved with density fitting;
its tables are several GB (measure it: the log of a worker says how long each
took, and `--status` says who holds what). The 8 GB Pis are the ones to try
first; a 4 GB one may not fit. If a worker is killed for memory the lease runs
out and another machine takes the task, so nothing is lost but the time.
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

A worker does ten tasks per process and then restarts (a long-lived process
slows down; the leak has not been found), waits a minute if the server cannot be
reached and gives up after an hour. A task that fails on three machines is left
alone and reported on the server's console. A result that cannot be delivered
is kept in `unsent-<name>.txt` on the worker. The server's `queue.log` has one
line a result: time, machine, class, file, index, seconds, which is how to tell
what a Pi is worth.

## What it does not do

It does not split one calculation across machines (the fit table would have to
cross the network each iteration), and it does not know how much memory a task
needs or a machine has. A machine that runs out of memory is killed by its
operating system, which the server sees only as a lease running out, so a task
too big for the small machines would be taken and lost by each of them in turn.
Classes are only `cpu` and `gpu`; keeping big tasks (hexamers) off the Pis
needs a third class or a memory field, which is not built and is the first
thing to add if it comes to that.
