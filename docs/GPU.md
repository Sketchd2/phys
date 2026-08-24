# GPU mapping

The reference implementation is single-threaded, scalar, and `f64`. That is
deliberate: it is the bit-exact oracle the GPU path is validated against. This
document describes the intended compute path, which is **designed but not
written**.

## What must survive the port

The engine's guarantees are not negotiable across backends:

1. **Bit-exact agreement with the CPU reference** for materialisation and
   restriction. Anything else and a scenario recorded on one machine cannot be
   replayed on another.
2. **Deterministic reductions.** No atomics whose ordering varies between runs.
   Sums use the same fixed pairwise tree as `math::det_sum`, which means a
   two-pass reduction with a fixed block shape rather than `atomicAdd`.
3. **Order-independent randomness.** Already satisfied: `Stream::nth_u64(i)` is
   a pure function of the address and index, so lane *i* computes draw *i* with
   no sequential dependency. This is why the RNG is addressed rather than
   streamed.

## Precision

`f64` on a 2060 runs at 1/32 rate — unusable for the hot loops. The split:

| Quantity | Precision | Why |
|---|---|---|
| Positions within a node | `f32` | Node-local coordinates span ≤ 6 orders; `f32` gives 7 digits |
| Velocities, accelerations | `f32` | Same |
| Node aggregates, conserved tuples | `f64` | Where exactness is claimed, and there are only ~10⁶ of them |
| Restriction reductions | `f64` | The conservation guarantee lives here |
| Time | `f64` | 10⁶⁰ dynamic range across tiers |

This split is exactly what the no-global-coordinates design already enables:
because positions are always node-relative, `f32` is sufficient for them, and
the `f64` work is confined to the small aggregate layer.

## Data layout

`Body` is 184 bytes on the CPU, dominated by the 8-species composition. On the
GPU it splits:

```wgsl
// hot: touched every step, 32 bytes
struct BodyHot {
    pos: vec3<f32>,     // 12   node-local metres
    mass: f32,          //  4
    vel: vec3<f32>,     // 12
    radius: f32,        //  4
};
// warm: touched by some solvers
struct BodyWarm { charge: f32, temperature: f32, internal: f32, spin: vec3<f32>, };
// cold: composition[8], kind, slot — parallel arrays, read on materialise/restrict only
```

32 bytes hot raises the resident ceiling from 19.6 M bodies to ~110 M on a 6 GB
card, and cuts the bandwidth per force pass by 5.75×.

## Kernels

| Kernel | Shape | Notes |
|---|---|---|
| `materialise` | 1 thread/body | `prolong`'s samplers are per-body pure functions of `(address, i)`. The projection passes are three reductions plus a broadcast. |
| `restrict` | tree reduction | Deterministic pairwise, fixed block shape, `f64` accumulator |
| `build_octree` | radix sort on Morton codes | Standard: sort, then build the hierarchy from the sorted order (Karras 2012) |
| `traverse_gravity` | 1 thread/body, warp-shared stack | The divergence-critical kernel; see below |
| `sph_density`, `sph_force` | 1 thread/body, cell lists | Regular and cache-friendly; the easy case |
| `md_force` | 1 thread/body, Verlet lists | Same |
| `statistical_step` | 1 thread/body | Trivially parallel; no interactions |

### The traversal kernel is the hard one

Barnes-Hut traversal diverges: each body opens a different set of cells, so
warps do wasted work. The standard mitigations, in order of value:

1. **Group traversal.** Bodies in the same leaf share almost their entire
   interaction list. Computing one list per group of 32 and applying it to all
   32 turns divergent traversal into a coalesced inner loop. This is the single
   biggest win and is what production codes do.
2. **Warp-synchronous stack in shared memory**, one stack per warp rather than
   per thread — the technique from Burtscher & Pingali.
3. **Sort bodies by Morton code** so that neighbouring threads handle
   neighbouring bodies, maximising list overlap.

## Multi-rate scheduling on the GPU

The multi-rate scheme is a natural fit rather than an obstacle. Nodes are
grouped into **rate classes** by their permitted timestep (powers of two apart),
and each class is dispatched as one kernel launch over a compacted body list. A
class that steps 1024× more often than another simply launches 1024× more often
over a much smaller list.

The causal constraint becomes a host-side check between launches: a class may
advance only while every disjoint neighbour class satisfies the light-travel
bound (`Tree::sibling_separations`). Because that bound is generous at large
separations, most classes never synchronise at all.

## Validation plan

The reference implementation is the oracle:

1. Run `tests/consistency.rs` against the GPU backend and require the *same*
   error bounds, not merely small ones.
2. Assert bit-exact agreement on materialisation for a fixed address set. Any
   divergence means the RNG or the reduction order has drifted.
3. Cross-check forces body-by-body against the CPU tree at θ = 0, where
   Barnes-Hut reduces to direct summation and there is a single right answer.
4. Run the full frame loop and compare the conserved tuple trajectory over 10⁴
   frames.

## What is not planned

- **Rendering** is out of scope here. The engine produces `Sighting`s with
  retarded positions, Doppler factors and fluxes; turning those into pixels is a
  separate concern, and deliberately so — the physics must not depend on the
  renderer.
- **Multi-GPU.** The causal decomposition would make it natural (regions with
  large separations are exactly the ones that need no synchronisation), but a
  2060 is the stated target.
