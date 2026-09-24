# Watching a test

The suite can assert that a thing happened. It cannot show it. This was the plan
for the second one; **it is built**, as Phase 4 item 2, and the text below is
kept as the design it was built to. What changed in the building is recorded at
the end.

What it is for is the class of question a number cannot answer: the box test
says momentum is conserved to 8×10⁻⁷ and angular momentum to 5×10⁻⁷ over eight
collisions, and says nothing about whether the ball bounced around inside the
box or tunnelled through a panel and came back. Both conserve. One is the
scenario and one is a bug, and the cheapest instrument that tells them apart is
a picture.

The same argument is already written into `src/render.rs`'s own module comment —
"so that 'the engine works' can be looked at rather than only asserted" — and
`examples/damage.rs` is the one place it has been taken up.

---

## What exists

**`src/render.rs`**, about four hundred lines and no dependencies:

- `Canvas::new(w, h)` with a depth buffer, `Canvas::sky(top, bottom)`.
- `Camera::framing(centre, radius, azimuth, elevation)` — pull back far enough
  to see all of something, from a corner so the depth reads.
- `Style`, with `daylight()` and `burned()`, and `show_litter`.
- `draw_structure(canvas, camera, bodies, topo, intact, style)` — members as
  shaded tubes, back to front, and loose parts as discs if asked.
- `write_png(canvas, path)` — stored deflate, forty lines instead of a
  dependency, valid PNG, runs anywhere the tests run.

**`examples/damage.rs`** grows a tree for ninety years, loads it to failure four
ways and shoots each. It also carries the one lesson already paid for: *one
camera for every frame in a set.* Auto-framing each shot rescales the subject
and hides exactly what the comparison exists to show — a tree that had lost a
third of its height looked identical to one that had not.

**`tests/topology.rs::renderer_draws_the_structure`** is the only test that
draws. It asserts that coverage is between 2% and 90% of the frame, writes a PNG
to the temp directory, checks the signature and `IEND`, and deletes it. That is
a test *of the renderer*, not a test *shown by* it.

**`phys-headless serve | watch`** is the other half of the picture and a
different program: `serve` writes a `view::Scene` per frame to a byte stream and
`watch` reads it back with no `World`, `Tree` or solver in its call graph. A
`Scene` already solves multi-node composition — each `NodeView` carries an
`offset` and a `scale` in the frame node's radii, and each `Speck` carries
position, velocity, mass, temperature and kind. What a `Scene` does *not* carry
is topology: no joints, no members. A client draws points.

---

## What is missing, measured

**1. A node with no topology draws nothing at all.** `draw_structure` iterates
`bodies.len().min(topo.base.len())`, so a node whose contents are loose produces
an empty sky — `show_litter` included, because the litter branch is inside that
same loop. Everything Phase 1 built is in this class: the ball and the
ninety-six panels of `a_ball_loose_in_a_box`, a node's hydro parcels, an
exchanging pair, a fragment in flight. The one drawable thing in the engine today is a
built structure standing still.

**2. One frame only.** Every position handed to the renderer is in one node's
own coordinates. Two promoted children are two frames, and the offset between
them is `Tree::separation`, which `a_branch_lands_on_the_next_tree` already
uses. Two caveats belong in the plan rather than in a surprise later:
*orientation is not composed anywhere* (backlog: "Derived gravity is in the
parent's axes"), so a composed scene is right in position and unrotated; and
`separation` returns a `Bounded`, whose `err` is the honest statement of when a
scene is below its own precision — the trap that once placed two children 50 m
apart at 2×10²⁰ m and drew them on top of each other.

**3. Fog is 60 metres, always.** `Style::fog` is a constant in the default, and
the depth cue is `z / fog`. At six metres everything is unfogged; at a kilometre
everything is fog. It wants to be a multiple of the framing radius.

**4. Colour carries three labels and no measurement.** `member`, `broken`,
`litter` — which is a label the caller supplies, not a property read off the
thing. Nothing can be coloured by temperature, speed or phase, so the exchange
result (a hot node beside a cold one equilibrating) has nothing to show even
once its bodies can be drawn at all. This is the axiom applied to the renderer:
**a diagnostic image should be a measurement, not a legend.** `Speck` already
carries temperature and speed for exactly this reason on the client side.

**5. Nothing writes a sequence.** One canvas, one file, one call. A scenario
unfolding is a sequence, and every caller that wants one writes its own loop,
its own naming and its own fixed camera — `damage.rs` being the only caller, and
having written all three.

---

## The design

Five pieces, in the order they unblock each other.

### `render::draw`, with the topology optional

`draw(canvas, camera, bodies, Option<&Topology>, &Paint, &Style)`. Members as
tubes where a topology says there are members; everything else as a disc at its
own radius. `draw_structure` becomes the `Some` case and keeps its name and its
signature, so `examples/damage.rs` and the existing test do not move.

This is the whole of step one and it is what makes the box, the fragments and
the parcels visible.

### `render::Paint`, a measured channel

```text
    Paint::Role                      what exists now: member / broken / litter
    Paint::Temperature { lo, hi }    kelvin, ramp
    Paint::Speed { lo, hi }          m/s
    Paint::Phase                     solid / liquid / gas fractions
```

`lo`/`hi` default to the frame's own extremes, printed alongside the image so a
colour can be read back as a number. Nothing here is a category the caller
asserts; each is read off the body being drawn.

### `render::Shot` and `render::Film`

`Shot` is a camera, a style, a paint and a size, made once. `Film::open(name)`
takes a directory and a counter and writes `name-0000.png` upward. One camera
for the whole film, by construction rather than by remembering — which is the
`damage.rs` lesson turned into a type.

### `film::of_node(&World, node, &Shot)` — a separate module

`render.rs` knows `math`, `state` and `topology`, and it should keep knowing
nothing else: it is a rasteriser, and a rasteriser that can reach a `World` will
eventually solve something. So the gathering — refine the node, take its
topology, walk its promoted children, place each one through
`Tree::separation`, flatten to one body list — belongs in a thin module beside
it that *is* allowed to know `World`.

It lives in the library rather than in `tests/`, because every integration test
is its own crate and a helper in one is invisible to the rest.

### Off unless asked

`PHYS_FILM=<dir>` or nothing happens. `cargo test` stays silent, fast and
identical; with the variable set, every instrumented test drops a filmstrip in
that directory. The directory is git-ignored.

---

## What this is deliberately not

**Not an assertion.** No golden images, no pixel comparison, no coverage
threshold beyond the one that already guards the renderer itself. An image
assertion is the archetype of a test that cannot fail — it passes for a hundred
reasons unrelated to what it claims — and this repo's rule is that a test is
verified against the defect it catches. These are diagnostics. The assertions
stay numeric.

**Not a client.** `view::Scene` is what a client gets and it carries no joints;
a viewer built on it is the real-time path in `docs/GPU.md`, and a `watch` that
draws instead of tabulating is a genuinely different program. The test viewer
reads engine state directly and on purpose: a test wants ground truth, not what
a client would have been sent.

**Not interactive.** No window, no input, no frame rate. PNGs in a directory.

---

## First consumers

Named, because a tool with no caller rots:

| test | what only a picture shows |
|---|---|
| `a_ball_loose_in_a_box` | the ball bouncing rather than passing through |
| `a_branch_lands_on_the_next_tree` | which tree it hit, and where |
| `a_hot_node_beside_a_cold_one_equilibrates` | the gradient, over frames |
| `examples/damage.rs` | already the shape of this; the regression check |
| Phase 4's terrain | patch boundaries, and whether they line up |
| `tests/accretion.rs` | **the first one taken up**: a ball of rock becoming a planet |
| Phase 5's beach | the one scenario nobody will believe a number about |

---

## Cost

Measured, and it is the reason `NodeFilm` defaults to 320×240 and takes an
`every`: `write_png` stores rather than compresses, so a frame costs its pixels.
A 200-frame film of the accretion scenario is **100 MB at 480×360 and 46 MB at
320×240**, so a thousand frames of anything is a quarter of a gigabyte and a
long film wants `every`.

---

## What changed in the building

Three things, each because a measurement said so.

**`film::of_node` does not refine**, where the design said "refine the node,
take its topology, walk its promoted children". Materialising sets
`last_disturbed`, spends the frame's byte budget and changes which nodes the
scheduler then finds materialised — a diagnostic that changes the LOD of the
thing it is drawing is measuring the picture rather than the world. So it draws
what is *there*: a node holding bodies draws its bodies, and a node holding none
draws as one disc at its own radius. That second case is what makes a ball of
matter becoming a planet watchable from its first frame, before anything has
resolved it.

**A promoted child replaces its stand-in rather than being drawn beside it.**
`children` runs parallel to `bodies`, so a resolved thing was otherwise drawn
twice: once where it is, and once as the smear it left behind.

**Orientation is composed**, by `Tree::axes_from`. A child's contents are in the
child's own axes, and before Phase 4 nothing in the tree composed a rotation
across more than one level — which is the same gap Phase 4 item 1 closed for
gravity, met again by the instrument.

And one measured defect the work found rather than designed around:
`draw_structure` bounded its loop by the topology's length, so **a node whose
contents are loose drew an empty sky** — which is everything Phase 1 built. A
cloud of 256 bodies now covers 0.02 to 0.90 of the frame where it used to cover
0.000.
