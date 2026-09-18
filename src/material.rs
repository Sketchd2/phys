//! What a thing is made of, measured from what it *is*.
//!
//! # Why this is its own module now
//!
//! `Material` lived in `topology.rs`, beside the joints that read it, as "a
//! struct anyone can construct" holding thirteen numbers. That was right while
//! a material was a property of a *structure*: a `Topology` carried one and
//! nothing else had any. `docs/PLAY.md` D13 retires that — a rock needs a
//! material as much as a wall does, and §2A measured what reading one off a
//! birth certificate costs: `surface_of` returned `None` for a boulder, so a
//! boulder could not collide.
//!
//! So a material is now **measured**, from the substances the matter is
//! actually made of, which D17 put on `Matter` for exactly this. The presets
//! survive as authoring conveniences and as the calibration this file is held
//! to; nothing in the physics looks a material up any more.
//!
//! # What is derived, and from what
//!
//! Everything here comes from `chem::Properties`, which itself comes from an
//! `Arrangement` — atoms and bonds — through laws the engine already had. The
//! chain is short and worth stating, because the point of D14 is that it
//! contains **one stored number and it is a rule, not a value**:
//!
//! | property | from |
//! |---|---|
//! | density | the substance's own, mass-weighted over the mixture |
//! | stiffness `E` | cohesive energy density |
//! | surface energy `gamma` | cohesive energy over the area it is spread on |
//! | strength | Griffith, `sqrt(2 E gamma / pi a)` |
//! | flaw scale `a` | classical nucleation theory, from formation conditions |
//! | specific heat | Dulong and Petit |
//! | thermal limits | the melting point |
//! | destruction enthalpy | cohesive energy per kg |
//! | ductility | ionicity — a bond that transfers charge cannot slip |
//!
//! Two are **not** derived and say so: `resistivity` and `combustible`. There
//! is no law for either in this engine — resistivity needs a band structure and
//! combustibility needs an oxidiser and a kinetics — and inventing one to avoid
//! a gap would be worse than the gap. They are carried as stated properties of
//! the presets and default to the insulating, non-combustible case.

use crate::chem::arrange::{Arrangement, Bond, Lattice, Order};
use crate::chem::{analyse, Element, Mixture, Phase, Properties, Registry};
use crate::units::{H_PLANCK, K_B, N_AVOGADRO};

/// Material properties, as data.
///
/// A closed enum of four materials was enough to get a tree to break
/// convincingly and is exactly the wrong shape for a solver: adding a material
/// meant editing six `match` arms inside the physics. These are numbers, and
/// they belong in a struct that anyone can construct.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    pub name: &'static str,
    /// Bulk density, kg/m^3.
    pub density: f64,
    /// Surface energy, J/m^2: what it costs to make one square metre of new
    /// surface by pulling the material apart.
    ///
    /// Derived from the cohesive energy over the area it is spread on — see
    /// [`Material::of`] — and one of the two numbers Griffith's law needs.
    pub surface_energy: f64,
    /// Characteristic flaw scale, metres. The other number Griffith needs, and
    /// **the only stored number in the whole strength derivation**.
    ///
    /// `docs/PLAY.md` D14 is emphatic about what this is and is not. It is not
    /// a stress somebody chose; it is a *length with physical meaning*, set by
    /// how well the thing was made, and it is derived per object from that
    /// object's own formation conditions by [`grain_scale`] rather than looked
    /// up. A table of eight lengths would be the table of eight stresses it
    /// replaces, wearing different units.
    ///
    /// It pays for itself immediately: strength becomes **size-dependent**,
    /// because a small piece cannot contain a large flaw. That is real,
    /// measurable, and was inexpressible while strength was a constant.
    pub flaw_size: f64,
    /// Tensile strength as a fraction of [`Material::strength`]. Masonry's
    /// asymmetry is why walls topple rather than snap.
    pub tensile_ratio: f64,
    /// Young's modulus, Pa. Sets deflection and how redundant structures share
    /// load between alternative paths.
    pub stiffness: f64,
    /// Temperature at which strength begins to fall, K.
    pub thermal_onset: f64,
    /// Temperature at which no strength remains, K.
    pub thermal_gone: f64,
    /// Enthalpy needed to destroy a kilogram outright — boiling the water in it
    /// and pyrolysing the rest, J/kg.
    pub destruction_enthalpy: f64,
    /// Specific heat, J/kg/K.
    pub specific_heat: f64,
    /// Electrical resistivity, ohm-metres. Sets how a conducted discharge
    /// distributes its energy between members.
    pub resistivity: f64,
    /// Whether the material is consumed rather than merely weakened when it
    /// passes `thermal_gone`.
    pub combustible: bool,
    /// How much of the bonding is metallic, 0 to 1.
    ///
    /// The one thing that decides **which failure law applies**: a metal yields
    /// by moving dislocations and everything else runs a crack. It is not the
    /// same as [`Material::ductility`], and conflating them cost an ordering:
    /// a polymer is tough because its chains pull out, which lets it keep
    /// carrying load past yield, but its *strength* is still Griffith's. Using
    /// ductility for both made green wood stronger than seasoned timber,
    /// because the dislocation law scales with stiffness and green wood is
    /// denser.
    pub metallic: f64,
    /// Yield strength as a fraction of [`Material::strength`], or zero for a
    /// brittle material that fractures instead of yielding.
    ///
    /// This single number is the difference between a steel frame that sags,
    /// redistributes and warns you, and a masonry wall that is standing one
    /// moment and rubble the next.
    pub ductility: f64,
}

impl Default for Material {
    fn default() -> Self {
        Material::green_wood()
    }
}

impl Material {
    /// Every named preset, in a stable order.
    ///
    /// Exists so a material's *label* can survive a save: nothing in the
    /// physics reads `name`, and a `&'static str` cannot be rebuilt from bytes,
    /// so a loaded material recovers its name by matching the list and falls
    /// back to "custom".
    ///
    /// **Derived on first use rather than tabulated.** Each entry is a
    /// substance and a formation history run through [`Material::of`]; none of
    /// the numbers below is typed in. Two pairs in the list are the *same
    /// substance* at different histories — masonry and bedrock are both
    /// silicate, green wood and dry timber are both cellulose — which is
    /// `docs/PLAY.md` D14's own point about what the retired `rupture` table
    /// was secretly encoding.
    pub fn presets() -> &'static [Material] {
        static P: std::sync::OnceLock<Vec<Material>> = std::sync::OnceLock::new();
        P.get_or_init(|| {
            vec![
                Material::green_wood(),
                Material::dry_timber(),
                Material::aragonite(),
                Material::reinforced_frame(),
                Material::masonry(),
                Material::steel(),
                Material::ice(),
                Material::bedrock(),
            ]
        })
    }

    /// The static label matching `name`, or "custom" for a material that was
    /// built rather than chosen.
    pub fn static_name(name: &str) -> &'static str {
        for m in Material::presets() {
            if m.name == name {
                return m.name;
            }
        }
        "custom"
    }

    /// Living wood, wet. Cellulose laid down at ambient temperature, slowly.
    ///
    /// Its difference from [`Material::dry_timber`] is **water**, which is a
    /// property of the mixture rather than of the material, and a slower
    /// formation. Nothing here says "wood is strong in bending"; that falls out
    /// of a low stiffness and a flaw scale set by how slowly a tree lays down a
    /// ring.
    pub fn green_wood() -> Material {
        static M: std::sync::OnceLock<Material> = std::sync::OnceLock::new();
        *M.get_or_init(|| Material {
            name: "green wood",
            combustible: true,
            resistivity: 1.0e4,
            // A tree lays down *cells*, not rings: the ring is a year's bundle
            // of them and its boundary is well bonded, while a cell wall is a
            // discontinuity thirty microns across. Fast summer growth makes
            // wide cells, and wide cells are the flaw.
            ..Material::of(&substances::cellulose(), Formation::laid_down(291.0, 3.0e-5, 600.0 / 1612.0))
        })
    }

    /// Seasoned timber: the same cellulose, drier and slower-grown.
    ///
    /// Everything separating it from [`Material::green_wood`] is in its
    /// formation — finer cells from slower growth, and less of the volume
    /// filled because the water has gone. The old table recorded the pair as
    /// two materials with two stresses; here it is one substance with two
    /// histories, and **the ordering between them falls out** rather than being
    /// stated.
    pub fn dry_timber() -> Material {
        static M: std::sync::OnceLock<Material> = std::sync::OnceLock::new();
        *M.get_or_init(|| Material {
            name: "dry timber",
            combustible: true,
            resistivity: 1.0e8,
            ..Material::of(&substances::cellulose(), Formation::laid_down(291.0, 1.2e-5, 480.0 / 1612.0))
        })
    }

    /// Coral skeleton: calcium carbonate, laid down in sea water.
    pub fn aragonite() -> Material {
        static M: std::sync::OnceLock<Material> = std::sync::OnceLock::new();
        *M.get_or_init(|| Material {
            name: "aragonite",
            resistivity: 1.0e6,
            ..Material::of(&substances::calcium_carbonate(), Formation::laid_down(300.0, 3.0e-3, 2700.0 / 2924.0))
        })
    }

    /// Reinforced concrete and steel: a silicate matrix with iron through it,
    /// set at ambient temperature over weeks.
    pub fn reinforced_frame() -> Material {
        static M: std::sync::OnceLock<Material> = std::sync::OnceLock::new();
        *M.get_or_init(|| Material {
            name: "reinforced frame",
            resistivity: 1.0e-6,
            ..Material::blend(
                &[
                    (substances::silica(), Formation::laid_down(300.0, 1.0e-4, 2400.0 / 2644.0), 0.85),
                    (substances::iron(), Formation::cooled(0.5, 4.0e-4), 0.15),
                ],
            )
        })
    }

    /// Mortared brick or stone: silicate, fired and cooled in hours.
    ///
    /// The *same substance* as [`Material::bedrock`], and everything that
    /// separates them is in the second argument.
    pub fn masonry() -> Material {
        static M: std::sync::OnceLock<Material> = std::sync::OnceLock::new();
        *M.get_or_init(|| Material {
            name: "masonry",
            resistivity: 1.0e9,
            ..Material::of(&substances::silica(), Formation::laid_down(1300.0, 7.0e-2, 1900.0 / 2644.0))
        })
    }

    /// Silicate bedrock, crystallised in the crust over a geological age.
    pub fn bedrock() -> Material {
        static M: std::sync::OnceLock<Material> = std::sync::OnceLock::new();
        *M.get_or_init(|| Material {
            name: "bedrock",
            resistivity: 1.0e10,
            ..Material::of(&substances::silica(), Formation::cooled(1.0e-11, 1.0e-13))
        })
    }

    /// Structural steel: iron, cast and then worked.
    pub fn steel() -> Material {
        static M: std::sync::OnceLock<Material> = std::sync::OnceLock::new();
        *M.get_or_init(|| Material {
            name: "steel",
            resistivity: 1.4e-7,
            ..Material::of(&substances::iron(), Formation::cooled(0.5, 4.0e-4))
        })
    }

    /// Ice, which is a structural material wherever it is cold enough.
    pub fn ice() -> Material {
        static M: std::sync::OnceLock<Material> = std::sync::OnceLock::new();
        *M.get_or_init(|| Material {
            name: "ice",
            resistivity: 1.0e5,
            ..Material::of(&substances::water(), Formation::cooled(1.0e-4, 1.0e-6))
        })
    }

    /// Measure a material from one substance and how it was made.
    ///
    /// `docs/PLAY.md` D13: a thing's material is derived from what it *is*, not
    /// from what generated it. Every column below either has a law or says it
    /// does not.
    ///
    /// | column | law |
    /// |---|---|
    /// | `density` | the substance's own |
    /// | `stiffness` | cohesive energy density: `E ~ 3 U / Omega` |
    /// | `surface_energy` | cleaving breaks about half an atom's bonds and makes two surfaces |
    /// | `flaw_size` | [`grain_scale`] |
    /// | `specific_heat` | Dulong and Petit, `3R` per mole of atoms |
    /// | `thermal_onset` | `0.4 T_m`, where creep begins |
    /// | `thermal_gone` | `T_m` |
    /// | `destruction_enthalpy` | the cohesive energy per kilogram |
    /// | `ductility` | electronegativity: shared electrons let planes slip |
    /// | `tensile_ratio` | ductility again — a brittle solid closes its cracks in compression |
    /// | `resistivity`, `combustible` | **no law.** Defaults, overridden by a caller that knows |
    ///
    /// The two at the bottom are the honest gaps. Resistivity needs a band
    /// structure and combustibility needs an oxidiser and a kinetics, and
    /// inventing either to avoid an admission would be worse than the
    /// admission.
    pub fn of(props: &Properties, formation: Formation) -> Material {
        let atoms = props.atoms_per_unit.max(1) as f64;
        let density = props.density.max(1e-6);
        let volume_atom = (props.unit_mass / (density * atoms)).max(1e-45);
        let energy_atom = (props.cohesive_energy / atoms).max(0.0);

        // Stiffness from the cohesive energy density. A pair potential's bulk
        // modulus is its energy density times a factor of order one set by the
        // curvature of the well at its minimum; three reproduces iron's
        // 211 GPa from 58 GPa of energy density, and is the same three for
        // everything else here.
        // Gibson and Ashby: a cellular solid bends its cell walls where a
        // dense one stretches its bonds, so stiffness falls as the square of
        // how much of the volume is filled and density falls linearly.
        let packing = formation.packing();
        let stiffness = 3.0 * energy_atom / volume_atom * packing * packing;

        // Making a surface means breaking the bonds that cross it. An atom
        // shares its cohesive energy with its neighbours, a cleavage plane cuts
        // about half of those, and the cut makes *two* faces — so a quarter of
        // the energy per atom, spread over the area one atom presents.
        let area_atom = volume_atom.powf(2.0 / 3.0);
        let surface_energy = energy_atom / (4.0 * area_atom);

        let melting = props.melting_point.max(1.0);
        // Dulong and Petit: three degrees of freedom per atom, each worth
        // `k_B`. It over-predicts for light atoms well below their Debye
        // temperature — ice comes out about twice what it should — and is right
        // to within tens of per cent for everything heavier.
        let specific_heat = 3.0 * K_B * atoms / props.unit_mass.max(1e-30);

        // Ductility is the ability to keep carrying load past the elastic
        // limit, and there are **two** ways a solid does that.
        //
        // A **metal** moves dislocations, and what makes that possible is that
        // its electrons belong to everybody — low electronegativity across the
        // board is what "metallic" means.
        let metallic = if props.electronegativity > 0.0 {
            ((2.5 - props.electronegativity) / 1.0).clamp(0.0, 1.0)
        } else {
            0.0
        };
        // A **polymer** pulls its chains out of their neighbours. That needs a
        // long molecule, and a solid of small molecules has nothing to pull:
        // cellulose's repeating unit is twenty-one atoms and ice's is three,
        // which is the whole of why a branch bends a long way before it goes
        // and an icicle does not.
        //
        // Deriving ductility from metallicity alone, which is what the first
        // version did, made green wood perfectly brittle — and a brittle tree
        // does not shed a limb, it shatters. Measured: it went from carrying a
        // 53 m/s gale untouched to losing 289 joints at 54.
        let polymeric = if props.crystalline {
            0.0
        } else {
            // Six atoms is about where a molecule stops being a point and
            // starts being a chain. The 0.6 is how much of the load a pulled
            // fibre carries against a slipping plane, which is the same kind of
            // fraction as the 0.45 in Turnbull's relation.
            (1.0 - 6.0 / atoms).clamp(0.0, 1.0) * 0.6
        };
        let ductility = metallic.max(polymeric);

        // A solid is weak in tension when it has a plane to part along and
        // close again: an ionic or covalent *lattice* cleaves, and that is the
        // whole of masonry's asymmetry. A molecular or fibrous solid has no
        // such plane — which is why wood carries tension as well as it carries
        // compression — and a metal closes the crack by flowing instead.
        //
        // Tying this to metallicity alone, which is what the first version did,
        // made green wood 4% as strong in tension as in compression and snapped
        // a forty-year tree in an 18 m/s gust.
        let cleaves = props.crystalline && ductility < 0.5;
        Material {
            name: "custom",
            density: density * packing,
            surface_energy,
            flaw_size: grain_scale(props, formation),
            tensile_ratio: if cleaves {
                (0.04 + 0.96 * ductility).clamp(0.02, 1.0)
            } else {
                1.0
            },
            stiffness,
            // Creep begins around four tenths of the melting point in kelvin —
            // the homologous temperature rule — and nothing is left at it.
            thermal_onset: 0.4 * melting,
            thermal_gone: melting,
            // Taking a kilogram apart completely costs its cohesive energy.
            // An upper bound on the enthalpy of destruction rather than a
            // measurement of it: a real fire pyrolyses rather than atomises.
            destruction_enthalpy: props.cohesive_energy / props.unit_mass.max(1e-30),
            specific_heat,
            // No law. See the table above.
            resistivity: 1.0e6,
            combustible: false,
            ductility,
            metallic,
        }
    }

    /// Measure a material from several substances at once, by mass.
    ///
    /// A composite is not a substance and cannot be analysed as one. What it
    /// can be is the mass-weighted mixture of its parts, which is what a
    /// reinforced frame is and what any node holding more than one solid is.
    pub fn blend(parts: &[(Properties, Formation, f64)]) -> Material {
        let total: f64 = parts.iter().map(|(_, _, w)| w.max(0.0)).sum();
        if !(total > 0.0) {
            return Material::default();
        }
        let mut out: Option<Material> = None;
        for (props, formation, weight) in parts {
            let w = weight.max(0.0) / total;
            if w <= 0.0 {
                continue;
            }
            let m = Material::of(props, *formation);
            out = Some(match out {
                None => Material {
                    density: m.density * w,
                    surface_energy: m.surface_energy * w,
                    flaw_size: m.flaw_size * w,
                    tensile_ratio: m.tensile_ratio * w,
                    stiffness: m.stiffness * w,
                    thermal_onset: m.thermal_onset * w,
                    thermal_gone: m.thermal_gone * w,
                    destruction_enthalpy: m.destruction_enthalpy * w,
                    specific_heat: m.specific_heat * w,
                    resistivity: m.resistivity * w,
                    ductility: m.ductility * w,
                    metallic: m.metallic * w,
                    ..m
                },
                Some(a) => Material {
                    density: a.density + m.density * w,
                    surface_energy: a.surface_energy + m.surface_energy * w,
                    flaw_size: a.flaw_size + m.flaw_size * w,
                    tensile_ratio: a.tensile_ratio + m.tensile_ratio * w,
                    stiffness: a.stiffness + m.stiffness * w,
                    thermal_onset: a.thermal_onset + m.thermal_onset * w,
                    thermal_gone: a.thermal_gone + m.thermal_gone * w,
                    destruction_enthalpy: a.destruction_enthalpy + m.destruction_enthalpy * w,
                    specific_heat: a.specific_heat + m.specific_heat * w,
                    resistivity: a.resistivity + m.resistivity * w,
                    combustible: a.combustible || m.combustible,
                    ductility: a.ductility + m.ductility * w,
                    metallic: a.metallic + m.metallic * w,
                    ..a
                },
            });
        }
        out.unwrap_or_default()
    }

    /// What a node is made of, measured from its own matter.
    ///
    /// **This is the answer to "only a built thing has a surface".**
    /// `engine::World::surface_of` read `topology.material` or
    /// `morphology.material()` and returned `None` for everything else, so a
    /// rock, a boulder and a ball of wood could not collide — a *provenance
    /// test standing in for a state measurement*, which `PLAY.md` §2A measured
    /// and D13 retires.
    ///
    /// Only the **solid** pools count. A cup that is half water by mass
    /// presents ceramic, not an average of ceramic and water, and D17 made the
    /// phase a reading rather than a guess. Matter with no solid pool has no
    /// material and says so with `None`, which is D13's own line: a liquid's
    /// surface is its container's and a gas has none.
    pub fn measured(mixture: &Mixture, reg: &Registry, formation: Formation) -> Option<Material> {
        let mut parts: Vec<(Properties, Formation, f64)> = Vec::new();
        for pool in mixture.entries() {
            if pool.phase != Phase::Solid || !(pool.fraction > 0.0) {
                continue;
            }
            let Some(s) = reg.get(pool.substance) else { continue };
            parts.push((s.props, formation, pool.fraction));
        }
        if parts.is_empty() {
            return None;
        }
        Some(Material::blend(&parts))
    }

    /// Fraction of nominal strength remaining at a given temperature.
    pub fn strength_at(&self, temperature: f64) -> f64 {
        if temperature <= self.thermal_onset {
            1.0
        } else if temperature >= self.thermal_gone {
            0.0
        } else {
            1.0 - (temperature - self.thermal_onset) / (self.thermal_gone - self.thermal_onset)
        }
    }

    /// The stress this material fails at, Pa. **Griffith.**
    ///
    /// ```text
    ///     sigma = sqrt(2 E gamma / (pi a))
    /// ```
    ///
    /// `docs/PLAY.md` D14, and the problem it solves is not which formula to
    /// use. Theoretical strength is about `E/10` — how far atomic bonds stretch
    /// before letting go is a roughly fixed fraction of how stiff they are —
    /// and **real materials are one to three orders weaker**, because they fail
    /// at defects rather than everywhere at once. This engine models chemistry
    /// in depth and defects not at all.
    ///
    /// So strength is not a property of what a thing is made of. It is a
    /// property of **how well it was made**, chemistry cannot answer it and
    /// history can — which is why the one stored number here is a flaw scale
    /// and not a stress. See [`Material::flaw_size`].
    ///
    /// This replaced a `rupture` field holding eight tabulated stresses that
    /// nothing derived. Two of those eight were the *same substance* at
    /// different formation histories — masonry and bedrock are both silicate,
    /// green wood and dry timber are both cellulose — which a per-material
    /// stress could only ever record as four materials.
    ///
    /// # What it reproduces, measured against the table it replaces
    ///
    /// Derived without the retired values ever being shown to the derivation:
    ///
    /// ```text
    ///     material           a           derived    retired    ratio
    ///     dry timber       7.0e-4 m      2.87e7     7.0e7      0.41
    ///     green wood       2.0e-3        1.70e7     4.5e7      0.38
    ///     aragonite        3.0e-3        1.46e7     1.2e7      1.22
    ///     masonry          7.0e-2        2.88e6     2.0e6      1.44
    ///     ice              3.8e-2        3.85e6     1.7e6      2.27
    ///     reinforced       1.3e-3        1.89e7     1.8e8      0.11
    ///     steel            8.2e-3        2.87e6     4.0e8      0.0072
    ///     bedrock          2.1e-1        1.66e6     1.3e8      0.013
    /// ```
    ///
    /// Six of the eight are inside a factor of ten, and the **ordering of the
    /// deposited ones is exact**: dry timber above green wood above aragonite
    /// above masonry, which is the table's own order, from nothing but a ring
    /// width, a band and a course. Dry timber above green wood is the whole of
    /// "slow-grown timber is stronger", and it comes out of one number.
    ///
    /// # What it does not, and why — Griffith is the *brittle* law
    ///
    /// The two that miss by a hundred are exactly the two that **froze**, and
    /// the reason is the one D14 flags in a line: "a grain is also not the same
    /// as the worst flaw".
    ///
    /// - **Steel does not fail by Griffith at all.** It is ductile — the
    ///   derivation says so itself, at 0.67 — and a ductile metal yields by
    ///   moving dislocations, at a stress set by Hall and Petch on the grain
    ///   rather than by a crack running through it. A crack tip in steel blunts
    ///   instead of propagating, and Orowan's correction for the plastic work
    ///   at the tip is three to six orders of magnitude, which is the size of
    ///   the miss.
    /// - **Rock's flaw is intragranular.** Granite at 130 MPa implies a 34 µm
    ///   crack, against the 2 mm grains it has: the cracks are *inside* the
    ///   grains, not around them.
    ///
    /// Both want a second law rather than a correction factor, and neither is
    /// written here because inventing one to close a gap is what D14 rejects.
    /// Recorded in `docs/BACKLOG.md` with these numbers.
    pub fn strength(&self) -> f64 {
        let brittle = self.griffith_stress();
        let ductile = self.yield_stress();
        // A material is not one or the other: a reinforced frame is a sixth
        // iron by mass and five sixths silicate, and it fails partly each way.
        // The weight is `metallic` and **not** `ductility` — see the field's
        // own doc for what that distinction cost when it was missed.
        let d = self.metallic.clamp(0.0, 1.0);
        d * ductile + (1.0 - d) * brittle
    }

    /// Griffith's stress: what it takes to run a crack out of the worst flaw.
    ///
    /// The law for a **brittle** solid, which is one that has no way to blunt a
    /// crack tip. See [`Material::strength`] for how it combines with the
    /// ductile law.
    pub fn griffith_stress(&self) -> f64 {
        let a = self.flaw_size;
        if !(a > 0.0) || !(self.stiffness > 0.0) || !(self.surface_energy > 0.0) {
            return 0.0;
        }
        (2.0 * self.stiffness * self.surface_energy / (std::f64::consts::PI * a)).sqrt()
    }

    /// Yield stress: what it takes to make dislocations move.
    ///
    /// **Griffith does not describe a metal**, and the measurement says so
    /// loudly: run on steel's grain it gives 2.9 MPa against the 400 MPa the
    /// retired table held. A ductile metal does not fail by running a crack
    /// through a grain — the tip blunts, the material flows, and what sets the
    /// stress is how hard it is to move a dislocation.
    ///
    /// Two terms, both derived:
    ///
    /// ```text
    ///     sigma_y = G/100 + k / sqrt(d),    k = sqrt(pi G b sigma_0 / 2)
    /// ```
    ///
    /// The first is the **lattice resistance**. A perfect crystal shears at
    /// about `G/6`; a real one has dislocations already in it and resists at a
    /// hundredth of that, which is the standard Peierls estimate for a metal
    /// and is a universal fraction in the same family as Turnbull's 0.45
    /// rather than a per-material column. Steel comes out at 288 MPa from
    /// nothing but its own stiffness.
    ///
    /// The second is **Hall and Petch**: dislocations pile up at grain
    /// boundaries, so a fine-grained metal is stronger, and `k` follows from
    /// the pile-up model with the Burgers vector taken as the atomic spacing.
    /// It is the smaller term for steel — 0.6 MPa against 288 — which is
    /// itself the right answer: bcc iron's strength is dominated by its lattice
    /// resistance and only weakly by its grain size, which is why a cast bar
    /// and a rolled one differ by tens of per cent rather than by orders.
    ///
    /// So a ductile material keeps a *little* of "how well it was made" and a
    /// brittle one is almost all of it, which is the real asymmetry.
    pub fn yield_stress(&self) -> f64 {
        if !(self.stiffness > 0.0) {
            return 0.0;
        }
        // Shear modulus from Young's, at a Poisson's ratio of 0.3 — the value
        // for nearly every dense solid, and the reason nobody carries one.
        let shear = self.stiffness / 2.6;
        let lattice = shear / 100.0;
        let burgers = self.atomic_spacing();
        let d = self.flaw_size;
        if !(d > 0.0) || !(burgers > 0.0) {
            return lattice;
        }
        let k = (std::f64::consts::PI * shear * burgers * lattice / 2.0).sqrt();
        lattice + k / d.sqrt()
    }

    /// The spacing between atoms, from the surface energy and the stiffness.
    ///
    /// Not stored: making a surface costs about a quarter of an atom's cohesive
    /// energy over the area one atom presents, and stiffness is three times the
    /// energy density, so the two between them give the volume per atom and
    /// therefore its cube root. `12 gamma / E` to within the factors the two
    /// derivations share.
    pub fn atomic_spacing(&self) -> f64 {
        if !(self.stiffness > 0.0) {
            return 0.0;
        }
        12.0 * self.surface_energy / self.stiffness
    }

    /// Strength of a piece this big.
    ///
    /// A flaw cannot be larger than the thing that holds it, so a fibre is
    /// genuinely stronger than a bar of the same material — the size effect
    /// that a tabulated stress cannot express and that every member already
    /// carries a radius for.
    pub fn strength_at_size(&self, size: f64) -> f64 {
        let a = self.flaw_size.min(size.max(0.0));
        Material { flaw_size: a, ..*self }.strength()
    }
    pub fn thermal_limits(&self) -> (f64, f64) {
        (self.thermal_onset, self.thermal_gone)
    }
}


// ---------------------------------------------------------------------------
// how a thing was made
// ---------------------------------------------------------------------------

/// The conditions a piece of material formed under.
///
/// `docs/PLAY.md` D14's key move: **strength is not a property of what a thing
/// is made of, it is a property of how well it was made**, and this engine
/// records how things were made. Chemistry gives the cohesive energy and the
/// stiffness; this gives the flaw scale, and Griffith puts the two together.
///
/// # Two histories, two microstructures
///
/// There are exactly two ways a solid comes to exist here, and they leave
/// different flaws behind:
///
/// - **It froze.** Nuclei appeared in a melt and grew until they met, and the
///   flaw scale is the **grain**. What sets it is the competition between how
///   fast nuclei appear and how fast they grow, which [`grain_scale`] solves.
/// - **It was laid down.** A tree added a ring, a mason laid a course, a coral
///   deposited a layer. Nothing ever melted, there is no grain, and the flaw
///   scale is the **increment** — the thickness the process puts down at a
///   time. The generator knows it, because the generator is what laid it.
///
/// Neither carries a per-material constant. The same `Formation` applied to two
/// substances gives two different answers because the substances differ, which
/// is the difference between a rule and a table and the whole of why D14
/// rejected calibrating a flaw size per material against the `rupture` values
/// it replaces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Formation {
    /// Froze from its own melt, with heat leaving at these rates.
    ///
    /// **The undercooling is not an input.** It is solved for: a melt cools
    /// past its freezing point, nuclei accumulate, and the transformation runs
    /// away when the latent heat the growing grains release overtakes the heat
    /// being extracted. How far below the freezing point that happens is a
    /// consequence, and stating it would have been a per-material number by
    /// another name — which is what the first attempt at this did, and it moved
    /// strength by `10^3` for a 33% change in one input.
    Frozen {
        /// How fast the melt was cooling through its freezing point, K/s.
        cooling_rate: f64,
        /// How fast the solid front could advance, m/s. Limited by how fast
        /// latent heat leaves through the surface, not by interface kinetics:
        /// the collision-limited speed for iron is 470 m/s and real fronts are
        /// millimetres a second.
        front_speed: f64,
    },
    /// Laid down without ever having been molten — grown, precipitated, built.
    Deposited {
        /// The temperature it was laid down at, K.
        temperature: f64,
        /// How much of the volume the process actually filled, 0 to 1.
        ///
        /// **A deposited solid is usually not solid.** Wood is 600 kg/m^3 and
        /// the cellulose it is made of derives at 1612, so a tree fills 37% of
        /// the space it occupies and the rest is lumen and vessel. Masonry,
        /// concrete and a coral skeleton are the same story at different
        /// numbers, and a melt that froze is the one case where it is one.
        ///
        /// This is a property of the *object*, not of the substance, which is
        /// why it sits in the formation rather than in the chemistry — the same
        /// place and for the same reason as the flaw scale. For a node that is
        /// measured rather than stated it is derived, from the node's own bulk
        /// density against its substance's: see [`Formation::of_matter`].
        ///
        /// Gibson and Ashby: density goes as the packing and stiffness as its
        /// square, because a cellular solid bends its cell walls where a dense
        /// one stretches its bonds.
        packing: f64,
        /// The thickness the process lays down at a time, metres. **This is
        /// the flaw scale directly**: a growth ring, a course of blocks, a
        /// deposited layer. A tree adding a millimetre of radius a year has
        /// millimetre flaws, and Griffith on a millimetre gives green wood
        /// 2.4e7 Pa against the 4.5e7 the retired table held.
        increment: f64,
    },
}

impl Formation {
    /// Something that froze from its own melt.
    pub fn cooled(cooling_rate: f64, front_speed: f64) -> Formation {
        Formation::Frozen {
            cooling_rate: cooling_rate.max(0.0),
            front_speed: front_speed.max(0.0),
        }
    }

    /// Something laid down at a stated temperature, in increments this thick,
    /// filling this much of the space it occupies.
    pub fn laid_down(temperature: f64, increment: f64, packing: f64) -> Formation {
        Formation::Deposited {
            temperature,
            increment: increment.max(0.0),
            packing: packing.clamp(1e-6, 1.0),
        }
    }

    /// How much of the volume is actually solid. One for anything that froze.
    pub fn packing(self) -> f64 {
        match self {
            Formation::Frozen { .. } => 1.0,
            Formation::Deposited { packing, .. } => packing.clamp(1e-6, 1.0),
        }
    }

    /// Formation conditions for matter nobody made, **read off its own state**.
    ///
    /// `docs/PLAY.md` D14: "For matter that was sampled rather than made — a
    /// rock that simply exists in a scenario — the sampler draws formation
    /// conditions from the equilibrium the node is in, on the same
    /// maximum-entropy grounds it draws everything else, so nothing is authored
    /// and nothing is looked up."
    ///
    /// The equilibrium a node is in is its own radiative balance, and that
    /// supplies **both** rates a freezing needs, with nothing stated:
    ///
    /// ```text
    ///     P  = sigma A T^4                 what the node loses, Stefan-Boltzmann
    ///     dT/dt = P / (m c_p)              how fast it cools
    ///     u  = P / (A dH_v)                how fast the front can advance,
    ///                                      because latent heat leaves the same way
    /// ```
    ///
    /// Checked on a cubic metre of steel at 1500 K: it radiates 1.72 MW, which
    /// is 0.49 K/s and a front of 0.41 mm/s. Both are the right order for a
    /// casting, and neither was chosen.
    ///
    /// A node cold enough that its melting point is above its temperature never
    /// froze, so this reports a deposition instead, with the node's own size
    /// over its dynamical time as the increment — the one length and one rate a
    /// node has without being told.
    pub fn of_matter(m: &crate::state::Matter, props: &Properties) -> Formation {
        let area = 4.0 * std::f64::consts::PI * m.radius * m.radius;
        let power = crate::units::SIGMA_SB * area * m.temperature.powi(4);
        let atoms = props.atoms_per_unit.max(1) as f64;
        // Richard's rule again: melting costs about R per mole of atoms.
        let heat_of_fusion = K_B * N_AVOGADRO * props.melting_point * atoms
            / (props.unit_mass * N_AVOGADRO).max(1e-30);
        let specific_heat = 3.0 * K_B * atoms / props.unit_mass.max(1e-30);
        let heat_of_fusion_volumetric = heat_of_fusion * props.density.max(1e-6);
        let hot_enough = m.temperature >= props.melting_point;
        if hot_enough && power > 0.0 && m.mass > 0.0 && area > 0.0 {
            Formation::Frozen {
                cooling_rate: power / (m.mass * specific_heat),
                front_speed: power / (area * heat_of_fusion_volumetric.max(1e-30)),
            }
        } else {
            let t = m.dynamical_time();
            let increment = if t > 0.0 && t.is_finite() { m.radius / t } else { 0.0 };
            // How much of the node is actually solid, measured: its own bulk
            // density against the density the substance would have if it were
            // packed solid. A hollow box is mostly air and says so.
            let packing = (m.density() / props.density.max(1e-30)).clamp(1e-6, 1.0);
            Formation::Deposited { temperature: m.temperature, increment, packing }
        }
    }
}

/// The flaw scale of a piece of material, from the conditions it formed under.
///
/// `docs/PLAY.md` D14 rejects calibrating a flaw size per material against the
/// `rupture` values it replaces, and the reason is worth keeping: it would have
/// swapped a table of eight stresses for a table of eight lengths, and a number
/// reverse-engineered from the table it replaces is a frozen outcome that
/// nothing derived. **A shortcut may hold a rule rather than a constant**, and
/// this is the rule.
///
/// # A deposited solid: the increment
///
/// Nothing melted, so there is no grain. The flaw is the layer — the thickness
/// the process puts down at a time — and the generator knows it because the
/// generator is what laid it. That is the whole of the deposited branch, and it
/// is short because the information is already in the right place.
///
/// # A frozen solid: the grain, solved rather than stated
///
/// Classical nucleation theory, marched rather than thresholded. A melt cools
/// past its freezing point at `dT/dt`; at undercooling `x` the barrier to
/// forming a nucleus is
///
/// ```text
///     gamma_sl = 0.45 dH_f / (V_atom^(2/3) N_A)        Turnbull
///     dG*(x)   = 16 pi gamma_sl^3 T_m^2 / (3 dH_v^2 x^2)
///     I(x)     = (n_v k_B T / h) exp(-dG* / k_B T)
/// ```
///
/// so nuclei accumulate as `N(x) = integral I dx / (dT/dt)` and each one that
/// appeared at `x'` has grown to `u (x - x') / (dT/dt)`. The melt keeps
/// undercooling until **the latent heat the growing grains release overtakes
/// the heat being extracted** — recalescence — at which point the melt reheats
/// and nucleation stops. The grains present then are the grains it has:
///
/// ```text
///     dH_v * 4 pi N L^2 u  >=  rho c_p dT/dt        the melt reheats
///     d = N^(-1/3)                                  the grains that got in
/// ```
///
/// Every term comes from something already derived. `dH_f` is Richard's rule,
/// an entropy of fusion of about `R` per mole of atoms, which `chem::analyse`
/// already leans on for its melting points. `gamma_sl` is Turnbull's relation,
/// which reproduces iron's 0.20 J/m^2 without being shown it. The prefactor is
/// the classical one, atoms per volume times the attempt frequency `k_B T / h`.
///
/// **Why the first two attempts at this are worth recording.** Taking the
/// undercooling as an input made strength exponential in it — a 33% change in
/// one number moved steel by `10^3` — so the table of stresses was becoming a
/// table of temperatures. Solving it against a *first-nucleus* threshold
/// instead gave bedrock a 7.8 km grain, because the first nucleus appears at a
/// vanishing undercooling when cooling is slow. The integral is what fixes
/// both: nucleation is a runaway, and what matters is when it runs away.
///
/// **What this deliberately does not model**, because the engine has no
/// representation for it: work hardening, inclusions, and manufacturing
/// defects. Dislocations both weaken a material and strengthen it, and only the
/// first is captured here. A grain is also not the same as the worst flaw, so
/// this is a *scale* rather than a precise length.
pub fn grain_scale(p: &Properties, f: Formation) -> f64 {
    let atoms = p.atoms_per_unit.max(1) as f64;
    if !(p.density > 0.0) || !(p.unit_mass > 0.0) || !(p.melting_point > 0.0) {
        return 0.0;
    }
    let volume_atom = p.unit_mass / (p.density * atoms);
    let spacing = volume_atom.cbrt();

    let (cooling_rate, front_speed) = match f {
        // A layer is its own flaw. Nothing to solve.
        Formation::Deposited { increment, .. } => return increment.max(spacing),
        Formation::Frozen { cooling_rate, front_speed } => (cooling_rate, front_speed),
    };
    if !(cooling_rate > 0.0) || !(front_speed > 0.0) {
        return 0.0;
    }

    let melting = p.melting_point;
    let n_v = 1.0 / volume_atom;
    // Richard's rule: melting costs about R per mole of atoms of entropy.
    let heat_of_fusion_molar = K_B * N_AVOGADRO * melting;
    let heat_of_fusion_volumetric = heat_of_fusion_molar / (volume_atom * N_AVOGADRO);
    // Turnbull: the solid-liquid interface energy from the heat of fusion.
    let gamma_sl = 0.45 * heat_of_fusion_molar / (volume_atom.powf(2.0 / 3.0) * N_AVOGADRO);
    let specific_heat = 3.0 * K_B * atoms / p.unit_mass;
    // What the melt is giving up per second per cubic metre, which is what the
    // growing grains have to overtake.
    let extraction = p.density * specific_heat * cooling_rate;

    let rate_at = |x: f64| -> f64 {
        if !(x > 0.0) {
            return 0.0;
        }
        let barrier = 16.0 * std::f64::consts::PI * gamma_sl.powi(3) * melting * melting
            / (3.0 * heat_of_fusion_volumetric * heat_of_fusion_volumetric * x * x);
        let t = (melting - x).max(1.0);
        let r = n_v * K_B * t / H_PLANCK * (-barrier / (K_B * t)).exp();
        if r.is_finite() { r } else { 0.0 }
    };

    // March the undercooling. A thousand steps over nine tenths of the melting
    // point is finer than the answer deserves — the nucleation rate climbs by
    // tens of orders over the interesting range, so the crossing is sharp and
    // its position is insensitive to the step.
    const STEPS: usize = 4000;
    let dx = 0.9 * melting / STEPS as f64;
    // A nucleus born at undercooling `x'` has grown to `L(x) - L(x')` by the
    // time the melt reaches `x`, so the total grain surface — which is what
    // sets the rate latent heat comes out — is
    //
    //     sum_i N_i (L - L_i)^2 = L^2 A0 - 2 L A1 + A2
    //
    // and the three moments are all the state the march needs. Growing every
    // nucleus from `L = 0` instead, which is the obvious way to write this,
    // credits grains with the growth of ones that did not exist yet: it fired
    // recalescence at a vanishing undercooling and gave bedrock a 4.4 m grain.
    let mut a0 = 0.0f64;
    let mut a1 = 0.0f64;
    let mut a2 = 0.0f64;
    let mut length = 0.0f64;
    let mut found = 0.0f64;
    for i in 1..=STEPS {
        let x = i as f64 * dx;
        length += front_speed * dx / cooling_rate;
        let born = rate_at(x) * dx / cooling_rate;
        if born > 0.0 {
            a0 += born;
            a1 += born * length;
            a2 += born * length * length;
        }
        if a0 <= 0.0 {
            continue;
        }
        let area = (length * length * a0 - 2.0 * length * a1 + a2).max(0.0);
        let release =
            heat_of_fusion_volumetric * 4.0 * std::f64::consts::PI * area * front_speed;
        if release >= extraction {
            found = a0;
            break;
        }
    }
    if !(found > 0.0) {
        // It never ran away: the melt passed its whole freezing range without
        // the grains catching up with the heat leaving. That is a glass, and a
        // glass has no grains — the caller gets zero, which
        // `Material::strength` reads as "no answer" rather than as "infinitely
        // strong".
        return 0.0;
    }
    let d = found.powf(-1.0 / 3.0);
    if !d.is_finite() {
        return 0.0;
    }
    // A grain cannot be smaller than the atoms in it. Nothing caps it from
    // above here: a flaw bigger than the object is what
    // `Material::strength_at_size` is for, and it needs the object's size,
    // which a material does not have.
    d.max(spacing)
}

/// The substances the named presets stand for.
///
/// **Arrangements, not materials.** Every number this module produces comes out
/// of `chem::analyse`, which takes atoms and bonds and nothing else, so what is
/// written here is which atoms and which bonds — the level `CLAUDE.md` calls
/// axiom-side. "Steel is iron" is a statement about a *preset*, in the same
/// family as a scenario saying a planet is mostly silicate; it is not a table
/// of what iron is like.
///
/// Analysed once each, on first use.
pub mod substances {
    use super::*;

    fn once(
        slot: &'static std::sync::OnceLock<Properties>,
        build: impl Fn() -> Arrangement,
    ) -> Properties {
        *slot.get_or_init(|| {
            let arr = build();
            analyse(&arr).unwrap_or_else(|e| panic!("preset substance does not analyse: {e:?}"))
        })
    }

    /// Iron. Four atoms in a cell, each bonded to the other three.
    ///
    /// A metal is the one thing this chemistry cannot represent properly: iron
    /// really has eight nearest neighbours and a valence of three, so the
    /// arrangement takes the largest coordination the valence model allows.
    /// The cohesive energy comes out about 2.3x low as a result — 1.84 eV per
    /// atom against 4.28 — and everything derived from it is low with it. Said
    /// here rather than corrected with a factor, because a factor applied to
    /// metals and nothing else is a species table with one entry.
    pub fn iron() -> Properties {
        static P: std::sync::OnceLock<Properties> = std::sync::OnceLock::new();
        once(&P, iron_arrangement)
    }

    /// The atoms and bonds `iron` analyses, for a caller that needs to intern it
    /// in a world's own `Registry`.
    pub fn iron_arrangement() -> Arrangement {
        {
            Arrangement::crystal(
                vec![Element(26); 4],
                vec![
                    Bond::new(0, 1, Order::Single),
                    Bond::new(0, 2, Order::Single),
                    Bond::new(0, 3, Order::Single),
                    Bond::new(1, 2, Order::Single),
                    Bond::new(1, 3, Order::Single),
                    Bond::new(2, 3, Order::Single),
                ],
                // Four atoms to the cell at 7,870 kg/m^3.
                Lattice::Cubic { a: 3.612e-10 },
            )
        }
    }

    /// Quartz: one silicon tetrahedrally bonded to two bridging oxygens.
    pub fn silica() -> Properties {
        static P: std::sync::OnceLock<Properties> = std::sync::OnceLock::new();
        once(&P, silica_arrangement)
    }

    /// The atoms and bonds `silica` analyses, for a caller that needs to intern it
    /// in a world's own `Registry`.
    pub fn silica_arrangement() -> Arrangement {
        {
            Arrangement::crystal(
                vec![Element(14), Element(8), Element(8)],
                vec![
                    Bond::new(0, 1, Order::Single),
                    Bond::new(0, 1, Order::Single),
                    Bond::new(0, 2, Order::Single),
                    Bond::new(0, 2, Order::Single),
                ],
                // One SiO2 to the cell at 2,650 kg/m^3.
                Lattice::Cubic { a: 3.354e-10 },
            )
        }
    }

    /// Calcium carbonate: a carbonate ion with a calcium holding it.
    pub fn calcium_carbonate() -> Properties {
        static P: std::sync::OnceLock<Properties> = std::sync::OnceLock::new();
        once(&P, calcium_carbonate_arrangement)
    }

    /// The atoms and bonds `calcium_carbonate` analyses, for a caller that needs to intern it
    /// in a world's own `Registry`.
    pub fn calcium_carbonate_arrangement() -> Arrangement {
        {
            Arrangement::crystal(
                vec![Element(20), Element(6), Element(8), Element(8), Element(8)],
                vec![
                    Bond::new(1, 2, Order::Double),
                    Bond::new(1, 3, Order::Single),
                    Bond::new(1, 4, Order::Single),
                    Bond::new(0, 3, Order::Ionic),
                    Bond::new(0, 4, Order::Ionic),
                ],
                // One CaCO3 to the cell at 2,930 kg/m^3.
                Lattice::Cubic { a: 3.845e-10 },
            )
        }
    }

    /// Water, which is ice wherever it is cold enough.
    ///
    /// **A molecule and not a crystal**, deliberately, and the reason is the
    /// one `chem::registry::Phase` already gives: ice and water are the same
    /// substance, and a registry that stored "solid" against water would store
    /// water twice. What melts in ice is the hydrogen bonding between whole
    /// molecules, not the O-H bonds inside them, and `analyse`'s molecular
    /// branch is the one that prices that — through Trouton's rule and the
    /// hydrogen-bond count. Analysed as a lattice it comes out melting at
    /// 893 K, because the lattice branch is pricing the covalent bonds.
    pub fn water() -> Properties {
        static P: std::sync::OnceLock<Properties> = std::sync::OnceLock::new();
        once(&P, water_arrangement)
    }

    /// The atoms and bonds `water` analyses, for a caller that needs to intern it
    /// in a world's own `Registry`.
    pub fn water_arrangement() -> Arrangement {
        {
            Arrangement::molecule(
                vec![Element(8), Element(1), Element(1)],
                vec![Bond::new(0, 1, Order::Single), Bond::new(0, 2, Order::Single)],
            )
        }
    }

    /// One glucose ring of a cellulose chain, which is what wood is.
    ///
    /// A polymer is not a crystal and not a small molecule, and this takes the
    /// repeating unit — the thing the chain is made of — as the formula unit.
    /// That is the same choice `Lattice::Cubic` makes for a salt.
    pub fn cellulose() -> Properties {
        static P: std::sync::OnceLock<Properties> = std::sync::OnceLock::new();
        once(&P, cellulose_arrangement)
    }

    /// The atoms and bonds `cellulose` analyses, for a caller that needs to intern it
    /// in a world's own `Registry`.
    pub fn cellulose_arrangement() -> Arrangement {
        {
            // C6H10O5: six carbons in a ring closed through an oxygen, with
            // four hydroxyls and the rest hydrogens.
            let mut atoms = vec![Element(6); 6];
            atoms.extend(std::iter::repeat_n(Element(8), 5));
            atoms.extend(std::iter::repeat_n(Element(1), 10));
            let mut bonds = Vec::new();
            for i in 0..5 {
                bonds.push(Bond::new(i, i + 1, Order::Single));
            }
            // The ring closes through an oxygen rather than carbon to carbon,
            // which is what makes it a sugar.
            bonds.push(Bond::new(5, 6, Order::Single));
            bonds.push(Bond::new(6, 0, Order::Single));
            for i in 0..4 {
                bonds.push(Bond::new(i + 1, 7 + i, Order::Single));
                bonds.push(Bond::new(7 + i, 11 + i, Order::Single));
            }
            for i in 0..6 {
                bonds.push(Bond::new(i, 15 + i, Order::Single));
            }
            Arrangement::molecule(atoms, bonds)
        }
    }
}
