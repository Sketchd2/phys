p = r"E:\Projects\phys\src\liquid.rs"
s = open(p, newline="").read()
nl = "\r\n" if "\r\n" in s else "\n"
s = s.replace("\r\n", "\n")
marker = "/// The quintic switch: 1 below `on`, 0 above `cut`"
assert marker in s
add = '''/// Everything the site-site pair loop of a flexible liquid needs, for an engine
/// that does it elsewhere (the GPU): the molecules' sites relative to their
/// centres, each site's charge and Coulomb width, and the pairs of molecules
/// inside the cut-off with the centres' separation and the switch
/// (`induction::Link`: `weight` is the switch, `dweight` its slope).
pub struct PairProblem<'a> {
    pub table: &'a SiteSite,
    /// The site type of each of a molecule's sites.
    pub types: &'a [usize],
    /// Per molecule, each site's place relative to the molecule's centre.
    pub rel: &'a [Vec<Vec3>],
    pub charge: &'a [Vec<f64>],
    pub sigma: &'a [Vec<f64>],
    pub links: &'a [crate::induction::Link],
}

/// What the pair loop gives back.
pub struct PairResult {
    pub energy: f64,
    pub virial: f64,
    /// Per molecule, the force on each site.
    pub fsite: Vec<Vec<Vec3>>,
    /// Per molecule, the switch's force on its centre.
    pub fcentre: Vec<Vec3>,
    /// Per molecule, per site, the energy's derivative with respect to the
    /// site's charge and to its Coulomb width.
    pub gq: Vec<Vec<f64>>,
    pub gs: Vec<Vec<f64>>,
}

/// A way of doing the pair loop faster than `Flex::forces` does, to the same
/// answer to round-off (its sums are added in another, fixed, order).
pub trait PairEngine: Send + Sync {
    fn pairs(&self, problem: &PairProblem) -> Option<PairResult>;
}

static PAIR_ENGINE: std::sync::Mutex<Option<std::sync::Arc<dyn PairEngine>>> = std::sync::Mutex::new(None);

/// Install a pair-loop engine for this process (`None` removes it).
pub fn set_pair_engine(engine: Option<std::sync::Arc<dyn PairEngine>>) {
    *PAIR_ENGINE.lock().expect("the pair engine") = engine;
}

/// The installed pair-loop engine, if any.
pub fn pair_engine() -> Option<std::sync::Arc<dyn PairEngine>> {
    PAIR_ENGINE.lock().expect("the pair engine").clone()
}

'''
s = s.replace(marker, add + marker, 1)
open(p, "w", newline="").write(s.replace("\n", nl))
print("liquid.rs patched")
