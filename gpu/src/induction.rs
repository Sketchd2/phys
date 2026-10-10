//! The induced-dipole solve of `phys::induction::solve` on the GPU, in double
//! precision — `docs/PLAY.md` E8c, the flexible liquid's cost.
//!
//! The reference spends about two thirds of a liquid step in its solve: the
//! screened tensors (erf-heavy, independent), twelve conjugate-gradient
//! iterations of a matrix-vector product, and the forces. Here the tensors, the
//! whole CG loop and the forces stay on the card, and one submission does a
//! step's solve and brings the dipoles, forces and virial back. The buffers are
//! kept from one call to the next (a liquid asks for nearly the same sizes every
//! step), so a call uploads positions, links and the starting dipoles, and
//! nothing is allocated.
//!
//! **Deterministic.** Every kernel is a gather in a fixed order and every
//! reduction a fixed tree, so a run gives the same bits every time; it differs
//! from the reference's chunked sums by round-off, which the tests bound.
//! `gexp` is built in the shader because WGSL has no `f64` exp.
//!
//! **What it takes.** A cluster of identical molecules (the same number of
//! polarisable sites and charges each), which is every liquid of one substance;
//! anything else is declined and the reference does it.

use phys::induction::{Cluster, Induced, InductionEngine, Link};
use phys::math::Vec3;
use std::sync::Mutex;

/// Conjugate-gradient iterations one submission runs at most.
const BATCH: usize = 24;
/// The reference's limit.
const MAX_ITERATIONS: u32 = 200;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    nm: u32,
    np: u32,
    nch: u32,
    nl: u32,
    n: u32,
    nc: u32,
    npairs: u32,
    max_iterations: u32,
}

const NAMES: [&str; 11] = ["k_terms_dd", "k_terms_cd", "k_terms_intra", "k_e0", "k_matvec_mu", "k_matvec_p", "k_cg_init", "k_cg_b", "k_link_force", "k_gather_pol", "k_gather_charge"];

/// The buffers for a liquid of `nm` molecules and room for `cap` links.
struct Work {
    nm: usize,
    np: usize,
    nch: usize,
    cap: usize,
    params: wgpu::Buffer,
    bufs: Vec<wgpu::Buffer>,
    bind: wgpu::BindGroup,
    staging: wgpu::Buffer,
    offsets: Vec<usize>,
}

// Indices into `Work::bufs`, in the shader's binding order from 1.
const POL_POS: usize = 0;
const POL_SIG: usize = 1;
const INV_ALPHA: usize = 2;
const CH_POS: usize = 3;
const CH_Q: usize = 4;
const CH_SIG: usize = 5;
const LF: usize = 6;
const LIJ: usize = 7;
const INC: usize = 8;
const E0: usize = 15;
const MU: usize = 16;
const SC: usize = 21;
const FPOL: usize = 24;
const FCH: usize = 25;
const LFORCE: usize = 26;
const LVIR: usize = 27;

pub struct GpuInduction {
    device: wgpu::Device,
    queue: wgpu::Queue,
    layout: wgpu::BindGroupLayout,
    pipelines: Vec<wgpu::ComputePipeline>,
    pub name: String,
    work: Mutex<Option<Work>>,
    /// Seconds on the card (submit to readback), in all, and the number of calls, for timing.
    pub timing: Mutex<(f64, f64, usize)>,
    profiled: std::sync::atomic::AtomicBool,
}

impl GpuInduction {
    /// The first GPU wgpu finds that offers double precision in shaders, or why not.
    pub fn new() -> Result<GpuInduction, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions { power_preference: wgpu::PowerPreference::HighPerformance, force_fallback_adapter: false, compatible_surface: None, ..Default::default() })).map_err(|e| format!("no adapter: {e}"))?;
        if !adapter.features().contains(wgpu::Features::SHADER_F64) {
            return Err("the adapter has no double precision in shaders".into());
        }
        let info = adapter.get_info();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor { label: Some("phys-gpu-induction"), required_features: wgpu::Features::SHADER_F64, required_limits: adapter.limits(), ..Default::default() })).map_err(|e| format!("no device: {e}"))?;
        // One layout for every kernel: binding 0 a uniform, 1..=28 storage.
        let mut entries = vec![wgpu::BindGroupLayoutEntry { binding: 0, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None }, count: None }];
        for b in 1..=28u32 {
            entries.push(wgpu::BindGroupLayoutEntry { binding: b, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: false }, has_dynamic_offset: false, min_binding_size: None }, count: None });
        }
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor { label: Some("induction"), entries: &entries });
        let pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor { label: Some("induction"), bind_group_layouts: &[Some(&layout)], immediate_size: 0 });
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("induction"), source: wgpu::ShaderSource::Wgsl(include_str!("induction.wgsl").into()) });
        let pipelines = NAMES
            .iter()
            .map(|n| device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: Some(n), layout: Some(&pl), module: &module, entry_point: Some(n), compilation_options: Default::default(), cache: None }))
            .collect();
        Ok(GpuInduction { device, queue, layout, pipelines, name: format!("gpu induction (f64, {} via {:?})", info.name, info.backend), work: Mutex::new(None), timing: Mutex::new((0.0, 0.0, 0)), profiled: std::sync::atomic::AtomicBool::new(false) })
    }

    fn make_work(&self, nm: usize, np: usize, nch: usize, cap: usize) -> Work {
        let (n, nc) = (nm * np, nm * nch);
        let npairs = np * (np - 1) / 2;
        let dd_n = cap * np * np;
        let cd_n = cap * 2 * np * nch;
        let in_n = nm * npairs;
        // Eight-byte units in each buffer, in the shader's binding order from 1.
        let sizes: [usize; 28] = [
            3 * n,                      // pol_pos
            n,                          // pol_sig
            n,                          // inv_alpha
            3 * nc,                     // ch_pos
            nc,                         // ch_q
            nc,                         // ch_sig
            8 * cap,                    // lf
            cap,                        // lij (two u32 a link)
            nm / 2 + cap + 2,           // inc (u32: nm + 1 starts, then two items a link)
            4 * dd_n,                   // dd_t
            3 * dd_n,                   // dd_r
            2 * cd_n,                   // cd_h
            3 * cd_n,                   // cd_r
            4 * in_n,                   // in_t
            3 * in_n,                   // in_r
            3 * n,                      // e0
            3 * n,                      // mu
            3 * n,                      // rv
            3 * n,                      // pv
            3 * n,                      // apv
            3 * n,                      // zv
            16,                         // sc
            3 * dd_n,                   // dd_f
            3 * cd_n,                   // cd_f
            3 * n,                      // fpol
            3 * nc,                     // fch
            3 * cap,                    // lforce
            cap,                        // lvir
        ];
        let usage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST;
        let bufs: Vec<wgpu::Buffer> = sizes.iter().map(|&s| self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("induction"), size: (s.max(1) * 8) as u64, usage, mapped_at_creation: false })).collect();
        let params = self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("params"), size: std::mem::size_of::<Params>() as u64, usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let mut entries = vec![wgpu::BindGroupEntry { binding: 0, resource: params.as_entire_binding() }];
        for (k, b) in bufs.iter().enumerate() {
            entries.push(wgpu::BindGroupEntry { binding: 1 + k as u32, resource: b.as_entire_binding() });
        }
        let bind = self.device.create_bind_group(&wgpu::BindGroupDescriptor { label: Some("induction"), layout: &self.layout, entries: &entries });
        // Where each output sits in the one staging buffer (doubles): mu, e0, fpol, fch, lforce, lvir, sc.
        let out_sizes = [3 * n, 3 * n, 3 * n, 3 * nc, 3 * cap, cap, 16];
        let mut offsets = Vec::new();
        let mut total = 0usize;
        for s in out_sizes {
            offsets.push(total);
            total += s.max(1);
        }
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("staging"), size: (total * 8) as u64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        Work { nm, np, nch, cap, params, bufs, bind, staging, offsets }
    }

    fn run(&self, cl: &Cluster, links: &[Link], warm: Option<&[Vec<Vec3>]>, tolerance: f64) -> Option<Induced> {
        let t_all = std::time::Instant::now();
        let nm = cl.pol.len();
        if nm == 0 || links.is_empty() {
            return None;
        }
        let (np, nch) = (cl.pol[0].len(), cl.charges[0].len());
        if np == 0 || cl.pol.iter().any(|m| m.len() != np) || cl.charges.iter().any(|m| m.len() != nch) {
            return None;
        }
        let (n, nc, nl) = (nm * np, nm * nch, links.len());
        let npairs = np * (np - 1) / 2;
        let mut guard = self.work.lock().expect("the induction work");
        let rebuild = match guard.as_ref() {
            Some(w) => w.nm != nm || w.np != np || w.nch != nch || w.cap < nl,
            None => true,
        };
        if rebuild {
            *guard = Some(self.make_work(nm, np, nch, nl + nl / 4 + 64));
        }
        let w = guard.as_ref().expect("the work");
        let q = &self.queue;
        let put = |i: usize, bytes: &[u8]| q.write_buffer(&w.bufs[i], 0, bytes);
        let f3 = |it: &mut dyn Iterator<Item = Vec3>| -> Vec<f64> { it.flat_map(|p| [p.x, p.y, p.z]).collect() };
        put(POL_POS, bytemuck::cast_slice(&f3(&mut cl.pol.iter().flat_map(|m| m.iter().map(|p| p.pos)))));
        let sigma: Vec<f64> = cl.pol.iter().flat_map(|m| m.iter().map(|p| phys::induction::sigma_of(p.alpha))).collect();
        put(POL_SIG, bytemuck::cast_slice(&sigma));
        let inv_alpha: Vec<f64> = cl.pol.iter().flat_map(|m| m.iter().map(|p| 1.0 / p.alpha)).collect();
        put(INV_ALPHA, bytemuck::cast_slice(&inv_alpha));
        if nc > 0 {
            put(CH_POS, bytemuck::cast_slice(&f3(&mut cl.charges.iter().flat_map(|m| m.iter().map(|c| c.pos)))));
            let chq: Vec<f64> = cl.charges.iter().flat_map(|m| m.iter().map(|c| c.q)).collect();
            put(CH_Q, bytemuck::cast_slice(&chq));
            let chs: Vec<f64> = cl.charges.iter().flat_map(|m| m.iter().map(|c| c.sigma)).collect();
            put(CH_SIG, bytemuck::cast_slice(&chs));
        }
        let mut lf = Vec::with_capacity(nl * 8);
        let mut lij: Vec<u32> = Vec::with_capacity(nl * 2);
        for l in links {
            lf.extend_from_slice(&[l.shift.x, l.shift.y, l.shift.z, l.weight, l.dweight, l.d.x, l.d.y, l.d.z]);
            lij.extend_from_slice(&[l.i as u32, l.j as u32]);
        }
        put(LF, bytemuck::cast_slice(&lf));
        put(LIJ, bytemuck::cast_slice(&lij));
        // Each molecule's links: [nm + 1 starts][items k * 2 + side].
        let mut count = vec![0u32; nm + 1];
        for l in links {
            count[l.i + 1] += 1;
            count[l.j + 1] += 1;
        }
        for m in 0..nm {
            count[m + 1] += count[m];
        }
        let mut inc = vec![0u32; nm + 1 + 2 * nl];
        inc[..nm + 1].copy_from_slice(&count);
        let mut fill = count.clone();
        for (k, l) in links.iter().enumerate() {
            inc[nm + 1 + fill[l.i] as usize] = (2 * k) as u32;
            fill[l.i] += 1;
            inc[nm + 1 + fill[l.j] as usize] = (2 * k + 1) as u32;
            fill[l.j] += 1;
        }
        put(INC, bytemuck::cast_slice(&inc));
        let mu0: Vec<f64> = match warm {
            Some(wm) if wm.iter().map(|m| m.len()).sum::<usize>() == n => wm.iter().flatten().flat_map(|p| [p.x, p.y, p.z]).collect(),
            _ => vec![0.0; 3 * n],
        };
        put(MU, bytemuck::cast_slice(&mu0));
        let mut sc0 = vec![0.0f64; 16];
        sc0[0] = tolerance;
        put(SC, bytemuck::cast_slice(&sc0));
        let params = Params { nm: nm as u32, np: np as u32, nch: nch as u32, nl: nl as u32, n: n as u32, nc: nc as u32, npairs: npairs as u32, max_iterations: MAX_ITERATIONS };
        q.write_buffer(&w.params, 0, bytemuck::bytes_of(&params));

        let (dd_n, cd_n, in_n) = (nl * np * np, nl * 2 * np * nch, nm * npairs);
        let groups = |count: usize, size: usize| count.div_ceil(size) as u32;
        // `PHYS_GPU_PROFILE=1`: once, run each kernel on its own and say how long the card took
        // (each includes a submission's fixed cost, so read them against one another).
        if std::env::var("PHYS_GPU_PROFILE").is_ok() && !self.profiled.swap(true, std::sync::atomic::Ordering::SeqCst) {
            let stages: Vec<(usize, u32, &str)> = vec![
                (0, groups(dd_n, 64), "terms dd"),
                (1, groups(cd_n, 64), "terms cd"),
                (2, groups(in_n, 64), "terms intra"),
                (3, n as u32, "e0"),
                (4, n as u32, "matvec"),
                (6, 1, "cg init"),
                (5, n as u32, "matvec p"),
                (7, 1, "cg b"),
                (8, groups(nl, 64), "link force"),
                (9, n as u32, "gather pol"),
                (10, nc as u32, "gather charge"),
            ];
            for rep in 0..2 {
                let mut line = String::new();
                for (which, g, label) in &stages {
                    let mut enc = self.device.create_command_encoder(&Default::default());
                    {
                        let mut pass = enc.begin_compute_pass(&Default::default());
                        pass.set_bind_group(0, &w.bind, &[]);
                        pass.set_pipeline(&self.pipelines[*which]);
                        pass.dispatch_workgroups(*g, 1, 1);
                    }
                    let t = std::time::Instant::now();
                    self.queue.submit(Some(enc.finish()));
                    self.device.poll(wgpu::PollType::wait_indefinitely()).expect("the GPU finished");
                    line += &format!("{label} {:.2} ms; ", t.elapsed().as_secs_f64() * 1000.0);
                }
                if rep == 1 {
                    eprintln!("gpu profile ({n} sites, {nl} links, {dd_n} dd and {cd_n} cd entries): {line}");
                }
            }
            // The profile ran kernels out of order; start the dipoles and flags afresh.
            put(MU, bytemuck::cast_slice(&mu0));
            put(SC, bytemuck::cast_slice(&sc0));
        }

        let mut first = true;
        let mut t_gpu = 0.0;
        loop {
            let mut enc = self.device.create_command_encoder(&Default::default());
            {
                let mut pass = enc.begin_compute_pass(&Default::default());
                pass.set_bind_group(0, &w.bind, &[]);
                let mut go = |pass: &mut wgpu::ComputePass, which: usize, g: u32| {
                    pass.set_pipeline(&self.pipelines[which]);
                    pass.dispatch_workgroups(g, 1, 1);
                };
                if first {
                    go(&mut pass, 0, groups(dd_n, 64));
                    go(&mut pass, 1, groups(cd_n, 64));
                    go(&mut pass, 2, groups(in_n, 64));
                    go(&mut pass, 3, n as u32);
                    go(&mut pass, 4, n as u32);
                    go(&mut pass, 6, 1);
                    first = false;
                }
                for _ in 0..BATCH {
                    go(&mut pass, 5, n as u32);
                    go(&mut pass, 7, 1);
                }
                go(&mut pass, 8, groups(nl, 64));
                go(&mut pass, 9, n as u32);
                if nc > 0 {
                    go(&mut pass, 10, nc as u32);
                }
            }
            let out: [(usize, usize); 7] = [(MU, 3 * n), (E0, 3 * n), (FPOL, 3 * n), (FCH, 3 * nc), (LFORCE, 3 * nl), (LVIR, nl), (SC, 16)];
            for (k, (buf, len)) in out.iter().enumerate() {
                if *len > 0 {
                    enc.copy_buffer_to_buffer(&w.bufs[*buf], 0, &w.staging, (w.offsets[k] * 8) as u64, (*len * 8) as u64);
                }
            }
            let t = std::time::Instant::now();
            self.queue.submit(Some(enc.finish()));
            let slice = w.staging.slice(..);
            slice.map_async(wgpu::MapMode::Read, |r| r.expect("a result mapped"));
            self.device.poll(wgpu::PollType::wait_indefinitely()).expect("the GPU finished");
            t_gpu += t.elapsed().as_secs_f64();
            let data: Vec<f64> = bytemuck::cast_slice(&slice.get_mapped_range().expect("the results range")).to_vec();
            w.staging.unmap();
            let sc = &data[w.offsets[6]..w.offsets[6] + 16];
            let iterations = sc[4] as usize;
            if sc[3] == 0.0 {
                continue;
            }
            let vec3s = |off: usize, count: usize| -> Vec<Vec3> { (0..count).map(|i| Vec3 { x: data[off + 3 * i], y: data[off + 3 * i + 1], z: data[off + 3 * i + 2] }).collect() };
            let mu = vec3s(w.offsets[0], n);
            let e0 = vec3s(w.offsets[1], n);
            let energy = -0.5 * mu.iter().zip(&e0).map(|(m, e)| m.dot(*e)).sum::<f64>();
            let fp = vec3s(w.offsets[2], n);
            let fc = vec3s(w.offsets[3], nc);
            let link_force = vec3s(w.offsets[4], nl);
            let virial: f64 = data[w.offsets[5]..w.offsets[5] + nl].iter().sum();
            let split = |v: &[Vec3], per: usize| -> Vec<Vec<Vec3>> { (0..nm).map(|i| v[i * per..(i + 1) * per].to_vec()).collect() };
            let mut tm = self.timing.lock().expect("the timing");
            tm.0 += t_gpu;
            tm.1 += t_all.elapsed().as_secs_f64();
            tm.2 += 1;
            return Some(Induced { energy, dipoles: split(&mu, np), force_pol: split(&fp, np), force_charge: split(&fc, nch), link_force, virial, iterations });
        }
    }
}

impl InductionEngine for GpuInduction {
    fn solve(&self, cl: &Cluster, links: &[Link], warm: Option<&[Vec<Vec3>]>, tolerance: f64) -> Option<Induced> {
        self.run(cl, links, warm, tolerance)
    }
}

/// Install the GPU's induction solve for this process, returning its name; or
/// why not, in which case the reference stays.
pub fn install_induction() -> Result<String, String> {
    let gpu = GpuInduction::new()?;
    let name = gpu.name.clone();
    phys::induction::set_induction_engine(Some(std::sync::Arc::new(gpu)));
    Ok(name)
}
