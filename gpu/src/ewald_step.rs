//! A whole Ewald force evaluation's long-range work on the card, in double precision —
//! `docs/PLAY.md` E8e, the owner's standing instruction (2026-10-18) to port what is faster.
//!
//! `GpuEwald` (`ewald.rs`) put the reciprocal sums on the card and gained 15 %: a
//! sub-millisecond kernel every few milliseconds, with the CPU doing the real-space work in
//! between, never raises the card's clock (540 of 2145 MHz was measured), so each submission cost
//! five times its benchmark (`docs/SOLUTIONS.md` 9.2). This keeps the **whole** of it resident:
//! the short-range tensors (the smeared tensor less the one at the Ewald width), the molecule's own
//! charge-dipole terms taken off, the permanent charges' reciprocal field, the conjugate-gradient
//! loop with the reciprocal dipole response added to every product, the forces, and the
//! reciprocal forces and arm sums. One submission does the solve (as `induction.rs` does for the
//! truncated law), a second the reciprocal forces once the dipoles have converged, and a few
//! small arrays come back.
//!
//! The shader is `ewald_step.wgsl`, generated from `ewald_step.wgsl.tpl`, itself the validated
//! `induction.wgsl.tpl` with the Ewald changes and the reciprocal kernels added; the flexible
//! liquid's engine is untouched.
//!
//! It answers to the reference's tolerance, not bit for bit (its sums are added in another,
//! fixed, order), and gives the same bits twice. A cluster of identical molecules, as there.

use phys::ewald::{set_step_engine, Ewald, EwaldStep, EwaldStepEngine, Finished, C};
use phys::induction::{Cluster, Induced, Link};
use phys::math::Vec3;
use std::sync::{Arc, Mutex};

const PI: f64 = std::f64::consts::PI;
/// Conjugate-gradient iterations one submission runs at most.
const BATCH: usize = 24;
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

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct RParams {
    nc: u32,
    nk: u32,
    n0: u32,
    n1: u32,
    n2: u32,
    pad: [u32; 3],
}

const NAMES: [&str; 23] = [
    "k_terms_dd", "k_terms_cd", "k_terms_intra", "k_e0", "k_matvec_mu", "k_matvec_p", "k_cg_init", "k_cg_b", "k_link_force", "k_gather_pol", "k_gather_charge", "k_terms_icd", "r_sq", "r_dsum_mu", "r_dsum_p", "r_stot", "r_field", "r_resp_mu", "r_resp_p", "r_gch", "r_gdpos",
    "r_arm_c", "r_arm_d",
];
const K_TERMS_DD: usize = 0;
const K_TERMS_CD: usize = 1;
const K_TERMS_INTRA: usize = 2;
const K_E0: usize = 3;
const K_MATVEC_MU: usize = 4;
const K_MATVEC_P: usize = 5;
const K_CG_INIT: usize = 6;
const K_CG_B: usize = 7;
const K_LINK_FORCE: usize = 8;
const K_GATHER_POL: usize = 9;
const K_GATHER_CHARGE: usize = 10;
const K_TERMS_ICD: usize = 11;
const R_SQ: usize = 12;
const R_DSUM_MU: usize = 13;
const R_DSUM_P: usize = 14;
const R_STOT: usize = 15;
const R_FIELD: usize = 16;
const R_RESP_MU: usize = 17;
const R_RESP_P: usize = 18;
const R_GCH: usize = 19;
const R_GDPOS: usize = 20;
const R_ARM_C: usize = 21;
const R_ARM_D: usize = 22;

// Group 0: indices into `Work::bufs0`, in the shader's binding order from 1.
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
// Group 1: indices into `Work::bufs1`, binding order from 1.
const TAB_C: usize = 0;
const TAB_D: usize = 1;
const KIDX: usize = 2;
const KV: usize = 3;
const FK: usize = 4;
const RQ: usize = 5;
const CARM: usize = 6;
const DARM: usize = 7;
const SQ: usize = 8;
const SD: usize = 9;
const TC: usize = 11;
const TD: usize = 12;
const O_GCH: usize = 13;
const O_GDP: usize = 14;

/// What a set of buffers was made for.
#[derive(PartialEq, Clone, Copy)]
struct Key {
    nm: usize,
    np: usize,
    nch: usize,
    cap: usize,
    nrc: usize,
    nk: usize,
    n: [usize; 3],
    /// The box and splitting width the wave vectors were made for.
    ewald: [u64; 4],
}

struct Work {
    key: Key,
    params: wgpu::Buffer,
    rparams: wgpu::Buffer,
    bufs0: Vec<wgpu::Buffer>,
    bufs1: Vec<wgpu::Buffer>,
    bind0: wgpu::BindGroup,
    bind1: wgpu::BindGroup,
    staging: wgpu::Buffer,
}

pub struct GpuEwaldStep {
    device: wgpu::Device,
    queue: wgpu::Queue,
    layout0: wgpu::BindGroupLayout,
    layout1: wgpu::BindGroupLayout,
    pipelines: Vec<wgpu::ComputePipeline>,
    pub name: String,
    work: Mutex<Option<Work>>,
    /// Seconds on the card (submit to readback) and the number of calls, for timing.
    pub timing: Mutex<(f64, f64, usize)>,
}

fn layout_entries(storage: u32) -> Vec<wgpu::BindGroupLayoutEntry> {
    let mut entries = vec![wgpu::BindGroupLayoutEntry { binding: 0, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None }, count: None }];
    for b in 1..=storage {
        entries.push(wgpu::BindGroupLayoutEntry { binding: b, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: false }, has_dynamic_offset: false, min_binding_size: None }, count: None });
    }
    entries
}

impl GpuEwaldStep {
    /// The first GPU wgpu finds that offers double precision in shaders, or why not.
    pub fn new() -> Result<GpuEwaldStep, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions { power_preference: wgpu::PowerPreference::HighPerformance, force_fallback_adapter: false, compatible_surface: None, ..Default::default() })).map_err(|e| format!("no adapter: {e}"))?;
        if !adapter.features().contains(wgpu::Features::SHADER_F64) {
            return Err("the adapter has no double precision in shaders".into());
        }
        let info = adapter.get_info();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor { label: Some("phys-gpu-ewald-step"), required_features: wgpu::Features::SHADER_F64, required_limits: adapter.limits(), ..Default::default() })).map_err(|e| format!("no device: {e}"))?;
        let layout0 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor { label: Some("step0"), entries: &layout_entries(31) });
        let layout1 = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor { label: Some("step1"), entries: &layout_entries(15) });
        let pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor { label: Some("step"), bind_group_layouts: &[Some(&layout0), Some(&layout1)], immediate_size: 0 });
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("ewald_step"), source: wgpu::ShaderSource::Wgsl(include_str!("ewald_step.wgsl").into()) });
        let pipelines = NAMES
            .iter()
            .map(|n| device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: Some(n), layout: Some(&pl), module: &module, entry_point: Some(n), compilation_options: Default::default(), cache: None }))
            .collect();
        Ok(GpuEwaldStep { device, queue, layout0, layout1, pipelines, name: format!("gpu ewald step (f64, {} via {:?})", info.name, info.backend), work: Mutex::new(None), timing: Mutex::new((0.0, 0.0, 0)) })
    }

    fn make_work(&self, key: Key) -> Work {
        let (nm, np, nch, cap, nrc, nk) = (key.nm, key.np, key.nch, key.cap, key.nrc, key.nk);
        let (n, nc) = (nm * np, nm * nch);
        let npairs = np * (np - 1) / 2;
        let dd_n = cap * np * np;
        let cd_n = cap * 2 * np * nch;
        let in_n = nm * npairs;
        let sizes0: [usize; 31] = [
            3 * n, n, n, 3 * nc, nc, nc, 8 * cap, cap, nm / 2 + cap + 2, 4 * dd_n, 3 * dd_n, 2 * cd_n, 3 * cd_n, 4 * in_n, 3 * in_n, 3 * n, 3 * n, 3 * n, 3 * n, 3 * n, 3 * n, 16, 3 * dd_n, 3 * cd_n, 3 * n, 3 * nc, 3 * cap, cap,
            2 * nm * np * nch, 3 * nm * np * nch, 3 * n,
        ];
        let tab = |ns: usize| 2 * ns * (key.n[0] + key.n[1] + key.n[2]);
        let sizes1: [usize; 15] = [tab(nrc), tab(n), (3 * nk).div_ceil(2), 3 * nk, nk, nrc, 3 * nrc, 3 * n, 2 * nk, 2 * nk, 2 * nk, 2 * nk, 2 * nk, 3 * nrc, 3 * n];
        let usage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST;
        let make = |sizes: &[usize]| -> Vec<wgpu::Buffer> { sizes.iter().map(|&s| self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("step"), size: (s.max(1) * 8) as u64, usage, mapped_at_creation: false })).collect() };
        let (bufs0, bufs1) = (make(&sizes0), make(&sizes1));
        let uniform = |size: usize| self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("params"), size: size as u64, usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let (params, rparams) = (uniform(std::mem::size_of::<Params>()), uniform(std::mem::size_of::<RParams>()));
        let group = |layout: &wgpu::BindGroupLayout, uni: &wgpu::Buffer, bufs: &[wgpu::Buffer]| {
            let mut entries = vec![wgpu::BindGroupEntry { binding: 0, resource: uni.as_entire_binding() }];
            for (k, b) in bufs.iter().enumerate() {
                entries.push(wgpu::BindGroupEntry { binding: 1 + k as u32, resource: b.as_entire_binding() });
            }
            self.device.create_bind_group(&wgpu::BindGroupDescriptor { label: Some("step"), layout, entries: &entries })
        };
        let bind0 = group(&self.layout0, &params, &bufs0);
        let bind1 = group(&self.layout1, &rparams, &bufs1);
        // Everything one call reads back: the solve's (mu, e0, fpol, fch, lforce, lvir, sc, sq),
        // then the reciprocal forces' (sd, tc, td, o_gch, o_gdp).
        let len = 3 * n * 3 + 3 * nc + 3 * cap + cap + 16 + 2 * nk + 3 * 2 * nk + 3 * nrc + 3 * n;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("staging"), size: (len * 8) as u64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        Work { key, params, rparams, bufs0, bufs1, bind0, bind1, staging }
    }

    /// Submit kernels, in order, each over the given number of workgroups, with the buffers to read
    /// back (group, index, count of doubles); returns the whole staging buffer and where each landed.
    fn submit(&self, w: &Work, kernels: &[(usize, usize)], reads: &[(u8, usize, usize)]) -> (Vec<f64>, Vec<usize>, f64) {
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_bind_group(0, &w.bind0, &[]);
            pass.set_bind_group(1, &w.bind1, &[]);
            for &(which, count) in kernels {
                pass.set_pipeline(&self.pipelines[which]);
                pass.dispatch_workgroups(count.max(1) as u32, 1, 1);
            }
        }
        let mut off = 0usize;
        let mut at = Vec::new();
        for &(g, buf, count) in reads {
            if count > 0 {
                let src = if g == 0 { &w.bufs0[buf] } else { &w.bufs1[buf] };
                enc.copy_buffer_to_buffer(src, 0, &w.staging, (off * 8) as u64, (count * 8) as u64);
            }
            at.push(off);
            off += count;
        }
        let t = std::time::Instant::now();
        self.queue.submit(Some(enc.finish()));
        let slice = w.staging.slice(..);
        slice.map_async(wgpu::MapMode::Read, |r| r.expect("a result mapped"));
        self.device.poll(wgpu::PollType::wait_indefinitely()).expect("the GPU finished");
        let secs = t.elapsed().as_secs_f64();
        let data: Vec<f64> = bytemuck::cast_slice(&slice.get_mapped_range().expect("the results range")).to_vec();
        w.staging.unmap();
        (data, at, secs)
    }

    #[allow(clippy::too_many_arguments)]
    fn run(&self, ew: &Ewald, cl: &Cluster, links: &[Link], cpos: &[Vec3], cq: &[f64], carm: &[Vec3], dpos: &[Vec3], darm: &[Vec3], warm: Option<&[Vec<Vec3>]>, tolerance: f64) -> Option<EwaldStep> {
        let t_all = std::time::Instant::now();
        let nm = cl.pol.len();
        if nm == 0 || links.is_empty() || cpos.is_empty() {
            return None;
        }
        let (np, nch) = (cl.pol[0].len(), cl.charges[0].len());
        if np == 0 || cl.pol.iter().any(|m| m.len() != np) || cl.charges.iter().any(|m| m.len() != nch) {
            return None;
        }
        let (n, nc, nl) = (nm * np, nm * nch, links.len());
        let (nrc, nk) = (cpos.len(), ew.nk());
        if dpos.len() != n || nk == 0 {
            return None;
        }
        let npairs = np * (np - 1) / 2;
        let nax = ew.nmax().map(|m| m + 1);
        let key = Key { nm, np, nch, cap: nl + nl / 4 + 64, nrc, nk, n: nax, ewald: [ew.s.to_bits(), ew.cell[0].to_bits(), ew.cell[1].to_bits(), ew.cell[2].to_bits()] };
        let mut guard = self.work.lock().expect("the step work");
        let rebuild = match guard.as_ref() {
            Some(w) => {
                let k = w.key;
                (k.nm, k.np, k.nch, k.nrc, k.nk, k.n, k.ewald) != (nm, np, nch, nrc, nk, nax, key.ewald) || k.cap < nl
            }
            None => true,
        };
        if rebuild {
            *guard = Some(self.make_work(key));
        }
        let w = guard.as_ref().expect("the work");
        let q = &self.queue;
        let put0 = |i: usize, bytes: &[u8]| {
            if !bytes.is_empty() {
                q.write_buffer(&w.bufs0[i], 0, bytes)
            }
        };
        let put1 = |i: usize, bytes: &[u8]| {
            if !bytes.is_empty() {
                q.write_buffer(&w.bufs1[i], 0, bytes)
            }
        };
        let f3 = |it: &mut dyn Iterator<Item = Vec3>| -> Vec<f64> { it.flat_map(|p| [p.x, p.y, p.z]).collect() };
        // ---- the induced-dipole solve's data, as induction.rs uploads it
        put0(POL_POS, bytemuck::cast_slice(&f3(&mut cl.pol.iter().flat_map(|m| m.iter().map(|p| p.pos)))));
        let sigma: Vec<f64> = cl.pol.iter().flat_map(|m| m.iter().map(|p| phys::induction::sigma_of(p.alpha))).collect();
        put0(POL_SIG, bytemuck::cast_slice(&sigma));
        let inv_alpha: Vec<f64> = cl.pol.iter().flat_map(|m| m.iter().map(|p| 1.0 / p.alpha)).collect();
        put0(INV_ALPHA, bytemuck::cast_slice(&inv_alpha));
        if nc > 0 {
            put0(CH_POS, bytemuck::cast_slice(&f3(&mut cl.charges.iter().flat_map(|m| m.iter().map(|c| c.pos)))));
            let chq: Vec<f64> = cl.charges.iter().flat_map(|m| m.iter().map(|c| c.q)).collect();
            put0(CH_Q, bytemuck::cast_slice(&chq));
            let chs: Vec<f64> = cl.charges.iter().flat_map(|m| m.iter().map(|c| c.sigma)).collect();
            put0(CH_SIG, bytemuck::cast_slice(&chs));
        }
        let mut lf = Vec::with_capacity(nl * 8);
        let mut lij: Vec<u32> = Vec::with_capacity(nl * 2);
        for l in links {
            lf.extend_from_slice(&[l.shift.x, l.shift.y, l.shift.z, l.weight, l.dweight, l.d.x, l.d.y, l.d.z]);
            lij.extend_from_slice(&[l.i as u32, l.j as u32]);
        }
        put0(LF, bytemuck::cast_slice(&lf));
        put0(LIJ, bytemuck::cast_slice(&lij));
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
        put0(INC, bytemuck::cast_slice(&inc));
        let mu0: Vec<f64> = match warm {
            Some(wm) if wm.iter().map(|m| m.len()).sum::<usize>() == n => wm.iter().flatten().flat_map(|p| [p.x, p.y, p.z]).collect(),
            _ => vec![0.0; 3 * n],
        };
        put0(MU, bytemuck::cast_slice(&mu0));
        // [tolerance, ., ., ., ., ., the Ewald width, 8 pi / V, a dipole's own coupling]
        let mut sc0 = vec![0.0f64; 16];
        sc0[0] = tolerance;
        sc0[6] = ew.s;
        sc0[7] = 8.0 * PI / ew.volume;
        sc0[8] = 4.0 / (3.0 * PI.sqrt() * ew.s.powi(3));
        put0(SC, bytemuck::cast_slice(&sc0));
        let params = Params { nm: nm as u32, np: np as u32, nch: nch as u32, nl: nl as u32, n: n as u32, nc: nc as u32, npairs: npairs as u32, max_iterations: MAX_ITERATIONS };
        q.write_buffer(&w.params, 0, bytemuck::bytes_of(&params));
        // ---- the reciprocal sums' data, as ewald.rs uploads it
        let rparams = RParams { nc: nrc as u32, nk: nk as u32, n0: nax[0] as u32, n1: nax[1] as u32, n2: nax[2] as u32, pad: [0; 3] };
        q.write_buffer(&w.rparams, 0, bytemuck::bytes_of(&rparams));
        if rebuild {
            let (idx, kv, fk) = ew.k_data();
            let kidx: Vec<i32> = idx.iter().flat_map(|a| [a[0], a[1], a[2]]).collect();
            put1(KIDX, bytemuck::cast_slice(&kidx));
            let kvf: Vec<f64> = kv.iter().flat_map(|a| [a[0], a[1], a[2]]).collect();
            put1(KV, bytemuck::cast_slice(&kvf));
            put1(FK, bytemuck::cast_slice(fk));
        }
        let tables = |pos: &[Vec3]| -> Vec<f64> {
            let mut out = Vec::new();
            for a in 0..3 {
                for p in pos {
                    let x = match a {
                        0 => p.x,
                        1 => p.y,
                        _ => p.z,
                    };
                    for c in ew.axis_table(x, a) {
                        out.push(c.re);
                        out.push(c.im);
                    }
                }
            }
            out
        };
        put1(TAB_C, bytemuck::cast_slice(&tables(cpos)));
        put1(TAB_D, bytemuck::cast_slice(&tables(dpos)));
        put1(RQ, bytemuck::cast_slice(cq));
        put1(CARM, bytemuck::cast_slice(&f3(&mut carm.iter().cloned())));
        put1(DARM, bytemuck::cast_slice(&f3(&mut darm.iter().cloned())));

        let (dd_n, cd_n, icd_n) = (nl * np * np, nl * 2 * np * nch, nm * np * nch);
        let in_n = nm * npairs;
        let groups = |count: usize, size: usize| count.div_ceil(size);
        let solve_reads = [(0u8, MU, 3 * n), (0, E0, 3 * n), (0, FPOL, 3 * n), (0, FCH, 3 * nc), (0, LFORCE, 3 * nl), (0, LVIR, nl), (0, SC, 16), (1, SQ, 2 * nk)];
        let mut first = true;
        let mut t_gpu = 0.0;
        let (data, at) = loop {
            let mut kernels: Vec<(usize, usize)> = Vec::new();
            if first {
                kernels.extend_from_slice(&[
                    (K_TERMS_DD, groups(dd_n, 64)),
                    (K_TERMS_CD, groups(cd_n, 64)),
                    (K_TERMS_INTRA, groups(in_n, 64)),
                    (K_TERMS_ICD, groups(icd_n, 64)),
                    (R_SQ, nk),
                    (R_FIELD, n),
                    (K_E0, n),
                    (R_DSUM_MU, nk),
                    (K_MATVEC_MU, n),
                    (R_RESP_MU, n),
                    (K_CG_INIT, 1),
                ]);
                first = false;
            }
            for _ in 0..BATCH {
                kernels.extend_from_slice(&[(R_DSUM_P, nk), (K_MATVEC_P, n), (R_RESP_P, n), (K_CG_B, 1)]);
            }
            kernels.extend_from_slice(&[(K_LINK_FORCE, groups(nl, 64)), (K_GATHER_POL, n)]);
            if nc > 0 {
                kernels.push((K_GATHER_CHARGE, nc));
            }
            let (data, at, secs) = self.submit(w, &kernels, &solve_reads);
            t_gpu += secs;
            let sc = &data[at[6]..at[6] + 16];
            if sc[3] == 0.0 {
                continue;
            }
            break (data, at);
        };
        let sc = &data[at[6]..at[6] + 16];
        let iterations = sc[4] as usize;
        let vec3s = |off: usize, count: usize| -> Vec<Vec3> { (0..count).map(|i| Vec3 { x: data[off + 3 * i], y: data[off + 3 * i + 1], z: data[off + 3 * i + 2] }).collect() };
        let mu = vec3s(at[0], n);
        let e0 = vec3s(at[1], n);
        let energy = -0.5 * mu.iter().zip(&e0).map(|(m, e)| m.dot(*e)).sum::<f64>();
        let fp = vec3s(at[2], n);
        let fc = vec3s(at[3], nc);
        let link_force = vec3s(at[4], nl);
        let virial: f64 = data[at[5]..at[5] + nl].iter().sum();
        let sq: Vec<C> = data[at[7]..at[7] + 2 * nk].chunks(2).map(|c| C { re: c[0], im: c[1] }).collect();
        let split = |v: &[Vec3], per: usize| -> Vec<Vec<Vec3>> { (0..nm).map(|i| v[i * per..(i + 1) * per].to_vec()).collect() };
        let induced = Induced { energy, dipoles: split(&mu, np), force_pol: split(&fp, np), force_charge: split(&fc, nch), link_force, virial, iterations };
        // ---- the reciprocal forces and arm sums, with the dipoles final
        let kernels = [(R_DSUM_MU, nk), (R_STOT, groups(nk, 64)), (R_GCH, nrc), (R_GDPOS, n), (R_ARM_C, nk), (R_ARM_D, nk)];
        let reads = [(1u8, SD, 2 * nk), (1, TC, 2 * nk), (1, TD, 2 * nk), (1, O_GCH, 3 * nrc), (1, O_GDP, 3 * n)];
        let (rdata, rat, secs) = self.submit(w, &kernels, &reads);
        t_gpu += secs;
        let complex = |off: usize| -> Vec<C> { rdata[off..off + 2 * nk].chunks(2).map(|c| C { re: c[0], im: c[1] }).collect() };
        let v3 = |off: usize, count: usize| -> Vec<Vec3> { (0..count).map(|i| Vec3 { x: rdata[off + 3 * i], y: rdata[off + 3 * i + 1], z: rdata[off + 3 * i + 2] }).collect() };
        let finished = Finished { d: complex(rat[0]), tc: complex(rat[1]), td: complex(rat[2]), grad_charge: v3(rat[3], nrc), grad_dipole_pos: v3(rat[4], n) };
        let mut tm = self.timing.lock().expect("the timing");
        tm.0 += t_gpu;
        tm.1 += t_all.elapsed().as_secs_f64();
        tm.2 += 1;
        Some(EwaldStep { sq, induced, finished })
    }
}

impl EwaldStepEngine for GpuEwaldStep {
    fn step(&self, ew: &Ewald, cl: &Cluster, links: &[Link], cpos: &[Vec3], cq: &[f64], carm: &[Vec3], dpos: &[Vec3], darm: &[Vec3], warm: Option<&[Vec<Vec3>]>, tolerance: f64) -> Option<EwaldStep> {
        self.run(ew, cl, links, cpos, cq, carm, dpos, darm, warm, tolerance)
    }
}

/// Install the GPU's whole Ewald step for this process, returning its name; or why not, in which
/// case the CPU's stays.
pub fn install_ewald_step() -> Result<String, String> {
    let gpu = GpuEwaldStep::new()?;
    let name = gpu.name.clone();
    set_step_engine(Some(Arc::new(gpu)));
    Ok(name)
}
