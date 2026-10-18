//! The reciprocal-space sums of `phys::ewald` on the GPU, in double precision —
//! `docs/PLAY.md` E8e, the owner's decision of 2026-10-18 (the Ewald port).
//!
//! In a force evaluation of the 216-molecule water liquid the reciprocal sums
//! were about half of 93 ms (`docs/SOLUTIONS.md` 9.1): the phases and structure
//! factor of the charges, the dipoles' phases, and the dipole response that the
//! conjugate-gradient solve asks for eleven times. They are dense and regular
//! (sites x 1535 wave vectors), which is what a card is for.
//!
//! **What moves and what stays.** The CPU builds each site's three axis tables of
//! `e^(i m theta)` (cos and sin, which WGSL has no `f64` version of) with
//! `Ewald::axis_table`, the very numbers its own phases are products of; the
//! card forms `e^(ik.r) = (tx ty) tz` from them and does the sums. A thread a wave
//! vector adds over the sites in index order, and a thread a site adds over the
//! wave vectors in index order, as the reference's loops do, so the only
//! difference from the reference is the card's arithmetic (it may fuse a multiply
//! and an add). The reference's scalar sums over `k` (the energy, the virial's
//! loop) stay on the CPU, which is given the per-`k` arrays back, so they are the
//! reference's own loops.
//!
//! **It answers to the reference's tolerance, not bit for bit**, like the
//! induction engine; `gpu/tests/ewald.rs` bounds the difference and checks that
//! two calls give the same bits.

use phys::ewald::{set_reciprocal_engine, Ewald, Finished, ReciprocalEngine, ReciprocalSession, C};
use phys::math::Vec3;
use std::sync::{Arc, Mutex};

const PI: f64 = std::f64::consts::PI;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    nc: u32,
    nd: u32,
    nk: u32,
    n0: u32,
    n1: u32,
    n2: u32,
    pad: [u32; 2],
}

const K_SQ: usize = 0;
const K_DSUM: usize = 1;
const K_STOT: usize = 2;
const K_FIELD: usize = 3;
const K_RESP: usize = 4;
const K_GCH: usize = 5;
const K_GDPOS: usize = 6;
const K_ARM_C: usize = 7;
const K_ARM_D: usize = 8;
const NAMES: [&str; 9] = ["k_sq", "k_dsum", "k_stot", "k_field", "k_resp", "k_gch", "k_gdpos", "k_arm_c", "k_arm_d"];

// Indices into `Work::bufs`, in the shader's binding order from 1.
const TAB_C: usize = 0;
const TAB_D: usize = 1;
const KIDX: usize = 2;
const KV: usize = 3;
const FK: usize = 4;
const Q: usize = 5;
const CARM: usize = 6;
const MU: usize = 7;
const DARM: usize = 8;
const SQ: usize = 9;
const SD: usize = 10;
const STOT: usize = 11;
const TC: usize = 12;
const TD: usize = 13;
const O_D: usize = 14;
const O_GCH: usize = 15;
const O_GDP: usize = 16;
const SC: usize = 17;

/// What a set of buffers was made for, and what is already in them that does not change.
#[derive(PartialEq, Clone, Copy)]
struct Key {
    nc: usize,
    nd: usize,
    nk: usize,
    n: [usize; 3],
    /// The box and the splitting width the wave vectors were made for.
    ewald: [u64; 4],
}

struct Work {
    key: Key,
    params: wgpu::Buffer,
    bufs: Vec<wgpu::Buffer>,
    bind: wgpu::BindGroup,
    staging: wgpu::Buffer,
}

struct Inner {
    device: wgpu::Device,
    queue: wgpu::Queue,
    layout: wgpu::BindGroupLayout,
    pipelines: Vec<wgpu::ComputePipeline>,
    work: Mutex<Option<Work>>,
    /// Seconds on the card (submit to readback) and the number of submissions, for timing.
    timing: Mutex<(f64, usize)>,
}

/// The reciprocal engine on the card.
pub struct GpuEwald {
    inner: Arc<Inner>,
    pub name: String,
}

impl GpuEwald {
    /// The first GPU wgpu finds that offers double precision in shaders, or why not.
    pub fn new() -> Result<GpuEwald, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions { power_preference: wgpu::PowerPreference::HighPerformance, force_fallback_adapter: false, compatible_surface: None, ..Default::default() })).map_err(|e| format!("no adapter: {e}"))?;
        if !adapter.features().contains(wgpu::Features::SHADER_F64) {
            return Err("the adapter has no double precision in shaders".into());
        }
        let info = adapter.get_info();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor { label: Some("phys-gpu-ewald"), required_features: wgpu::Features::SHADER_F64, required_limits: adapter.limits(), ..Default::default() })).map_err(|e| format!("no device: {e}"))?;
        let mut entries = vec![wgpu::BindGroupLayoutEntry { binding: 0, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None }, count: None }];
        for b in 1..=18u32 {
            entries.push(wgpu::BindGroupLayoutEntry { binding: b, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: false }, has_dynamic_offset: false, min_binding_size: None }, count: None });
        }
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor { label: Some("ewald"), entries: &entries });
        let pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor { label: Some("ewald"), bind_group_layouts: &[Some(&layout)], immediate_size: 0 });
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("ewald"), source: wgpu::ShaderSource::Wgsl(include_str!("ewald.wgsl").into()) });
        let pipelines = NAMES
            .iter()
            .map(|n| device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: Some(n), layout: Some(&pl), module: &module, entry_point: Some(n), compilation_options: Default::default(), cache: None }))
            .collect();
        let inner = Inner { device, queue, layout, pipelines, work: Mutex::new(None), timing: Mutex::new((0.0, 0)) };
        Ok(GpuEwald { inner: Arc::new(inner), name: format!("gpu ewald (f64, {} via {:?})", info.name, info.backend) })
    }

    /// For timing: the milliseconds a submission of these kernels (by name, each over the given
    /// number of workgroups) and a read of 8 bytes takes, on average over `reps`, with the buffers
    /// of the last session. Panics if there has been none.
    pub fn bench(&self, kernels: &[(&str, usize)], reps: usize) -> f64 {
        let guard = self.inner.work.lock().expect("the ewald work");
        let w = guard.as_ref().expect("a session first");
        let ks: Vec<(usize, usize)> = kernels.iter().map(|(n, c)| (NAMES.iter().position(|x| x == n).expect("a kernel"), *c)).collect();
        let t = std::time::Instant::now();
        for _ in 0..reps {
            self.inner.run(w, &ks, &[(SC, 1)]);
        }
        t.elapsed().as_secs_f64() * 1000.0 / reps as f64
    }

    /// Seconds on the card in all, and the number of submissions.
    pub fn timing(&self) -> (f64, usize) {
        *self.inner.timing.lock().expect("the timing")
    }
}

impl Inner {
    fn make_work(&self, key: Key) -> Work {
        let (nc, nd, nk) = (key.nc, key.nd, key.nk);
        let tab = |ns: usize| 2 * ns * (key.n[0] + key.n[1] + key.n[2]);
        // Eight-byte units, in the shader's binding order from 1 (the wave vectors' integers are
        // i32, three to a vector: half as many eight-byte units, rounded up).
        let sizes: [usize; 18] = [tab(nc), tab(nd), (3 * nk).div_ceil(2), 3 * nk, nk, nc, 3 * nc, 3 * nd, 3 * nd, 2 * nk, 2 * nk, 2 * nk, 2 * nk, 2 * nk, 3 * nd, 3 * nc, 3 * nd, 2];
        let usage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST;
        let bufs: Vec<wgpu::Buffer> = sizes.iter().map(|&s| self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("ewald"), size: (s.max(1) * 8) as u64, usage, mapped_at_creation: false })).collect();
        let params = self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("params"), size: std::mem::size_of::<Params>() as u64, usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let mut entries = vec![wgpu::BindGroupEntry { binding: 0, resource: params.as_entire_binding() }];
        for (k, b) in bufs.iter().enumerate() {
            entries.push(wgpu::BindGroupEntry { binding: 1 + k as u32, resource: b.as_entire_binding() });
        }
        let bind = self.device.create_bind_group(&wgpu::BindGroupDescriptor { label: Some("ewald"), layout: &self.layout, entries: &entries });
        // Room for every read one call makes: two structure factors and the arm sums
        // (4 nk complex), and the per-site vectors.
        let staging_len = 4 * 2 * nk + 3 * nd + 3 * nc + 3 * nd;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("staging"), size: (staging_len.max(1) * 8) as u64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        Work { key, params, bufs, bind, staging }
    }

    fn put(&self, w: &Work, i: usize, bytes: &[u8]) {
        if !bytes.is_empty() {
            self.queue.write_buffer(&w.bufs[i], 0, bytes);
        }
    }

    /// Run kernels, in order, in one submission, each over the given number of workgroups, and bring
    /// back the given buffers (the first `count` doubles of each).
    fn run(&self, w: &Work, kernels: &[(usize, usize)], reads: &[(usize, usize)]) -> Vec<Vec<f64>> {
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_bind_group(0, &w.bind, &[]);
            for &(which, count) in kernels {
                pass.set_pipeline(&self.pipelines[which]);
                pass.dispatch_workgroups(count.max(1) as u32, 1, 1);
            }
        }
        let mut off = 0u64;
        let mut at = Vec::new();
        for &(buf, count) in reads {
            if count > 0 {
                enc.copy_buffer_to_buffer(&w.bufs[buf], 0, &w.staging, off, (count * 8) as u64);
            }
            at.push(off as usize / 8);
            off += (count * 8) as u64;
        }
        let t = std::time::Instant::now();
        self.queue.submit(Some(enc.finish()));
        let slice = w.staging.slice(..);
        slice.map_async(wgpu::MapMode::Read, |r| r.expect("a result mapped"));
        self.device.poll(wgpu::PollType::wait_indefinitely()).expect("the GPU finished");
        let data: Vec<f64> = bytemuck::cast_slice(&slice.get_mapped_range().expect("the results range")).to_vec();
        w.staging.unmap();
        let mut tm = self.timing.lock().expect("the timing");
        tm.0 += t.elapsed().as_secs_f64();
        tm.1 += 1;
        reads.iter().zip(at).map(|(&(_, count), o)| data[o..o + count].to_vec()).collect()
    }
}

fn complex(v: &[f64]) -> Vec<C> {
    v.chunks(2).map(|c| C { re: c[0], im: c[1] }).collect()
}

fn vec3s(v: &[f64]) -> Vec<Vec3> {
    v.chunks(3).map(|c| Vec3 { x: c[0], y: c[1], z: c[2] }).collect()
}

fn flat3(v: &[Vec3]) -> Vec<f64> {
    v.iter().flat_map(|p| [p.x, p.y, p.z]).collect()
}

/// One force evaluation's session: the sites' tables are on the card.
struct Session {
    inner: Arc<Inner>,
    nc: usize,
    nd: usize,
    nk: usize,
    sq_done: bool,
}

impl ReciprocalEngine for GpuEwald {
    fn begin(&self, ew: &Ewald, cpos: &[Vec3], cq: &[f64], carm: &[Vec3], dpos: &[Vec3], darm: &[Vec3]) -> Option<Box<dyn ReciprocalSession>> {
        let (nc, nd, nk) = (cpos.len(), dpos.len(), ew.nk());
        if nc == 0 || nk == 0 {
            return None;
        }
        let n = ew.nmax().map(|m| m + 1);
        let key = Key { nc, nd, nk, n, ewald: [ew.s.to_bits(), ew.cell[0].to_bits(), ew.cell[1].to_bits(), ew.cell[2].to_bits()] };
        let inner = &self.inner;
        let mut guard = inner.work.lock().expect("the ewald work");
        let fresh = guard.as_ref().map(|w| w.key != key).unwrap_or(true);
        if fresh {
            *guard = Some(inner.make_work(key));
        }
        let w = guard.as_ref().expect("the work");
        if fresh {
            let (idx, kv, fk) = ew.k_data();
            let kidx: Vec<i32> = idx.iter().flat_map(|a| [a[0], a[1], a[2]]).collect();
            inner.put(w, KIDX, bytemuck::cast_slice(&kidx));
            let kvf: Vec<f64> = kv.iter().flat_map(|a| [a[0], a[1], a[2]]).collect();
            inner.put(w, KV, bytemuck::cast_slice(&kvf));
            inner.put(w, FK, bytemuck::cast_slice(fk));
            let sc = [8.0 * PI / ew.volume, 4.0 / (3.0 * PI.sqrt() * ew.s.powi(3))];
            inner.put(w, SC, bytemuck::cast_slice(&sc));
            let params = Params { nc: nc as u32, nd: nd as u32, nk: nk as u32, n0: n[0] as u32, n1: n[1] as u32, n2: n[2] as u32, pad: [0; 2] };
            inner.queue.write_buffer(&w.params, 0, bytemuck::bytes_of(&params));
        }
        // The sites' tables: axis 0's block, then axis 1's, then axis 2's, site-major.
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
        inner.put(w, TAB_C, bytemuck::cast_slice(&tables(cpos)));
        inner.put(w, Q, bytemuck::cast_slice(cq));
        inner.put(w, CARM, bytemuck::cast_slice(&flat3(carm)));
        if nd > 0 {
            inner.put(w, TAB_D, bytemuck::cast_slice(&tables(dpos)));
            inner.put(w, DARM, bytemuck::cast_slice(&flat3(darm)));
        }
        Some(Box::new(Session { inner: self.inner.clone(), nc, nd, nk, sq_done: false }))
    }
}

impl ReciprocalSession for Session {
    fn charge_sum(&mut self) -> Vec<C> {
        let guard = self.inner.work.lock().expect("the ewald work");
        let w = guard.as_ref().expect("the work");
        let out = self.inner.run(w, &[(K_SQ, self.nk)], &[(SQ, 2 * self.nk)]);
        self.sq_done = true;
        complex(&out[0])
    }

    fn field_at_dipoles(&mut self) -> Vec<Vec3> {
        assert!(self.sq_done, "the charges' structure factor first");
        let guard = self.inner.work.lock().expect("the ewald work");
        let w = guard.as_ref().expect("the work");
        let out = self.inner.run(w, &[(K_FIELD, self.nd)], &[(O_D, 3 * self.nd)]);
        vec3s(&out[0])
    }

    fn dipole_response(&mut self, mu: &[Vec3]) -> Vec<Vec3> {
        let guard = self.inner.work.lock().expect("the ewald work");
        let w = guard.as_ref().expect("the work");
        self.inner.put(w, MU, bytemuck::cast_slice(&flat3(mu)));
        let out = self.inner.run(w, &[(K_DSUM, self.nk), (K_RESP, self.nd)], &[(O_D, 3 * self.nd)]);
        vec3s(&out[0])
    }

    fn finish(&mut self, mu: &[Vec3]) -> Finished {
        let guard = self.inner.work.lock().expect("the ewald work");
        let w = guard.as_ref().expect("the work");
        let (nc, nd, nk) = (self.nc, self.nd, self.nk);
        if nd > 0 {
            self.inner.put(w, MU, bytemuck::cast_slice(&flat3(mu)));
        }
        let mut kernels = vec![(K_DSUM, nk), (K_STOT, nk.div_ceil(64)), (K_GCH, nc), (K_ARM_C, nk)];
        if nd > 0 {
            kernels.push((K_GDPOS, nd));
            kernels.push((K_ARM_D, nk));
        }
        let reads = [(SD, if nd > 0 { 2 * nk } else { 0 }), (TC, 2 * nk), (TD, if nd > 0 { 2 * nk } else { 0 }), (O_GCH, 3 * nc), (O_GDP, 3 * nd)];
        let out = self.inner.run(w, &kernels, &reads);
        if std::env::var_os("PHYS_PROFILE").is_some() {
            let tm = self.inner.timing.lock().expect("the timing");
            eprintln!("  gpu ewald: {} submissions so far, {:.2} ms on the card each", tm.1, tm.0 * 1000.0 / tm.1 as f64);
        }
        Finished {
            d: if nd > 0 { complex(&out[0]) } else { Vec::new() },
            tc: complex(&out[1]),
            td: if nd > 0 { complex(&out[2]) } else { vec![C::default(); nk] },
            grad_charge: vec3s(&out[3]),
            grad_dipole_pos: vec3s(&out[4]),
        }
    }
}

/// Install the GPU's reciprocal sums for this process, returning its name; or why not, in which
/// case the CPU's stay.
pub fn install_ewald() -> Result<String, String> {
    let gpu = GpuEwald::new()?;
    let name = gpu.name.clone();
    set_reciprocal_engine(Some(Arc::new(gpu)));
    Ok(name)
}
