//! Hartree-Fock and MP2's heavy passes on the GPU, in double precision, against
//! whitened three-index integrals kept resident — `docs/PLAY.md` Phase 6, E8a.
//! This implements the core crate's [`phys::electrons::hf::FockEngine`]; a
//! program that wants it calls [`install_fock`].
//!
//! # What is resident, and why that is the point
//!
//! A Hartree-Fock iteration needs the integrals `(mn|P)` contracted with the
//! occupied orbitals, and the CPU reads the whole table to do it, whitens what
//! it read (a forward substitution for every entry, `nk^2/2` operations each),
//! and does both again next iteration. Here the whitening is done once, when a
//! table is first used, and `B[mn][k]` (for `m >= n`, `nk` doubles a pair) lives
//! in the card's memory from then on. An iteration is then
//!
//! * the contraction `x[i][m][k] = sum_n B[mn][k] c[n][i]` (one pass over `B`),
//! * the exchange matrix `K = sum_ik x x^T` (a product, lower triangle of
//!   tiles only),
//! * and the Coulomb matrix: the density fitted onto the directions `dt[k] =
//!   sum_p D[p] B[p][k]` and `J[p] = sum_k B[p][k] dt[k]` (two more passes),
//!
//! with energy `1/2 dt . dt`. The same table serves every field of a
//! counterpoise correction, a cluster, a finite-field polarisability.
//!
//! # Precision, size and limits
//!
//! All double: the answers are the CPU's to round-off. `B` is split into
//! segments of at most 1.5 GB (a storage binding is limited to just under 2
//! GB), each pass run once a segment; if the card cannot hold the whole of it
//! `load` says so and the CPU does what it always did. A dispatch is kept to a
//! few billion multiply-adds so that Windows does not reset the display driver.
//! More than 32 occupied orbitals are done in blocks of 32.

use phys::electrons::hf::FockEngine;
use phys::electrons::linalg::Matrix;
use std::sync::Mutex;
use wgpu::util::DeviceExt;

/// Multiply-adds one dispatch does, at most.
const WORK_PER_DISPATCH: f64 = 3e9;
/// The most bytes of `B` a segment holds.
const SEGMENT_BYTES: usize = 1_500_000_000;
/// Orbitals the contraction shader holds at a time.
const BLOCK: usize = 32;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    n: u32,
    nk: u32,
    count: u32,
    p0: u32,
    p1: u32,
    m0: u32,
    rows: u32,
    pad: u32,
}

struct Resident {
    key: usize,
    n: usize,
    nk: usize,
    npairs: usize,
    seg_pairs: usize,
    segs: Vec<wgpu::Buffer>,
}

/// The GPU, its shaders, and the integrals currently resident.
pub struct GpuFock {
    device: wgpu::Device,
    queue: wgpu::Queue,
    contract: wgpu::ComputePipeline,
    kernel_k: wgpu::ComputePipeline,
    dt: wgpu::ComputePipeline,
    jp: wgpu::ComputePipeline,
    name: String,
    resident: Mutex<Option<Resident>>,
}

fn pidx(a: usize, b: usize) -> usize {
    if a >= b {
        a * (a + 1) / 2 + b
    } else {
        b * (b + 1) / 2 + a
    }
}

impl GpuFock {
    /// The first GPU wgpu finds that offers double precision in shaders, or why not.
    pub fn new() -> Result<GpuFock, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions { power_preference: wgpu::PowerPreference::HighPerformance, force_fallback_adapter: false, compatible_surface: None, ..Default::default() })).map_err(|e| format!("no adapter: {e}"))?;
        if !adapter.features().contains(wgpu::Features::SHADER_F64) {
            return Err("the adapter has no double precision in shaders".into());
        }
        let info = adapter.get_info();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor { label: Some("phys-gpu-fock"), required_features: wgpu::Features::SHADER_F64, required_limits: adapter.limits(), ..Default::default() })).map_err(|e| format!("no device: {e}"))?;
        let pipeline = |label: &str, src: &str| {
            let module = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some(label), source: wgpu::ShaderSource::Wgsl(src.to_string().into()) });
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: Some(label), layout: None, module: &module, entry_point: Some("main"), compilation_options: Default::default(), cache: None })
        };
        let contract = pipeline("contract", include_str!("fock_contract.wgsl"));
        let kernel_k = pipeline("exchange", include_str!("fock_k.wgsl"));
        let dt = pipeline("fit", include_str!("fock_dt.wgsl"));
        let jp = pipeline("coulomb", include_str!("fock_j.wgsl"));
        Ok(GpuFock { device, queue, contract, kernel_k, dt, jp, name: format!("gpu Fock (f64, {} via {:?})", info.name, info.backend), resident: Mutex::new(None) })
    }

    /// The adapter's name, for a message.
    pub fn name_string(&self) -> String {
        self.name.clone()
    }

    fn storage(&self, label: &str, bytes: usize) -> wgpu::Buffer {
        self.device.create_buffer(&wgpu::BufferDescriptor { label: Some(label), size: bytes.max(8) as u64, usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false })
    }

    fn upload(&self, label: &str, values: &[f64]) -> wgpu::Buffer {
        self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some(label), contents: bytemuck::cast_slice(values), usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST })
    }

    /// One dispatch, run to completion before the next.
    fn run(&self, pipeline: &wgpu::ComputePipeline, params: &Params, buffers: &[&wgpu::Buffer], groups: (u32, u32)) {
        let dev = &self.device;
        let p = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("params"), contents: bytemuck::bytes_of(params), usage: wgpu::BufferUsages::UNIFORM });
        let layout = pipeline.get_bind_group_layout(0);
        let mut entries = vec![wgpu::BindGroupEntry { binding: 0, resource: p.as_entire_binding() }];
        for (k, b) in buffers.iter().enumerate() {
            entries.push(wgpu::BindGroupEntry { binding: k as u32 + 1, resource: b.as_entire_binding() });
        }
        let bind = dev.create_bind_group(&wgpu::BindGroupDescriptor { label: Some("fock"), layout: &layout, entries: &entries });
        let mut enc = dev.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("fock") });
        {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: Some("fock"), timestamp_writes: None });
            pass.set_pipeline(pipeline);
            pass.set_bind_group(0, &bind, &[]);
            pass.dispatch_workgroups(groups.0, groups.1, 1);
        }
        self.queue.submit(Some(enc.finish()));
        dev.poll(wgpu::PollType::wait_indefinitely()).expect("the GPU finished a pass");
    }

    fn read(&self, buf: &wgpu::Buffer, len: usize) -> Vec<f64> {
        let dev = &self.device;
        let bytes = (len * 8) as u64;
        let readback = dev.create_buffer(&wgpu::BufferDescriptor { label: Some("readback"), size: bytes, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let mut enc = dev.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("readback") });
        enc.copy_buffer_to_buffer(buf, 0, &readback, 0, bytes);
        self.queue.submit(Some(enc.finish()));
        let slice = readback.slice(..);
        slice.map_async(wgpu::MapMode::Read, |r| r.expect("a result mapped"));
        dev.poll(wgpu::PollType::wait_indefinitely()).expect("the results arrived");
        let out: Vec<f64> = bytemuck::cast_slice(&slice.get_mapped_range().expect("the results' range")).to_vec();
        readback.unmap();
        out
    }

    /// `x[i][m][k] = sum_n B(m, n)[k] c[n][i]` for the `count` columns of `c`
    /// (`n x count`, row-major), on a fresh buffer, which is returned.
    fn contract_all(&self, r: &Resident, c: &[f64], count: usize) -> wgpu::Buffer {
        let (n, nk) = (r.n, r.nk);
        let cb = self.upload("c", c);
        let x = self.storage("x", count * n * nk * 8);
        let rows_per = ((WORK_PER_DISPATCH / (n as f64 * nk as f64 * count as f64)) as usize).clamp(1, n);
        for (s, seg) in r.segs.iter().enumerate() {
            let (p0, p1) = (s * r.seg_pairs, ((s + 1) * r.seg_pairs).min(r.npairs));
            let mut m0 = 0;
            while m0 < n {
                let rows = rows_per.min(n - m0);
                let params = Params { n: n as u32, nk: nk as u32, count: count as u32, p0: p0 as u32, p1: p1 as u32, m0: m0 as u32, rows: rows as u32, pad: 0 };
                self.run(&self.contract, &params, &[seg, &cb, &x], (nk.div_ceil(64) as u32, rows as u32));
                m0 += rows;
            }
        }
        x
    }
}

/// Install the GPU's Fock engine for this process, returning its name; or why
/// not, in which case the CPU's path stays.
pub fn install_fock() -> Result<String, String> {
    let gpu = GpuFock::new()?;
    let name = gpu.name.clone();
    phys::electrons::hf::set_fock_engine(Some(std::sync::Arc::new(gpu)));
    Ok(name)
}

impl FockEngine for GpuFock {
    fn load(&self, key: usize, n: usize, nk: usize, source: &dyn Fn(&mut dyn FnMut(usize, usize, &[f64]))) -> bool {
        *self.resident.lock().expect("the resident integrals") = None;
        let npairs = n * (n + 1) / 2;
        let seg_pairs = (SEGMENT_BYTES / (nk * 8)).max(1);
        let nseg = npairs.div_ceil(seg_pairs);
        let scope = self.device.push_error_scope(wgpu::ErrorFilter::OutOfMemory);
        let mut segs = Vec::new();
        for s in 0..nseg {
            let pairs = seg_pairs.min(npairs - s * seg_pairs);
            segs.push(self.storage("B", pairs * nk * 8));
        }
        if pollster::block_on(scope.pop()).is_some() {
            return false;
        }
        let mut written = 0usize;
        source(&mut |m, nn, v| {
            let p = pidx(m, nn);
            let (s, off) = (p / seg_pairs, (p % seg_pairs) * nk * 8);
            self.queue.write_buffer(&segs[s], off as u64, bytemuck::cast_slice(v));
            written += 1;
            if written % 256 == 0 {
                self.queue.submit(std::iter::empty());
            }
        });
        self.queue.submit(std::iter::empty());
        self.device.poll(wgpu::PollType::wait_indefinitely()).expect("the integrals arrived");
        // Pairs of functions whose product is negligible are not in the table;
        // their integrals are zero, which is what a new buffer holds.
        assert!(written <= npairs, "the table gave {written} entries for {npairs} pairs");
        *self.resident.lock().expect("the resident integrals") = Some(Resident { key, n, nk, npairs, seg_pairs, segs });
        true
    }

    fn is_loaded(&self, key: usize) -> bool {
        self.resident.lock().expect("the resident integrals").as_ref().map(|r| r.key == key).unwrap_or(false)
    }

    fn fock(&self, key: usize, c: &[f64], stride: usize, count: usize, d: &Matrix) -> (Matrix, f64, Matrix) {
        let guard = self.resident.lock().expect("the resident integrals");
        let r = guard.as_ref().filter(|r| r.key == key).expect("integrals resident for this table");
        let (n, nk, npairs) = (r.n, r.nk, r.npairs);
        // The density, packed by pair and doubled off the diagonal.
        let mut dp = vec![0.0f64; npairs];
        for m in 0..n {
            for nn in 0..=m {
                dp[pidx(m, nn)] = if m == nn { d.get(m, m) } else { d.get(m, nn) + d.get(nn, m) };
            }
        }
        let dpb = self.upload("dp", &dp);
        let dtb = self.storage("dt", nk * 8);
        let jpb = self.storage("jp", npairs * 8);
        let kmat = self.storage("kmat", n * n * 8);
        // Coulomb: the density fitted, then the matrix.
        for (s, seg) in r.segs.iter().enumerate() {
            let (p0, p1) = (s * r.seg_pairs, ((s + 1) * r.seg_pairs).min(npairs));
            let params = Params { n: n as u32, nk: nk as u32, count: 0, p0: p0 as u32, p1: p1 as u32, m0: 0, rows: 0, pad: 0 };
            self.run(&self.dt, &params, &[seg, &dpb, &dtb], (nk.div_ceil(64) as u32, 1));
        }
        for (s, seg) in r.segs.iter().enumerate() {
            let (p0, p1) = (s * r.seg_pairs, ((s + 1) * r.seg_pairs).min(npairs));
            let pairs = p1 - p0;
            let gx = pairs.min(65535);
            let params = Params { n: n as u32, nk: nk as u32, count: 0, p0: p0 as u32, p1: p1 as u32, m0: 0, rows: 0, pad: 0 };
            self.run(&self.jp, &params, &[seg, &dtb, &jpb], (gx as u32, pairs.div_ceil(gx) as u32));
        }
        // Exchange, the orbitals a block at a time.
        let tiles = n.div_ceil(64);
        for first in (0..count).step_by(BLOCK) {
            let cnt = BLOCK.min(count - first);
            let mut cb = vec![0.0f64; n * cnt];
            for a in 0..n {
                for i in 0..cnt {
                    cb[a * cnt + i] = c[a * stride + first + i];
                }
            }
            let x = self.contract_all(r, &cb, cnt);
            let mut r0 = 0;
            while r0 < tiles {
                let mut r1 = r0;
                let mut work = 0.0;
                while r1 < tiles && (work == 0.0 || work + (r1 + 1) as f64 * 4096.0 * nk as f64 * cnt as f64 <= WORK_PER_DISPATCH) {
                    work += (r1 + 1) as f64 * 4096.0 * nk as f64 * cnt as f64;
                    r1 += 1;
                }
                let params = Params { n: n as u32, nk: nk as u32, count: cnt as u32, p0: 0, p1: 0, m0: r0 as u32, rows: (r1 - r0) as u32, pad: 0 };
                self.run(&self.kernel_k, &params, &[&x, &kmat], (tiles as u32, (r1 - r0) as u32));
                r0 = r1;
            }
        }
        let dt = self.read(&dtb, nk);
        let jp = self.read(&jpb, npairs);
        let kraw = self.read(&kmat, n * n);
        let ej = 0.5 * dt.iter().map(|v| v * v).sum::<f64>();
        let mut j = Matrix::zeros(n);
        let mut k = Matrix::zeros(n);
        for m in 0..n {
            for nn in 0..=m {
                let v = jp[pidx(m, nn)];
                j.set(m, nn, v);
                j.set(nn, m, v);
                let kv = kraw[m * n + nn];
                k.set(m, nn, kv);
                k.set(nn, m, kv);
            }
        }
        (j, ej, k)
    }

    fn half(&self, key: usize, c: &[f64], stride: usize, first: usize, count: usize) -> Vec<f64> {
        let guard = self.resident.lock().expect("the resident integrals");
        let r = guard.as_ref().filter(|r| r.key == key).expect("integrals resident for this table");
        let (n, nk) = (r.n, r.nk);
        let mut out = Vec::with_capacity(count * n * nk);
        for f in (0..count).step_by(BLOCK) {
            let cnt = BLOCK.min(count - f);
            let mut cb = vec![0.0f64; n * cnt];
            for a in 0..n {
                for i in 0..cnt {
                    cb[a * cnt + i] = c[a * stride + first + f + i];
                }
            }
            let x = self.contract_all(r, &cb, cnt);
            out.extend(self.read(&x, cnt * n * nk));
        }
        out
    }

    fn name(&self) -> &str {
        &self.name
    }
}
