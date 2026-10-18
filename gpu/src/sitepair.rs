//! The site-site pair loop of a flexible liquid on the GPU, in double
//! precision — `docs/PLAY.md` E8c. Implements the core's
//! [`phys::liquid::PairEngine`]; `sitepair.wgsl` is the shader.
//!
//! One thread a link computes the 25 site-site terms and the forces and
//! parameter derivatives they give, one workgroup a site gathers them over the
//! links its molecule is in, in a fixed order, so the answer is the same every
//! run and agrees with the CPU loop to round-off. The buffers are kept between
//! calls. A law that is not a site-site table, or has widths for some types and
//! not others, is declined and the CPU does it.

use phys::liquid::{PairEngine, PairProblem, PairResult};
use phys::math::Vec3;
use std::sync::Mutex;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    n: u32,
    ns: u32,
    nl: u32,
    ntypes: u32,
}

const NAMES: [&str; 2] = ["k_pair_link", "k_pair_gather"];

struct Work {
    n: usize,
    ns: usize,
    ntypes: usize,
    cap: usize,
    params: wgpu::Buffer,
    bufs: Vec<wgpu::Buffer>,
    bind: wgpu::BindGroup,
    staging: wgpu::Buffer,
    offsets: Vec<usize>,
}

// Indices into `Work::bufs` (binding - 1).
const TYPES: usize = 0;
const REL: usize = 1;
const QV: usize = 2;
const SV: usize = 3;
const PAIR: usize = 4;
const OVK: usize = 5;
const BDAMP: usize = 6;
const SBASE: usize = 7;
const LK: usize = 8;
const LIJ: usize = 9;
const INC: usize = 10;
const LSCAL: usize = 17;
const OUT_FSITE: usize = 18;
const OUT_FC: usize = 19;
const OUT_GQ: usize = 20;
const OUT_GS: usize = 21;

pub struct GpuPairs {
    device: wgpu::Device,
    queue: wgpu::Queue,
    layout: wgpu::BindGroupLayout,
    pipelines: Vec<wgpu::ComputePipeline>,
    pub name: String,
    work: Mutex<Option<Work>>,
    /// Seconds on the card (submit to readback), in all, and the number of calls.
    pub timing: Mutex<(f64, f64, usize)>,
}

impl GpuPairs {
    pub fn new() -> Result<GpuPairs, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions { power_preference: wgpu::PowerPreference::HighPerformance, force_fallback_adapter: false, compatible_surface: None, ..Default::default() })).map_err(|e| format!("no adapter: {e}"))?;
        if !adapter.features().contains(wgpu::Features::SHADER_F64) {
            return Err("the adapter has no double precision in shaders".into());
        }
        let info = adapter.get_info();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor { label: Some("phys-gpu-pairs"), required_features: wgpu::Features::SHADER_F64, required_limits: adapter.limits(), ..Default::default() })).map_err(|e| format!("no device: {e}"))?;
        let mut entries = vec![wgpu::BindGroupLayoutEntry { binding: 0, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None }, count: None }];
        for b in 1..=22u32 {
            entries.push(wgpu::BindGroupLayoutEntry { binding: b, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: false }, has_dynamic_offset: false, min_binding_size: None }, count: None });
        }
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor { label: Some("pairs"), entries: &entries });
        let pl = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor { label: Some("pairs"), bind_group_layouts: &[Some(&layout)], immediate_size: 0 });
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("pairs"), source: wgpu::ShaderSource::Wgsl(include_str!("sitepair.wgsl").into()) });
        let pipelines = NAMES
            .iter()
            .map(|n| device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: Some(n), layout: Some(&pl), module: &module, entry_point: Some(n), compilation_options: Default::default(), cache: None }))
            .collect();
        Ok(GpuPairs { device, queue, layout, pipelines, name: format!("gpu pair loop (f64, {} via {:?})", info.name, info.backend), work: Mutex::new(None), timing: Mutex::new((0.0, 0.0, 0)) })
    }

    fn make_work(&self, n: usize, ns: usize, ntypes: usize, cap: usize) -> Work {
        // Eight-byte units in each buffer, in the shader's binding order from 1.
        let sizes: [usize; 22] = [
            ns / 2 + 1,                  // types (u32)
            n * ns * 3,                  // rel
            n * ns,                      // qv
            n * ns,                      // sv
            ntypes * ntypes * 4 + 1,     // pair, then the Ewald width (0 for none)
            ntypes * ntypes,             // ovk
            ntypes * ntypes,             // bdamp
            ntypes,                      // sbase
            cap * 6,                     // lk
            cap,                         // lij (u32 pairs)
            n / 2 + cap + 2,             // inc (u32)
            cap * ns * 3,                // fa
            cap * ns * 3,                // fb
            cap * ns,                    // gqa
            cap * ns,                    // gsa
            cap * ns,                    // gqb
            cap * ns,                    // gsb
            cap * 5,                     // lscal
            n * ns * 3,                  // out_fsite
            n * 3,                       // out_fc
            n * ns,                      // out_gq
            n * ns,                      // out_gs
        ];
        let usage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST;
        let bufs: Vec<wgpu::Buffer> = sizes.iter().map(|&s| self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("pairs"), size: (s.max(1) * 8) as u64, usage, mapped_at_creation: false })).collect();
        let params = self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("params"), size: std::mem::size_of::<Params>() as u64, usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let mut entries = vec![wgpu::BindGroupEntry { binding: 0, resource: params.as_entire_binding() }];
        for (k, b) in bufs.iter().enumerate() {
            entries.push(wgpu::BindGroupEntry { binding: 1 + k as u32, resource: b.as_entire_binding() });
        }
        let bind = self.device.create_bind_group(&wgpu::BindGroupDescriptor { label: Some("pairs"), layout: &self.layout, entries: &entries });
        // The staging buffer: out_fsite, out_fc, out_gq, out_gs, then the per-link scalars (energy and virial are summed on the CPU).
        let out_sizes = [n * ns * 3, n * 3, n * ns, n * ns, cap * 5];
        let mut offsets = Vec::new();
        let mut total = 0usize;
        for s in out_sizes {
            offsets.push(total);
            total += s.max(1);
        }
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor { label: Some("staging"), size: (total * 8) as u64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        Work { n, ns, ntypes, cap, params, bufs, bind, staging, offsets }
    }

    fn run(&self, p: &PairProblem) -> Option<PairResult> {
        let t_all = std::time::Instant::now();
        let table = p.table;
        let ntypes = table.charge.len();
        let n = p.rel.len();
        let ns = p.types.len();
        let nl = p.links.len();
        if n == 0 || ns == 0 || nl == 0 {
            return None;
        }
        // A table with widths for some types and not others reads differently on the CPU.
        if !table.sigma.is_empty() && table.sigma.len() != ntypes {
            return None;
        }
        let mut guard = self.work.lock().expect("the pair work");
        let rebuild = match guard.as_ref() {
            Some(w) => w.n != n || w.ns != ns || w.ntypes != ntypes || w.cap < nl,
            None => true,
        };
        if rebuild {
            *guard = Some(self.make_work(n, ns, ntypes, nl + nl / 4 + 64));
        }
        let w = guard.as_ref().expect("the work");
        let q = &self.queue;
        let put = |i: usize, bytes: &[u8]| q.write_buffer(&w.bufs[i], 0, bytes);
        let ty: Vec<u32> = p.types.iter().map(|&t| t as u32).collect();
        put(TYPES, bytemuck::cast_slice(&ty));
        let rel: Vec<f64> = p.rel.iter().flat_map(|m| m.iter().flat_map(|v| [v.x, v.y, v.z])).collect();
        put(REL, bytemuck::cast_slice(&rel));
        let qv: Vec<f64> = p.charge.iter().flatten().cloned().collect();
        put(QV, bytemuck::cast_slice(&qv));
        let sv: Vec<f64> = p.sigma.iter().flatten().cloned().collect();
        put(SV, bytemuck::cast_slice(&sv));
        let mut pair: Vec<f64> = table.pair.iter().flat_map(|a| a.iter().cloned()).collect();
        pair.push(p.ewald_width.unwrap_or(0.0));
        put(PAIR, bytemuck::cast_slice(&pair));
        let mut ovk = vec![0.0f64; ntypes * ntypes];
        for (k, v) in table.overlap.iter().enumerate().take(ntypes * ntypes) {
            ovk[k] = *v;
        }
        put(OVK, bytemuck::cast_slice(&ovk));
        let mut bd = vec![0.0f64; ntypes * ntypes];
        for (k, v) in table.damp.iter().enumerate().take(ntypes * ntypes) {
            bd[k] = *v;
        }
        put(BDAMP, bytemuck::cast_slice(&bd));
        let sbase: Vec<f64> = (0..ntypes).map(|t| table.sigma.get(t).copied().unwrap_or(-1.0)).collect();
        put(SBASE, bytemuck::cast_slice(&sbase));
        let mut lk = Vec::with_capacity(nl * 6);
        let mut lij: Vec<u32> = Vec::with_capacity(nl * 2);
        for l in p.links {
            lk.extend_from_slice(&[l.d.x, l.d.y, l.d.z, l.weight, l.dweight, 0.0]);
            lij.extend_from_slice(&[l.i as u32, l.j as u32]);
        }
        put(LK, bytemuck::cast_slice(&lk));
        put(LIJ, bytemuck::cast_slice(&lij));
        let mut count = vec![0u32; n + 1];
        for l in p.links {
            count[l.i + 1] += 1;
            count[l.j + 1] += 1;
        }
        for m in 0..n {
            count[m + 1] += count[m];
        }
        let mut inc = vec![0u32; n + 1 + 2 * nl];
        inc[..n + 1].copy_from_slice(&count);
        let mut fill = count.clone();
        for (k, l) in p.links.iter().enumerate() {
            inc[n + 1 + fill[l.i] as usize] = (2 * k) as u32;
            fill[l.i] += 1;
            inc[n + 1 + fill[l.j] as usize] = (2 * k + 1) as u32;
            fill[l.j] += 1;
        }
        put(INC, bytemuck::cast_slice(&inc));
        q.write_buffer(&w.params, 0, bytemuck::bytes_of(&Params { n: n as u32, ns: ns as u32, nl: nl as u32, ntypes: ntypes as u32 }));

        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_bind_group(0, &w.bind, &[]);
            pass.set_pipeline(&self.pipelines[0]);
            pass.dispatch_workgroups(nl.div_ceil(64) as u32, 1, 1);
            pass.set_pipeline(&self.pipelines[1]);
            pass.dispatch_workgroups((n * ns) as u32, 1, 1);
        }
        let out: [(usize, usize); 5] = [(OUT_FSITE, n * ns * 3), (OUT_FC, n * 3), (OUT_GQ, n * ns), (OUT_GS, n * ns), (LSCAL, nl * 5)];
        for (k, (buf, len)) in out.iter().enumerate() {
            enc.copy_buffer_to_buffer(&w.bufs[*buf], 0, &w.staging, (w.offsets[k] * 8) as u64, (*len * 8) as u64);
        }
        let t = std::time::Instant::now();
        self.queue.submit(Some(enc.finish()));
        let slice = w.staging.slice(..);
        slice.map_async(wgpu::MapMode::Read, |r| r.expect("a result mapped"));
        self.device.poll(wgpu::PollType::wait_indefinitely()).expect("the GPU finished");
        let t_gpu = t.elapsed().as_secs_f64();
        let data: Vec<f64> = bytemuck::cast_slice(&slice.get_mapped_range().expect("the results range")).to_vec();
        w.staging.unmap();
        let v3 = |off: usize, i: usize| Vec3 { x: data[off + 3 * i], y: data[off + 3 * i + 1], z: data[off + 3 * i + 2] };
        let fsite: Vec<Vec<Vec3>> = (0..n).map(|m| (0..ns).map(|a| v3(w.offsets[0], m * ns + a)).collect()).collect();
        let fcentre: Vec<Vec3> = (0..n).map(|m| v3(w.offsets[1], m)).collect();
        let gq: Vec<Vec<f64>> = (0..n).map(|m| data[w.offsets[2] + m * ns..w.offsets[2] + (m + 1) * ns].to_vec()).collect();
        let gs: Vec<Vec<f64>> = (0..n).map(|m| data[w.offsets[3] + m * ns..w.offsets[3] + (m + 1) * ns].to_vec()).collect();
        let (mut energy, mut virial) = (0.0, 0.0);
        for k in 0..nl {
            energy += data[w.offsets[4] + 5 * k];
            virial += data[w.offsets[4] + 5 * k + 4];
        }
        let mut tm = self.timing.lock().expect("the timing");
        tm.0 += t_gpu;
        tm.1 += t_all.elapsed().as_secs_f64();
        tm.2 += 1;
        Some(PairResult { energy, virial, fsite, fcentre, gq, gs })
    }
}

impl PairEngine for GpuPairs {
    fn pairs(&self, problem: &PairProblem) -> Option<PairResult> {
        self.run(problem)
    }
}

/// Install the GPU's pair loop for this process, returning its name; or why
/// not, in which case the CPU's stays.
pub fn install_pairs() -> Result<String, String> {
    let gpu = GpuPairs::new()?;
    let name = gpu.name.clone();
    phys::liquid::set_pair_engine(Some(std::sync::Arc::new(gpu)));
    Ok(name)
}
