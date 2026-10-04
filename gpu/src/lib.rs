//! The non-local correlation's double sum on the GPU, in single precision —
//! `docs/PLAY.md` Phase 6, E7.
//!
//! # Why this exists, and what it is not
//!
//! The double sum is the one part of the van der Waals density functional
//! that grows as the square of the grid's points, and on the owner's machine
//! it was most of a dimer's cost: 20 s an evaluation at (50, 12), minutes at
//! (75, 18). It is shaped for a GPU — billions of independent pairs over a
//! small table — but this card does double precision at 1/32 of single, no
//! faster than the twelve-thread CPU, and wgpu's shading language has no
//! double at all. So the owner had it built in single precision, each row
//! summed with Kahan compensation and the rows added in double on the CPU,
//! with its precision to be measured after E8 before anything relies on it.
//! Until then the CPU's [`phys::electrons::vdw::CpuRows`] stays the reference
//! and the default.
//!
//! The core crate carries no GPU code (it has no dependencies and builds for
//! wasm32); this crate implements its [`RowEngine`] seam, and a program that
//! wants it calls [`install`].
//!
//! # Windows and long dispatches
//!
//! Windows resets a display driver whose GPU work runs past about two
//! seconds. The columns are therefore added in slices, one dispatch each,
//! sized to keep a dispatch short, with each row's running sums and their
//! compensation carried from one to the next.

use phys::electrons::vdw::{KernelTable, RowEngine, RowPoint, RowSums, ASYMPTOTE_C, ASYMPTOTIC_FROM};
use std::sync::{Arc, Mutex};
use wgpu::util::DeviceExt;

/// Pairs one dispatch adds, at most: about a tenth of a second on the
/// owner's RTX 2060, far inside Windows' two-second limit.
const PAIRS_PER_DISPATCH: u64 = 1 << 31;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    n: u32,
    j0: u32,
    j1: u32,
    m: u32,
    d_max: f32,
    scale: f32,
    asym_from: f32,
    asym_c: f32,
    with_c: u32,
    _pad: [u32; 3],
}

/// The GPU and the compiled shader.
pub struct GpuRows {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    name: String,
    /// The kernel table last uploaded, by its address and size, so that a
    /// field's every iteration does not upload it again.
    table: Mutex<Option<(usize, usize, Arc<wgpu::Buffer>)>>,
}

impl GpuRows {
    /// The first GPU wgpu finds, or why there is none.
    pub fn new() -> Result<GpuRows, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            force_fallback_adapter: false,
            compatible_surface: None,
            ..Default::default()
        }))
        .map_err(|e| format!("no adapter: {e}"))?;
        let info = adapter.get_info();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("phys-gpu"),
            required_limits: adapter.limits(),
            ..Default::default()
        }))
        .map_err(|e| format!("no device: {e}"))?;
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("rows"),
            source: wgpu::ShaderSource::Wgsl(include_str!("rows.wgsl").into()),
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("rows"),
            layout: None,
            module: &module,
            entry_point: Some("main"),
            compilation_options: Default::default(),
            cache: None,
        });
        Ok(GpuRows { device, queue, pipeline, name: format!("gpu (f32, {} via {:?})", info.name, info.backend), table: Mutex::new(None) })
    }

    fn table_buffer(&self, table: &KernelTable) -> Arc<wgpu::Buffer> {
        let key = (table.values.as_ptr() as usize, table.values.len());
        let mut cached = self.table.lock().expect("the table cache");
        if let Some((a, l, b)) = cached.as_ref() {
            if (*a, *l) == key {
                return b.clone();
            }
        }
        let values: Vec<f32> = table.values.iter().map(|v| *v as f32).collect();
        let buffer = Arc::new(self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("kernel table"),
            contents: bytemuck::cast_slice(&values),
            usage: wgpu::BufferUsages::STORAGE,
        }));
        *cached = Some((key.0, key.1, buffer.clone()));
        buffer
    }
}

/// Install the GPU's rows for this process, returning its name; or why not,
/// in which case the CPU's stay.
pub fn install() -> Result<String, String> {
    let gpu = GpuRows::new()?;
    let name = gpu.name.clone();
    phys::electrons::vdw::set_row_engine(Some(Arc::new(gpu)));
    Ok(name)
}

impl RowEngine for GpuRows {
    fn rows(&self, points: &[RowPoint], table: &KernelTable, with_c: bool) -> Vec<RowSums> {
        let n = points.len();
        if n == 0 {
            return Vec::new();
        }
        // Positions relative to the points' centre, so that single precision
        // spends its digits on separations rather than on where the molecule
        // happens to sit.
        let mut centre = [0.0f64; 3];
        for p in points {
            for k in 0..3 {
                centre[k] += p.r[k] / n as f64;
            }
        }
        let pts: Vec<[f32; 4]> = points.iter().map(|p| [(p.r[0] - centre[0]) as f32, (p.r[1] - centre[1]) as f32, (p.r[2] - centre[2]) as f32, p.wn as f32]).collect();
        let qs: Vec<f32> = points.iter().map(|p| p.q as f32).collect();
        let dev = &self.device;
        let pts_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("points"), contents: bytemuck::cast_slice(&pts), usage: wgpu::BufferUsages::STORAGE });
        let qs_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("q"), contents: bytemuck::cast_slice(&qs), usage: wgpu::BufferUsages::STORAGE });
        let acc_size = (n * 10 * 4) as u64;
        let acc = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("sums"), contents: bytemuck::cast_slice(&vec![0.0f32; n * 10]), usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC });
        let readback = dev.create_buffer(&wgpu::BufferDescriptor { label: Some("readback"), size: acc_size, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let params_buf = dev.create_buffer(&wgpu::BufferDescriptor { label: Some("params"), size: std::mem::size_of::<Params>() as u64, usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let table_buf = self.table_buffer(table);
        let layout = self.pipeline.get_bind_group_layout(0);
        let bind = dev.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("rows"),
            layout: &layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: pts_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: qs_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: table_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: acc.as_entire_binding() },
            ],
        });
        let u_max = table.d_max / (1.0 + table.d_max);
        let columns = ((PAIRS_PER_DISPATCH / n as u64).max(256) as usize).min(n);
        let groups = n.div_ceil(256) as u32;
        let mut j0 = 0usize;
        while j0 < n {
            let j1 = (j0 + columns).min(n);
            let params = Params {
                n: n as u32,
                j0: j0 as u32,
                j1: j1 as u32,
                m: table.m as u32,
                d_max: table.d_max as f32,
                scale: ((table.m - 1) as f64 / u_max) as f32,
                asym_from: ASYMPTOTIC_FROM as f32,
                asym_c: ASYMPTOTE_C as f32,
                with_c: with_c as u32,
                _pad: [0; 3],
            };
            self.queue.write_buffer(&params_buf, 0, bytemuck::bytes_of(&params));
            let mut encoder = dev.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("rows") });
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor { label: Some("rows"), timestamp_writes: None });
                pass.set_pipeline(&self.pipeline);
                pass.set_bind_group(0, &bind, &[]);
                pass.dispatch_workgroups(groups, 1, 1);
            }
            self.queue.submit(Some(encoder.finish()));
            // One dispatch at a time, so that no submission holds the GPU
            // for longer than one slice.
            dev.poll(wgpu::PollType::wait_indefinitely()).expect("the GPU finished its slice");
            j0 = j1;
        }
        let mut encoder = dev.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("readback") });
        encoder.copy_buffer_to_buffer(&acc, 0, &readback, 0, acc_size);
        self.queue.submit(Some(encoder.finish()));
        let slice = readback.slice(..);
        slice.map_async(wgpu::MapMode::Read, |r| r.expect("the sums mapped"));
        dev.poll(wgpu::PollType::wait_indefinitely()).expect("the sums arrived");
        let data: Vec<f32> = bytemuck::cast_slice(&slice.get_mapped_range().expect("the sums' range")).to_vec();
        readback.unmap();
        (0..n).map(|i| {
            let s = &data[i * 10..i * 10 + 5];
            RowSums { a: s[0] as f64, b: s[1] as f64, c: if with_c { [s[2] as f64, s[3] as f64, s[4] as f64] } else { [0.0; 3] } }
        }).collect()
    }

    fn name(&self) -> &str {
        &self.name
    }
}
