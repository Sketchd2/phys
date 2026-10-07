//! Dense products on the GPU in single precision, for Hartree-Fock and MP2 —
//! `docs/PLAY.md` Phase 6, E8a.
//!
//! The work that grows fastest in the correlated methods is one shape, `out = a
//! b^T`: the exchange matrix (`n^2 * nk` multiply-adds for each occupied
//! orbital), the transformation of the fitted integrals to virtual orbitals
//! (`virtuals * n * nk` for each), and MP2's pair integrals (`virtuals^2 * nk`
//! for each pair of occupied orbitals). For six waters the transformation alone
//! is some 4e14 multiply-adds a field, hours on the CPU and a few minutes here.
//! This implements the core crate's [`phys::electrons::hf::ProductEngine`] seam;
//! a program that wants it calls [`install_product`]. Without it nothing here
//! runs and every product is the CPU's, in double precision.
//!
//! # Precision
//!
//! The card does double at a thirty-second of single, so this is single
//! precision, accumulating in single across the whole inner dimension. That
//! gives relative errors of a few parts in a million at the longest inner
//! dimensions here; where the products enter as a sum of many signed terms (the
//! MP2 energy) it is harmless, and in the exchange matrix its effect on the
//! converged energy is second order. It is measured in `tests/product.rs` against the CPU, and
//! the owner's rule for single precision stands: nothing relies on it before its
//! precision is on record.
//!
//! # Long dispatches
//!
//! A display driver that runs a shader for more than about two seconds is reset
//! on Windows, so a product is done a slice of rows at a time, each slice sized
//! to keep a dispatch to a fraction of a second.

use phys::electrons::hf::ProductEngine;
use std::sync::Arc;
use wgpu::util::DeviceExt;

/// Multiply-adds one dispatch does, at most: a few tenths of a second on the
/// owner's RTX 2060, far inside Windows' two-second limit.
const WORK_PER_DISPATCH: f64 = 1.5e11;

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct Params {
    m: u32,
    n: u32,
    kd: u32,
    row0: u32,
    rows: u32,
    pad: [u32; 3],
}

/// The GPU and the compiled product shader.
pub struct GpuProduct {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    name: String,
}

impl GpuProduct {
    /// The first GPU wgpu finds, or why there is none.
    pub fn new() -> Result<GpuProduct, String> {
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
            label: Some("phys-gpu-product"),
            required_limits: adapter.limits(),
            ..Default::default()
        }))
        .map_err(|e| format!("no device: {e}"))?;
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("gemm"), source: wgpu::ShaderSource::Wgsl(include_str!("gemm.wgsl").into()) });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: Some("gemm"), layout: None, module: &module, entry_point: Some("main"), compilation_options: Default::default(), cache: None });
        Ok(GpuProduct { device, queue, pipeline, name: format!("gpu products (f32, {} via {:?})", info.name, info.backend) })
    }
}

/// Install the GPU's products for Hartree-Fock and MP2 in this process,
/// returning its name; or why not, in which case the CPU's stay.
pub fn install_product() -> Result<String, String> {
    let gpu = GpuProduct::new()?;
    let name = gpu.name.clone();
    phys::electrons::hf::set_product_engine(Some(Arc::new(gpu)));
    Ok(name)
}

impl ProductEngine for GpuProduct {
    fn product_nt(&self, a: &[f64], b: &[f64], m: usize, n: usize, kd: usize, out: &mut [f64]) {
        assert!(a.len() >= m * kd && b.len() >= n * kd && out.len() >= m * n);
        let dev = &self.device;
        let af: Vec<f32> = a[..m * kd].iter().map(|&x| x as f32).collect();
        let bf: Vec<f32> = b[..n * kd].iter().map(|&x| x as f32).collect();
        let a_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("a"), contents: bytemuck::cast_slice(&af), usage: wgpu::BufferUsages::STORAGE });
        let b_buf = dev.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("b"), contents: bytemuck::cast_slice(&bf), usage: wgpu::BufferUsages::STORAGE });
        drop((af, bf));
        // Rows of the result per dispatch: a multiple of the tile, enough
        // that the work is near the budget.
        let rows_per = (((WORK_PER_DISPATCH / (n as f64 * kd as f64)) as usize).max(64) / 64 * 64).min(m.div_ceil(64) * 64);
        let c_size = (rows_per.min(m) * n * 4) as u64;
        let c_buf = dev.create_buffer(&wgpu::BufferDescriptor { label: Some("c"), size: c_size, usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC, mapped_at_creation: false });
        let readback = dev.create_buffer(&wgpu::BufferDescriptor { label: Some("readback"), size: c_size, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let params_buf = dev.create_buffer(&wgpu::BufferDescriptor { label: Some("params"), size: std::mem::size_of::<Params>() as u64, usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let layout = self.pipeline.get_bind_group_layout(0);
        let bind = dev.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("gemm"),
            layout: &layout,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: a_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: b_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: c_buf.as_entire_binding() },
            ],
        });
        let mut row0 = 0usize;
        while row0 < m {
            let rows = rows_per.min(m - row0);
            let params = Params { m: m as u32, n: n as u32, kd: kd as u32, row0: row0 as u32, rows: rows as u32, pad: [0; 3] };
            self.queue.write_buffer(&params_buf, 0, bytemuck::bytes_of(&params));
            let mut encoder = dev.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("gemm") });
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor { label: Some("gemm"), timestamp_writes: None });
                pass.set_pipeline(&self.pipeline);
                pass.set_bind_group(0, &bind, &[]);
                pass.dispatch_workgroups(n.div_ceil(64) as u32, rows.div_ceil(64) as u32, 1);
            }
            encoder.copy_buffer_to_buffer(&c_buf, 0, &readback, 0, (rows * n * 4) as u64);
            self.queue.submit(Some(encoder.finish()));
            let slice = readback.slice(..(rows * n * 4) as u64);
            slice.map_async(wgpu::MapMode::Read, |r| r.expect("the product mapped"));
            dev.poll(wgpu::PollType::wait_indefinitely()).expect("the GPU finished its slice");
            {
                let data = slice.get_mapped_range().expect("the product's range");
                let values: &[f32] = bytemuck::cast_slice(&data);
                for (o, v) in out[row0 * n..(row0 + rows) * n].iter_mut().zip(values) {
                    *o = *v as f64;
                }
            }
            readback.unmap();
            row0 += rows;
        }
    }

    fn name(&self) -> &str {
        &self.name
    }
}
