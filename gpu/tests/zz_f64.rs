//! Scratch: does double precision work in a shader on this adapter, and how fast is a product in it
//! against the CPU's. Prints, asserts nothing.

use phys::electrons::hf::product_cpu;
use std::time::Instant;
use wgpu::util::DeviceExt;

#[test]
fn double_precision_product_against_the_cpu() {
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions { power_preference: wgpu::PowerPreference::HighPerformance, force_fallback_adapter: false, compatible_surface: None, ..Default::default() })).expect("an adapter");
    println!("ADAPTER {} f64 offered {}", adapter.get_info().name, adapter.features().contains(wgpu::Features::SHADER_F64));
    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor { label: Some("f64"), required_features: wgpu::Features::SHADER_F64, required_limits: adapter.limits(), ..Default::default() })).expect("a device with f64");
    let src = include_str!("../src/gemm.wgsl").replace("f32", "f64").replace("0.0;", "0.0lf;");
    let module = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("gemm64"), source: wgpu::ShaderSource::Wgsl(src.into()) });
    let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: Some("gemm64"), layout: None, module: &module, entry_point: Some("main"), compilation_options: Default::default(), cache: None });
    println!("the f64 shader compiled");
    for (m, n, kd) in [(1254usize, 513usize, 4000usize), (428, 428, 1772)] {
        let mut seed = 99u64;
        let mut fill = |len: usize| -> Vec<f64> {
            (0..len).map(|_| { seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); ((seed >> 11) as f64 / (1u64 << 53) as f64) - 0.5 }).collect()
        };
        let a = fill(m * kd);
        let b = fill(n * kd);
        let t = Instant::now();
        let mut cpu = vec![0.0; m * n];
        product_cpu(&a, &b, m, n, kd, &mut cpu);
        let t_cpu = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let a_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("a"), contents: bytemuck::cast_slice(&a), usage: wgpu::BufferUsages::STORAGE });
        let b_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("b"), contents: bytemuck::cast_slice(&b), usage: wgpu::BufferUsages::STORAGE });
        let c_buf = device.create_buffer(&wgpu::BufferDescriptor { label: Some("c"), size: (m * n * 8) as u64, usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC, mapped_at_creation: false });
        let readback = device.create_buffer(&wgpu::BufferDescriptor { label: Some("rb"), size: (m * n * 8) as u64, usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST, mapped_at_creation: false });
        let params: [u32; 8] = [m as u32, n as u32, kd as u32, 0, m as u32, 0, 0, 0];
        let p_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("p"), contents: bytemuck::cast_slice(&params), usage: wgpu::BufferUsages::UNIFORM });
        let layout = pipeline.get_bind_group_layout(0);
        let bind = device.create_bind_group(&wgpu::BindGroupDescriptor { label: Some("g"), layout: &layout, entries: &[
            wgpu::BindGroupEntry { binding: 0, resource: p_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 1, resource: a_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 2, resource: b_buf.as_entire_binding() },
            wgpu::BindGroupEntry { binding: 3, resource: c_buf.as_entire_binding() },
        ] });
        let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("e") });
        {
            let mut pass = enc.begin_compute_pass(&wgpu::ComputePassDescriptor { label: Some("p"), timestamp_writes: None });
            pass.set_pipeline(&pipeline);
            pass.set_bind_group(0, &bind, &[]);
            pass.dispatch_workgroups(n.div_ceil(64) as u32, m.div_ceil(64) as u32, 1);
        }
        enc.copy_buffer_to_buffer(&c_buf, 0, &readback, 0, (m * n * 8) as u64);
        queue.submit(Some(enc.finish()));
        let slice = readback.slice(..);
        slice.map_async(wgpu::MapMode::Read, |r| r.expect("mapped"));
        device.poll(wgpu::PollType::wait_indefinitely()).expect("done");
        let got: Vec<f64> = bytemuck::cast_slice(&slice.get_mapped_range().expect("range")).to_vec();
        let t_gpu = t.elapsed().as_secs_f64();
        let worst = cpu.iter().zip(&got).map(|(x, y)| (x - y).abs()).fold(0.0f64, f64::max);
        let flops = 2.0 * m as f64 * n as f64 * kd as f64;
        println!("PRODUCT {m} x {n} x {kd}: CPU (all threads, double) {t_cpu:.3} s = {:.1} GFLOP/s; GPU double {t_gpu:.3} s = {:.1} GFLOP/s incl. transfers; worst difference {worst:.2e}", flops / t_cpu / 1e9, flops / t_gpu / 1e9);
    }
}
