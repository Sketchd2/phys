fn main() {
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        force_fallback_adapter: false,
        ..Default::default()
    }));
    match adapter {
        Ok(a) => println!("{:?}", a.get_info()),
        Err(e) => println!("none: {e}"),
    }
}
