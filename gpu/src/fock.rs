//! Hartree-Fock and MP2's heavy passes on the GPU, in double precision, against
//! whitened three-index integrals held in three tiers — `docs/PLAY.md` Phase 6,
//! E8a. This implements the core crate's [`phys::electrons::hf::FockEngine`]; a
//! program that wants it calls [`install_fock`].
//!
//! # What is kept, and why that is the point
//!
//! A Hartree-Fock iteration needs the integrals `(mn|P)` contracted with the
//! occupied orbitals, and the CPU reads the whole table to do it, whitens what
//! it read (a forward substitution for every entry, `nk^2/2` operations each),
//! and does both again next iteration. Here the whitening is done once, when a
//! table is first used, and `B[mn][k]` (for `m >= n`, `nk` doubles a pair) is
//! kept from then on. An iteration is then
//!
//! * the contraction `x[i][m][k] = sum_n B[mn][k] c[n][i]` (one pass over `B`
//!   for each block of orbitals), which also yields the density fitted onto the
//!   directions, `dt[k] = sum_p D[p] B[p][k]`, on the same visit to each segment,
//! * the exchange matrix `K = sum_ik x x^T` (a product, lower triangle of
//!   tiles only),
//! * and the Coulomb matrix `J[p] = sum_k B[p][k] dt[k]` (one more pass),
//!
//! with energy `1/2 dt . dt`. The same table serves every field of a
//! counterpoise correction, a cluster, a finite-field polarisability.
//!
//! # Three tiers, and what moves between them
//!
//! `B` is cut into segments of a fixed 256 MB (so that nothing about the
//! arithmetic depends on the machine), in pair order, and each segment lives in
//! exactly one tier, the lowest-numbered in the fastest:
//!
//! * **the card's memory** — the first segments, uploaded once and read in
//!   place at hundreds of gigabytes a second;
//! * **system memory** — the next, uploaded each pass through a two-buffer ring
//!   on the card;
//! * **a file** — the rest (`PHYS_SPILL_DIR`, as the raw table's spill), read
//!   by a thread that runs ahead of the card.
//!
//! A pass visits the segments in index order always, so every sum is added in
//! the same order whatever the tiers are and **the answer does not depend on how
//! much memory the machine had**: the same bits from a card that holds the lot
//! and from one that holds a tenth (`gpu/tests/fock.rs` asserts equality, not
//! closeness). The streamed segments follow the resident ones, so the reader
//! thread is already loading the first of them while the card works on the
//! resident segments, and from then on the next segment is read while the
//! current one is used.
//!
//! **What decides the split.** Everything that is not `B` has to fit as well:
//! the exchange work space (`x`, a block of orbitals at a time — more orbitals a
//! block means fewer passes over what is streamed), the exchange matrix, the
//! Coulomb vectors, and the ring if anything streams. Those are taken first
//! and what is left of the card holds `B`; the card's free memory is asked of
//! `nvidia-smi` where there is one (`PHYS_GPU_GB` overrides it), because wgpu
//! does not report memory and a Windows driver will quietly page an allocation
//! that does not fit into system memory at a tenth of the speed. System memory
//! is 60% of what is free when the table is loaded, after the raw table has
//! taken its share; the raw table is let go of once this has its own copy
//! ([`phys::electrons::scf::Fitted::release_raw`]). If the file is not enough,
//! or there is none, `load` says no and the CPU does it as it always did.
//!
//! **Why the passes are not the cost to be cut.** A hexamer's `B` is 37 GB; a
//! pass over it from system memory is three seconds at PCIe's 12 GB/s. The
//! exchange product is `~10^13` double-precision operations an iteration, a
//! minute on this card and a good fraction of that on a CPU. So the block of
//! orbitals is a trade of passes for memory and nothing cleverer is attempted:
//! the arithmetic is the bound, and at the sizes where the card holds the lot
//! (a dimer) the tiers do nothing at all.
//!
//! # Precision, limits
//!
//! All double: the answers are the CPU's to round-off. A dispatch is kept to a
//! few billion multiply-adds so that Windows does not reset the display driver.

use phys::electrons::hf::FockEngine;
use phys::electrons::linalg::Matrix;
use phys::electrons::scf::{free_physical_memory, SpillTo};
use std::sync::Mutex;
use wgpu::util::DeviceExt;

/// Multiply-adds one dispatch does, at most.
const WORK_PER_DISPATCH: f64 = 3e9;
/// The bytes a segment of `B` holds, unless a test asks for a different cut.
const SEGMENT_BYTES: usize = 256 << 20;
/// Orbitals the contraction shader holds at a time, at most.
const BLOCK: usize = 32;
/// The exchange work space `x` is given at most this much of the card.
const X_BYTES: usize = 1_200 << 20;
/// Left unallocated on the card for the display and the driver.
const MARGIN: usize = 600 << 20;

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

/// What the tiers may use. [`Budget::from_environment`] is what a program gets;
/// a test builds one to force a split on a table that would fit whole.
#[derive(Debug, Clone)]
pub struct Budget {
    /// Bytes of the card's memory for everything this engine allocates.
    pub vram: usize,
    /// Bytes of system memory for segments, or `None` for 60% of what is free
    /// when the table is loaded.
    pub ram: Option<usize>,
    /// Where segments that fit in neither go, and how much of the disk.
    pub disk: Option<SpillTo>,
    /// The bytes of a segment of `B`.
    pub segment_bytes: usize,
    /// The most bytes of the exchange work space.
    pub x_bytes: usize,
}

impl Budget {
    /// The card's free memory less a margin (`PHYS_GPU_GB` instead where set;
    /// 2 GB where nothing can say), `PHYS_STORE_RAM_GB` for system memory
    /// where set, and `PHYS_SPILL_DIR` for the file.
    pub fn from_environment() -> Budget {
        let gb = |name: &str| std::env::var(name).ok().and_then(|v| v.trim().parse::<f64>().ok()).map(|g| (g * 1e9) as usize);
        let vram = gb("PHYS_GPU_GB").unwrap_or_else(|| free_card_memory().map(|b| b.saturating_sub(MARGIN)).unwrap_or(2_000_000_000));
        Budget { vram, ram: gb("PHYS_STORE_RAM_GB"), disk: SpillTo::from_environment(), segment_bytes: SEGMENT_BYTES, x_bytes: X_BYTES }
    }
}

/// The first NVIDIA card's free memory in bytes, from `nvidia-smi`.
fn free_card_memory() -> Option<usize> {
    let out = std::process::Command::new("nvidia-smi").args(["--query-gpu=memory.free", "--format=csv,noheader,nounits"]).output().ok()?;
    let mib: usize = String::from_utf8_lossy(&out.stdout).lines().next()?.trim().parse().ok()?;
    Some(mib << 20)
}

/// How a table of `n` functions and `nk` directions is split, from the budget.
#[derive(Debug, Clone, PartialEq)]
pub struct Plan {
    pub seg_pairs: usize,
    pub segments: usize,
    pub on_card: usize,
    pub in_memory: usize,
    pub on_disk: usize,
    /// Orbitals of `x` a block holds.
    pub block: usize,
}

/// The split of a table, or `None` if it does not fit anywhere.
pub fn plan(n: usize, nk: usize, budget: &Budget, free_ram: Option<u64>) -> Option<Plan> {
    let npairs = n * (n + 1) / 2;
    let seg_pairs = (budget.segment_bytes / (nk * 8)).max(1);
    let segments = npairs.div_ceil(seg_pairs);
    let seg_bytes = seg_pairs * nk * 8;
    let total = npairs * nk * 8;
    // What is not B: the exchange matrix, the density and the Coulomb vector
    // by pair, the fitted density.
    let fixed = n * n * 8 + 2 * npairs * 8 + nk * 8;
    let x_one = n * nk * 8;
    let block = (budget.x_bytes / x_one).clamp(1, BLOCK);
    let (on_card, block) = if total + fixed + block * x_one <= budget.vram {
        (segments, block)
    } else {
        let ring = 2 * seg_bytes;
        let mut b = block;
        loop {
            let used = fixed + b * x_one + ring;
            if used <= budget.vram {
                break ((budget.vram - used) / seg_bytes, b);
            }
            if b == 1 {
                return None;
            }
            b -= 1;
        }
    };
    let on_card = on_card.min(segments);
    let rest = segments - on_card;
    let ram = budget.ram.unwrap_or_else(|| (free_ram.unwrap_or(0) as f64 * 0.6) as usize);
    let in_memory = rest.min(ram / seg_bytes);
    let on_disk = rest - in_memory;
    if on_disk > 0 {
        let to = budget.disk.as_ref()?;
        if (on_disk * seg_bytes) as u64 > to.cap_bytes {
            return None;
        }
    }
    Some(Plan { seg_pairs, segments, on_card, in_memory, on_disk, block })
}

fn pidx(a: usize, b: usize) -> usize {
    if a >= b {
        a * (a + 1) / 2 + b
    } else {
        b * (b + 1) / 2 + a
    }
}

fn read_at(file: &std::fs::File, buf: &mut [u8], offset: u64) -> std::io::Result<()> {
    #[cfg(unix)]
    {
        std::os::unix::fs::FileExt::read_exact_at(file, buf, offset)
    }
    #[cfg(windows)]
    {
        let mut done = 0;
        while done < buf.len() {
            let n = std::os::windows::fs::FileExt::seek_read(file, &mut buf[done..], offset + done as u64)?;
            if n == 0 {
                return Err(std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "the store is short"));
            }
            done += n;
        }
        Ok(())
    }
}

fn write_at(file: &std::fs::File, buf: &[u8], offset: u64) -> std::io::Result<()> {
    #[cfg(unix)]
    {
        std::os::unix::fs::FileExt::write_all_at(file, buf, offset)
    }
    #[cfg(windows)]
    {
        let mut done = 0;
        while done < buf.len() {
            done += std::os::windows::fs::FileExt::seek_write(file, &buf[done..], offset + done as u64)?;
        }
        Ok(())
    }
}

/// The file the last tier lives in, deleted when it is dropped.
struct DiskStore {
    path: std::path::PathBuf,
    file: std::fs::File,
}

impl Drop for DiskStore {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

enum Tier {
    Card(wgpu::Buffer),
    Memory(Vec<f64>),
    Disk,
}

/// What a pass hands over for a streamed segment.
enum Chunk<'a> {
    Borrowed(&'a [f64]),
    Owned(Vec<f64>),
}

struct Resident {
    key: usize,
    n: usize,
    nk: usize,
    npairs: usize,
    plan: Plan,
    tiers: Vec<Tier>,
    disk: Option<DiskStore>,
    /// Two buffers the streamed segments are uploaded into by turns.
    ring: Vec<wgpu::Buffer>,
    /// The first segment that is not on the card.
    first_streamed: usize,
}

impl Resident {
    fn pairs_of(&self, s: usize) -> (usize, usize) {
        (s * self.plan.seg_pairs, ((s + 1) * self.plan.seg_pairs).min(self.npairs))
    }

    fn disk_offset(&self, s: usize) -> u64 {
        ((s - self.plan.on_card - self.plan.in_memory) * self.plan.seg_pairs * self.nk * 8) as u64
    }
}

/// The GPU, its shaders, and the integrals currently kept.
pub struct GpuFock {
    device: wgpu::Device,
    queue: wgpu::Queue,
    contract: wgpu::ComputePipeline,
    kernel_k: wgpu::ComputePipeline,
    dt: wgpu::ComputePipeline,
    jp: wgpu::ComputePipeline,
    name: String,
    budget: Budget,
    resident: Mutex<Option<Resident>>,
}

impl GpuFock {
    /// The first GPU wgpu finds that offers double precision in shaders, with
    /// the budget the environment gives, or why not.
    pub fn new() -> Result<GpuFock, String> {
        GpuFock::with_budget(Budget::from_environment())
    }

    /// As [`GpuFock::new`] with a budget of the caller's.
    pub fn with_budget(budget: Budget) -> Result<GpuFock, String> {
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
        Ok(GpuFock { device, queue, contract, kernel_k, dt, jp, name: format!("gpu Fock (f64, {} via {:?})", info.name, info.backend), budget, resident: Mutex::new(None) })
    }

    /// The adapter's name, for a message.
    pub fn name_string(&self) -> String {
        self.name.clone()
    }

    /// How the table now kept is split, if one is.
    pub fn current_plan(&self) -> Option<Plan> {
        self.resident.lock().expect("the resident integrals").as_ref().map(|r| r.plan.clone())
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

    /// Visit every segment of `B` in index order, handing `f` the segment's
    /// number and a buffer holding it: the card's own for the first, the ring
    /// for the rest. The streamed segments are read ahead by a thread that
    /// starts before the first is needed, so it is already at work while the
    /// resident ones are used.
    fn pass(&self, r: &Resident, f: &mut dyn FnMut(usize, &wgpu::Buffer)) {
        let segments = r.plan.segments;
        let first = r.first_streamed;
        std::thread::scope(|sc| {
            let (tx, rx) = std::sync::mpsc::sync_channel::<Chunk>(2);
            let (give_back, recycled) = std::sync::mpsc::channel::<Vec<f64>>();
            if first < segments {
                sc.spawn(move || {
                    for s in first..segments {
                        let chunk = match &r.tiers[s] {
                            Tier::Memory(v) => Chunk::Borrowed(v),
                            Tier::Disk => {
                                let (p0, p1) = r.pairs_of(s);
                                let len = (p1 - p0) * r.nk;
                                let mut buf = recycled.try_recv().unwrap_or_default();
                                buf.resize(len, 0.0);
                                let store = r.disk.as_ref().expect("the store's file");
                                read_at(&store.file, bytemuck::cast_slice_mut(&mut buf[..]), r.disk_offset(s)).unwrap_or_else(|e| panic!("the store {} cannot be read: {e}", store.path.display()));
                                Chunk::Owned(buf)
                            }
                            Tier::Card(_) => unreachable!("a segment on the card is not streamed"),
                        };
                        if tx.send(chunk).is_err() {
                            return;
                        }
                    }
                });
            }
            for s in 0..first {
                if let Tier::Card(buf) = &r.tiers[s] {
                    f(s, buf);
                }
            }
            for s in first..segments {
                let chunk = rx.recv().expect("the reader thread ended early");
                let ring = &r.ring[(s - first) % 2];
                match &chunk {
                    Chunk::Borrowed(v) => self.queue.write_buffer(ring, 0, bytemuck::cast_slice(v)),
                    Chunk::Owned(v) => self.queue.write_buffer(ring, 0, bytemuck::cast_slice(v)),
                }
                f(s, ring);
                if let Chunk::Owned(v) = chunk {
                    let _ = give_back.send(v);
                }
            }
        });
    }

    /// `x[i][m][k] += sum_n B(m, n)[k] c[n][i]` over the pairs of segment `s`.
    fn contract_segment(&self, r: &Resident, s: usize, seg: &wgpu::Buffer, cb: &wgpu::Buffer, x: &wgpu::Buffer, cnt: usize) {
        let (n, nk) = (r.n, r.nk);
        let (p0, p1) = r.pairs_of(s);
        let rows_per = ((WORK_PER_DISPATCH / (n as f64 * nk as f64 * cnt as f64)) as usize).clamp(1, n);
        let mut m0 = 0;
        while m0 < n {
            let rows = rows_per.min(n - m0);
            let params = Params { n: n as u32, nk: nk as u32, count: cnt as u32, p0: p0 as u32, p1: p1 as u32, m0: m0 as u32, rows: rows as u32, pad: 0 };
            self.run(&self.contract, &params, &[seg, cb, x], (nk.div_ceil(64) as u32, rows as u32));
            m0 += rows;
        }
    }

    /// The columns `first .. first + cnt` of `c`, packed `n x cnt`.
    fn block_of(&self, c: &[f64], stride: usize, n: usize, first: usize, cnt: usize) -> wgpu::Buffer {
        let mut cb = vec![0.0f64; n * cnt];
        for a in 0..n {
            for i in 0..cnt {
                cb[a * cnt + i] = c[a * stride + first + i];
            }
        }
        self.upload("c", &cb)
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
        // The last table goes first: its memory is part of what this one is
        // sized against.
        *self.resident.lock().expect("the resident integrals") = None;
        let Some(plan) = plan(n, nk, &self.budget, free_physical_memory()) else {
            eprintln!("  whitened integrals: {:.1} GB do not fit the card, memory and disk allowed; the CPU does this one", (n * (n + 1) / 2 * nk * 8) as f64 / 1e9);
            return false;
        };
        let npairs = n * (n + 1) / 2;
        let seg_bytes = plan.seg_pairs * nk * 8;
        let first_streamed = plan.on_card;
        let scope = self.device.push_error_scope(wgpu::ErrorFilter::OutOfMemory);
        let mut tiers = Vec::with_capacity(plan.segments);
        for s in 0..plan.segments {
            let pairs = plan.seg_pairs.min(npairs - s * plan.seg_pairs);
            tiers.push(if s < plan.on_card {
                Tier::Card(self.storage("B", pairs * nk * 8))
            } else if s < plan.on_card + plan.in_memory {
                Tier::Memory(vec![0.0f64; pairs * nk])
            } else {
                Tier::Disk
            });
        }
        let ring: Vec<wgpu::Buffer> = if first_streamed < plan.segments { (0..2).map(|_| self.storage("ring", seg_bytes)).collect() } else { Vec::new() };
        if pollster::block_on(scope.pop()).is_some() {
            eprintln!("  whitened integrals: the card refused its share; the CPU does this one");
            return false;
        }
        let disk = if plan.on_disk > 0 {
            let to = self.budget.disk.as_ref().expect("planned a disk tier with a disk");
            static COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
            let made = std::fs::create_dir_all(&to.dir).and_then(|_| {
                let path = to.dir.join(format!("phys-store-{}-{}.bin", std::process::id(), COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)));
                let file = std::fs::OpenOptions::new().read(true).write(true).create(true).truncate(true).open(&path)?;
                let store = DiskStore { path, file };
                store.file.set_len((plan.on_disk * seg_bytes) as u64)?;
                Ok(store)
            });
            match made {
                Ok(d) => Some(d),
                Err(e) => {
                    eprintln!("  whitened integrals: the store file could not be made in {}: {e}; the CPU does this one", to.dir.display());
                    return false;
                }
            }
        } else {
            None
        };
        let mut r = Resident { key, n, nk, npairs, plan: plan.clone(), tiers, disk, ring, first_streamed };
        // Entries arrive in the table's order, which is not pair order, so a
        // file's writes are batched and sorted by offset before they are made.
        let mut pending: Vec<(u64, Vec<f64>)> = Vec::new();
        let mut pending_bytes = 0usize;
        let mut written = 0usize;
        let mut failed: Option<String> = None;
        let flush = |pending: &mut Vec<(u64, Vec<f64>)>, store: &Option<DiskStore>, failed: &mut Option<String>| {
            pending.sort_by_key(|e| e.0);
            if let Some(store) = store {
                for (off, v) in pending.iter() {
                    if let Err(e) = write_at(&store.file, bytemuck::cast_slice(v), *off) {
                        *failed = Some(format!("{}: {e}", store.path.display()));
                    }
                }
            }
            pending.clear();
        };
        {
            let (tiers, disk) = (&mut r.tiers, &r.disk);
            let disk_first = plan.on_card + plan.in_memory;
            source(&mut |m, nn, v| {
                let p = pidx(m, nn);
                let (s, at) = (p / plan.seg_pairs, (p % plan.seg_pairs) * nk);
                match &mut tiers[s] {
                    Tier::Card(buf) => {
                        self.queue.write_buffer(buf, (at * 8) as u64, bytemuck::cast_slice(v));
                        if written % 256 == 0 {
                            self.queue.submit(std::iter::empty());
                        }
                    }
                    Tier::Memory(vec) => vec[at..at + nk].copy_from_slice(v),
                    Tier::Disk => {
                        pending.push((((s - disk_first) * plan.seg_pairs * nk + at) as u64 * 8, v.to_vec()));
                        pending_bytes += nk * 8;
                        if pending_bytes >= (128 << 20) {
                            flush(&mut pending, disk, &mut failed);
                            pending_bytes = 0;
                        }
                    }
                }
                written += 1;
            });
            flush(&mut pending, disk, &mut failed);
        }
        if let Some(e) = failed {
            eprintln!("  whitened integrals: writing the store failed ({e}); the CPU does this one");
            return false;
        }
        self.queue.submit(std::iter::empty());
        self.device.poll(wgpu::PollType::wait_indefinitely()).expect("the integrals arrived");
        // Pairs of functions whose product is negligible are not in the table;
        // their integrals are zero, which is what new storage holds.
        assert!(written <= npairs, "the table gave {written} entries for {npairs} pairs");
        eprintln!(
            "  whitened integrals: {:.1} GB in {} segments of {} MB: {} on the card, {} in memory, {} on disk; orbitals {} to a block",
            (npairs * nk * 8) as f64 / 1e9,
            plan.segments,
            seg_bytes >> 20,
            plan.on_card,
            plan.in_memory,
            plan.on_disk,
            plan.block
        );
        *self.resident.lock().expect("the resident integrals") = Some(r);
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
        drop(dp);
        let dtb = self.storage("dt", nk * 8);
        let jpb = self.storage("jp", npairs * 8);
        let kmat = self.storage("kmat", n * n * 8);
        // Exchange, the orbitals a block at a time. The first block's visit to
        // each segment also fits the density (the Coulomb matrix's right-hand
        // side), so that costs no pass of its own.
        let tiles = n.div_ceil(64);
        for (b, first) in (0..count).step_by(r.plan.block).enumerate() {
            let cnt = r.plan.block.min(count - first);
            let cb = self.block_of(c, stride, n, first, cnt);
            let x = self.storage("x", cnt * n * nk * 8);
            self.pass(r, &mut |s, seg| {
                self.contract_segment(r, s, seg, &cb, &x, cnt);
                if b == 0 {
                    let (p0, p1) = r.pairs_of(s);
                    let params = Params { n: n as u32, nk: nk as u32, count: 0, p0: p0 as u32, p1: p1 as u32, m0: 0, rows: 0, pad: 0 };
                    self.run(&self.dt, &params, &[seg, &dpb, &dtb], (nk.div_ceil(64) as u32, 1));
                }
            });
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
        // Coulomb, once the density is fitted.
        self.pass(r, &mut |s, seg| {
            let (p0, p1) = r.pairs_of(s);
            let pairs = p1 - p0;
            let gx = pairs.min(65535);
            let params = Params { n: n as u32, nk: nk as u32, count: 0, p0: p0 as u32, p1: p1 as u32, m0: 0, rows: 0, pad: 0 };
            self.run(&self.jp, &params, &[seg, &dtb, &jpb], (gx as u32, pairs.div_ceil(gx) as u32));
        });
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
        for f in (0..count).step_by(r.plan.block) {
            let cnt = r.plan.block.min(count - f);
            let cb = self.block_of(c, stride, n, first + f, cnt);
            let x = self.storage("x", cnt * n * nk * 8);
            self.pass(r, &mut |s, seg| self.contract_segment(r, s, seg, &cb, &x, cnt));
            out.extend(self.read(&x, cnt * n * nk));
        }
        out
    }

    fn name(&self) -> &str {
        &self.name
    }
}
