#!/usr/bin/env bash
cd /e/Projects/phys/gpu && CARGO_TARGET_DIR=../target/dev5 cargo build --release --bin phys-flexbulk-gpu > /tmp/bg.out 2>&1
cd /e/Projects/phys/runs/e8c/timing
for b in phys-flexbulk phys-flexbulk-gpu; do echo "== $b" >> /tmp/bg.out; PHYS_GPU_PROFILE=1 /e/Projects/phys/target/dev5/release/$b.exe law-final.txt 298 216 1.0 1 --time 100 2>&1 | grep -E "induced|TIMING|GPU solve|gpu profile" | cut -c1-700 >> /tmp/bg.out; done
echo finished >> /tmp/bg.out
