#!/bin/bash
# G1: MP2 pairs drawn from the wall-law liquid (population-weighted), then their radial slopes.
# Waits for the virial job to leave the GPU.
cd /e/Projects/phys/runs/e8d/g1
while tasklist | grep -q phys-cluster-gpu; do sleep 30; done
. /e/Projects/phys/runs/many/env.sh
B=E:/Projects/phys/target/release
$B/phys-pairs-liquid-gpu.exe water wall-liquid.snap 100 --mp2 >> pairs.log 2>&1
$B/phys-es-gpu.exe water pairs-water-liq-gpu-mp2.txt --deriv --max-oo 5.2 >> slopes.log 2>&1
