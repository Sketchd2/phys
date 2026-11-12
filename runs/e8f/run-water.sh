#!/bin/sh
# Water first: every close liquid-drawn pair and the wall pairs, two streams.
B=../../target/release/phys-overlap
( for every in 2 1; do $B water ../e8c/avg/pairs-water-gpu-bias-mp2.txt --max-sep 7 --every $every --out overlap-water-bias.txt; done; echo bias finished ) > overlap-bias.log 2>&1 &
( for every in 4 2 1; do $B water ../e8c/wall-rigid/pairs-water-gpu-wall-mp2.txt --every $every --out overlap-water-wall.txt; done; echo wall finished ) > overlap-wall.log 2>&1 &
wait
./run-overlap.sh > overlap.log 2>&1
