# Many-body probe (E8d) — working notes, 2026-10-14

Classical (phys-manybody, law runs/e8c/law-bias.txt, 216 centres, nearest-N clusters of bulk-water-298.snap):
see classical-all.txt. Law's induction non-additivity ~0 at every size 3-10.
Quantum: phys-cluster-gpu --method hf (HF only; new flag) — timings in t4.log.

Quantum route: 4-mer HF whole-cluster field alone = 2348 s (disk-bound), so direct clusters >4 infeasible.
HF trimer = 236 s cluster + ~46 s/pair; reproduces stored MP2 file HF numbers (3-body -0.0065 kcal/mol centre 0).
Queue (queue.sh): phys-cluster-gpu --manybody --method hf; centres 0..210 to size 6, then centre 0 to 10.
Files: mb-hf-c<centre>.txt (resumable), logs mb-hf-c<centre>.log.
