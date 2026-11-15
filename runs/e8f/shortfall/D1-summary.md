# D1: the basis bracket at the S22 dimer minimum (kcal/mol), gathered by the fleet

Gap = uncorrected minus counterpoise total; "middle" is their mean. Accurate values are S22's CCSD(T)/CBS.

| Dimer | Basis | Functions | Counterpoise | Uncorrected | Gap | Gap / accurate | Middle vs accurate |
|---|---|---|---|---|---|---|---|
| methane | engine's | — | -0.475 | -0.660 | 0.185 | 35% | -0.568 (+7.1%) |
| methane | + g on C, f on H (A) | 714 | -0.478 | -0.589 | 0.111 | 21% | -0.533 (+0.7%) |
| methane | + A and a second g/f, f on C, d on H (AE) | 892 | -0.483 | -0.538 | 0.055 | 10.3% | -0.510 (-3.7%) |
| ammonia | engine's | — | -3.049 | -3.236 | 0.187 | 5.9% | -3.143 (-0.9%) |
| ammonia | A | 628 | -3.074 | -3.162 | 0.089 | 2.8% | -3.118 (-1.7%) |
| ammonia | AE | 774 | -3.104 | -3.161 | 0.057 | 1.8% | -3.133 (-1.2%) |
| water | ADE (HANDOFF section 5) | 648 | -4.870 | -5.021 | 0.151 | 3.0% | -4.946 (-1.5%) |

Exponents for methane and ammonia were scaled from water's by heavy atom (shortfall/more2.sh).
Against D1's done-when (water under 2% wide and middle within 2%; methane within 5%): ammonia's bracket
closes to 1.8% with its middle 1.2% off; methane's is still 10% wide but its middle is 3.7% off, inside 5%;
water's last measured step is 3.0% wide. As with water, the uncorrected end settles first (methane -0.589 -> -0.538,
ammonia -3.162 -> -3.161) while the counterpoise end keeps climbing.
