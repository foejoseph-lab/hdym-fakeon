# bubble/ — the one-loop bubble of a fakeon scalar, two ways

The instrument check, and the p⃗ = 0 case of rung 2 on the 31 Oct ladder. The question asked here is a check, not a result: can the
one-loop bubble of scalar fields be computed, with the fakeon prescription (Anselmi's average continuation),
the Lee–Wick prescription and the principal-value (FW) prescription, by two independent methods that agree with each
other and with the numbers Anselmi and Piva published? Everything downstream (QED, Yang–Mills, the h/χ/φ
multiplet of quadratic gravity) reuses these two methods, so they are checked here on the simplest integral, in
D = 2, 3, 4, against every closed form 1801.00915, 2001.01942 and 1703.04584 give. Results and their reading:
`notes/2026-10-30_rung1_fakeon_scalar_bubble.md`.

The two methods, for B_D(p²) = ∫d^Dk/(2π)^D 1/((k² − m₁²)((k − p)² − m₂²)) with complex masses allowed:

- **Method P — Feynman parameters, continued.** B_D ∝ ∫₀¹dx Δ^{(D−4)/2}, Δ = m₁²(1−x) + m₂²x − p²x(1−x), is analytic in the
  Euclidean region p² < 0 and is continued to p² > 0 along a path in the complex p² plane; the roots of Δ cross
  the integration segment on the way and the monodromy of the integrand (2πi(1−x_j) for the log, 2πi Res_j for the
  pole) is added, the integers tracked by continuity. f₊ and f₋ are the two paths; f_AV = (f₊ + f₋)/2 is the
  fakeon/Lee–Wick amplitude, f_d = (f₊ − f₋)/2 the would-be discontinuity (1801.00915 §4).
- **Method M — momentum space, nonanalytic Wick rotation.** In the rest frame the k⁰ integral is done by residues; the
  |k| integral has a pinch where p⁰ = ω₁(k) + ω₂(k) for real k, and there the fakeon/Lee–Wick prescription is the
  principal value (1703.04584 eq. 3.11), the Feynman prescription adds iπ times the residue, and the
  principal-value (FW) mixed term has its pinch at the pseudothreshold instead (2001.01942 eq. 2.3).

## Scripts

| script | what it does | outputs in `results/` |
|---|---|---|
| `bubble_checks.py --fakeon` | ordinary masses: f_AV against 1801.00915 (5.4), B_PV − B_F against 2001.01942 (6.6), Im B_F against (4.1), methods P = M in D = 2, 3, 4, the fakeon as the E → 0 limit of a Lee–Wick pair (6.3) | `fakeon.log`, `fakeon.json` |
| `bubble_checks.py --lw` | the Lee–Wick pair m² = ±iM²: (5.8) in D = 4, the D = 2 closed form, (5.12) in D = 3, f_d above the LW threshold, f₊ = f₋ below it, methods P = M | `lw.log`, `lw.json` |
| `bubble_checks.py --fw` | principal-value (FW) mixed term: the nonlocal divergence (2.6), the pseudothreshold imaginary part (4.3), the finite part (3.1), the D = 3, 2 imaginary parts of §7 | `fw.log`, `fw.json` |
| `bubble_checks.py --fig13` | 1703.04584 propagator (2.1) with m = 3, M = 1: Im J(p⁰) across the LW threshold, Re J = 0, the cusp at p⁰ = √2, the half-difference of the two continuations against the rigid-contour residue of their §3 | `fig13.log`, `fig13.json` |
| `bubble_frame.py` (2 Nov; rung 2, technique of the day) | the complex-conjugate-pair bubble at p⃗ ≠ 0 with the Lee–Wick energy contour and real spatial loop momenta — Nakanishi 1971's computation (restated in hep-th/0609206 eqs. 3.3–3.6): the boosted-frame method M against the rest frame at the same s, the moving threshold s_min(|p⃗|) = 4(Re E_{p⃗/2})² − p⃗², continuity across it, the |p⃗| → 0 limit against the three rest-frame prescriptions, method P as the anchor; predictions P1–P5 in the header | `frame.log`, `frame.json` |
| `bubble_cert.py` | trusted checker: the published closed forms in ball arithmetic (Arb via python-flint), compared with the generator's numbers in the `.json` files | `cert.log` |

**Spyder:**
```
import sys
%cd C:\GitHub\hdym-fakeon\bubble
!"{sys.executable}" -B -u bubble_checks.py --fakeon > results\fakeon.log
!"{sys.executable}" -B -u bubble_checks.py --lw     > results\lw.log
!"{sys.executable}" -B -u bubble_checks.py --fw     > results\fw.log
!"{sys.executable}" -B -u bubble_checks.py --fig13  > results\fig13.log
!"{sys.executable}" -B -u bubble_cert.py            > results\cert.log
!"{sys.executable}" -B -u bubble_frame.py --json results/frame.json > results\frame.log     [2 Nov; 8-12 min predicted from 1131 s on the 2-core sandbox]
```

Wall clock on one core: `--lw` about 90 s (four continuation sweeps per dimension, two dimensions), the rest
seconds; each log ends with its elapsed seconds and `status=PASS|FAIL`. Dependencies: mpmath (generator),
python-flint (checker; `pip install python-flint`). Nothing here reads from hdym. Precision knob: `--dps 30`
(default 20; the end-point offset of the continuation is 10^−(dps−6)).

## Conventions

Metric (+,−,−,−), propagator 1/(k² − m²), p = (p⁰, 0), t = p² in units of M² (Lee–Wick) or m² (ordinary masses).
Published formulas are quoted by equation number from 1703.04584, 1801.00915 and 2001.01942; the reduced
integrals f₀, g, f_AV of 1801.00915 §5 are what the logs print, and B or J (with the partial-fraction weights of
the propagator) where a paper prints those. In D = 4 everything is UV-subtracted at the Euclidean point t = −1 in
both methods, so only differences are compared; D = 3 and D = 2 are finite.
