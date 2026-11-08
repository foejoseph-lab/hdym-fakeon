# matter_loop — one-loop matter on a classical background (rung 1 of the ladder)

Concept folder: everything that computes the one-loop effective action of quantized matter on a classical
metric lives here, whatever the year or technique — Utiyama–DeWitt 1962 (cutoff), Schwinger proper time,
DeWitt's 1965 recursion, 't Hooft–Veltman 1974 §3 (dimensional regularization), Salvio 2018's matter
terms. Techniques accumulate inside; no sibling folders per paper.

**Asks:** redo Utiyama & DeWitt, J. Math. Phys. 3, 608 (1962) — the one-loop divergences of a quantized
scalar on a classical gravitational background, and the three counterterms they force (cosmological, Newton,
quadratic-in-Riemann) — with the technique of 1962 (momentum cutoff), then Schwinger proper time (1951, on the
perturbative graphs), then 1972–74 (dimensional regularization), then DeWitt's 1965 covariant recursion on a generic
background, and compare. The paper itself is paywalled and unseen; the computation is
reconstructed from its abstract and scored against 't Hooft & Veltman 1974 eq. (3.40) and Salvio 2018
eqs. 3.15–3.16. Setup: `../notes/2026-10-31_rung1_utiyama_dewitt_setup.md`. Results and findings:
`../notes/2026-10-31_rung1_utiyama_dewitt_results.md`.

**Run** (sympy; 1–8 min each, single-threaded, exact rationals):

    python3 rung1_run.py --regulator cutoff     --measure plain                  --log results/rung1_cutoff_plain.log
    python3 rung1_run.py --regulator cutoff     --measure covariant              --log results/rung1_cutoff_covariant.log
    python3 rung1_run.py --regulator propertime --measure plain     --heatkernel --log results/rung1_propertime_plain.log
    python3 rung1_run.py --regulator propertime --measure covariant --heatkernel --log results/rung1_propertime_covariant.log
    python3 rung1_dimreg.py --measure plain      --log results/rung1_dimreg_plain.log
    python3 rung1_dimreg.py --measure covariant  --log results/rung1_dimreg_covariant.log
    python3 rung1_dewitt.py --samples 6 --jobs 6             --log results/rung1_dewitt_samples.log
    python3 rung1_dewitt.py --samples 6 --jobs 6 --quartic 0 --log results/rung1_dewitt_quartic0.log

Added 1 Nov 2026 (the three items that closed rung 1, four runs; single-threaded sympy, no `--jobs`; predicted wall clock on the
author's machine in brackets, measured on the author's machine 1 Nov: 234 s, 301 s, 154 s, 149 s):

    python3 rung1_run.py --regulator cutoff --measure plain     --placement unshifted --log results/rung1_cutoff_plain_unshifted.log       [3-5 min]
    python3 rung1_run.py --regulator cutoff --measure covariant --placement unshifted --log results/rung1_cutoff_covariant_unshifted.log   [3-5 min; measured 301 s]
    python3 rung1_finite.py --measure plain      --log results/rung1_finite_plain.log       [2-4 min]
    python3 rung1_finite.py --measure covariant  --log results/rung1_finite_covariant.log   [2-4 min]

`--placement unshifted` puts the hard cutoff on the propagator momentum k itself (no Feynman parameter; large-k
expansion, exact in Λ⁴, Λ², log Λ²) and prints the difference from the shifted placement per order. `rung1_finite.py`
is the finite p⁴ log p² part of the massless matter loop under all three regulators (dim reg with n-dimensional
angular averages), its spin-2/spin-0 decomposition, the linearized Ward-identity test of the regulator constant, and
the 1/p⁴ resummation with the tree term derived from (√g R)⁽²⁾. Predictions are in both headers; results note §8.
END statuses of the committed logs: `LOG-OK_POWER-ORDERS-NONCOVARIANT` for both placement runs;
`LOG-OK_FINITE-NONTRANSVERSE[cutoff,propertime]` for the plain-measure finite run and
`LOG-OK_FINITE-NONTRANSVERSE[cutoff]` for the covariant one (the cutoff's finite part is not transverse with either
measure; proper time's is transverse only with DeWitt's measure; dim reg's always).

Added 2 Nov 2026 (the back pass to Sakharov 1967: the massive orders; two runs; single-threaded sympy; predicted wall clock on the
author's machine in brackets, from the sandbox's 426 s and 796 s at 2 cores and the 1 Nov ratio):

    python3 rung1_massive.py --measure plain      --log results/rung1_massive_plain.log       [6-8 min]
    python3 rung1_massive.py --measure covariant  --log results/rung1_massive_covariant.log   [11-14 min]

`rung1_massive.py` is the same graphs with m ≠ 0, expanded about p = 0: the p⁰ (m⁴√g) and p² (m²R) orders under all three
regulators, each as Λ⁴ + Λ² + log Λ² + log m² + constant, the constant fitted to the covariant density and Ward-tested, the Λ² term
beside the m² log term (the "−1/12 discipline"); the induced Planck mass printed per regulator. Predictions P0–P5 in the header
(heat-kernel/MS-bar values for the constants); results note §9. Expected END statuses: `LOG-OK_FINITE-NONTRANSVERSE[cutoff,propertime]`
(plain) and `LOG-OK_FINITE-NONTRANSVERSE[cutoff]` (covariant). The dim-reg tadpole now carries the k_⊥² of the flat trace term
(`ud1962_lib.perp_symbol`; results note §9 for the slip it fixed, kept in `results/superseded/rung1_massive_plain_sandbox_4d_tadpole_symbol.log`).

Added 4 Nov 2026 (rung 3, step 1, method M — 't Hooft–Veltman 1974 eq. (3.40) as code on bundles; summary and plan
`../notes/2026-11-04_rung3_thooft_veltman_summary.md`; one run, seconds):

    python3 hv1974_master.py --samples 8 --log results/hv1974_master.log

`hv1974_master.py` evaluates (3.40) pointwise on random algebraic curvature tensors for a scalar (any ξ), a vector in Feynman gauge
and Maxwell (vector minus two ghosts), fits R², Ric², Riem² and reduces by Gauss–Bonnet. Standard curvature sign ('t Hooft–Veltman's
R is minus the standard one; header). Scalar reproduces (5.23) and rung 1's ξ-polynomial; Maxwell comes out conformal, (1/10)(Ric² − R²/3),
twelve times the conformal scalar. Review of 4 Nov folded in. Method P (same evening):

    python3 rung3_dewitt_bundle.py --samples 4 --jobs 4 --log results/rung3_dewitt_bundle.log      [~16 min wall at sandbox speed; 4 cores]

(8 Nov: the script body is under a main guard — the 4 Nov `--jobs 4` run forked forever on Windows — a non-empty log is refused, and
a line is printed per finished background.) `rung3_dewitt_bundle.py` is DeWitt's recursion on a bundle (matrix a_k, the connection in the transport equation so a₀ is the parallel
propagator, the endomorphism as a field), rung 1's normal-coordinate machinery copied with credit. One-background smoke in the sandbox:
Tr a₂ equals rung 1's ξ-polynomial (scalar) and the Gilkey coefficients (vector) exactly; the rule-12 check is that half of it differs
from (3.40) by exactly (d/360)·Gauss–Bonnet and agrees after. Fresh review folded in; what the backgrounds cannot see is in the header.

Step 2 (8 Nov): the same script carries the graviton bundle (10-dimensional, the induced connection, the endomorphism derived in
`../graviton_loop/` as a field) and 't Hooft–Veltman's ghost; predictions from method M by a₂ = 2M + (d/180)(1, −4, 1); slow — an
overnight on the author's machine:

    python3 rung3_dewitt_bundle.py --samples 4 --jobs 4 --bundles graviton,ghost --log results/rung3_dewitt_graviton.log

`--measure plain` integrates over φ (operator K = √g·D); `--measure covariant` over χ = g^{1/4}φ (DeWitt's
covariant measure). `--heatkernel` derives a₀, a₂, a₄ at second order in h from Tr e^{−tD} and compares with
Vassilevich hep-th/0306138 eqs. 4.26–4.28.

`rung1_dewitt.py` is technique (b2): DeWitt's recursion for the heat-kernel coefficients, run in Riemann normal
coordinates on random generic backgrounds (exact rationals, no expansion in h; `--jobs` over the samples; the
`--quartic 0` run is a control, see the results note §7). Setup and predictions:
`../notes/2026-10-31_rung1_b2_dewitt_recursion_setup.md`.

**Files.** `ud1962_lib.py` — conventions (header), vertices, angular averages, the three regulators'
integrals, metric expansion to second order, the tadpole as an operator symbol, fitting; dated block of 1 Nov: the
unshifted-cutoff bubble, Barnes–Rivers projector forms, the massless finite part with n-dimensional dot products.
`rung1_run.py` — cutoff (both placements) and proper time. `rung1_dimreg.py` — dimensional regularization, pole only.
`rung1_finite.py` — finite p⁴ part, Ward test, resummation. `rung1_massive.py` — the massive p⁰, p² orders (2 Nov). `rung1_dewitt.py` — the covariant
recursion. `results/` — the six graph logs, the two recursion logs, the two placement logs and the two finite-part logs (and `results/superseded/`, the development scripts and logs, kept because one of
them records a prediction that failed — see the results note §4.1); each ends
`# END <name> elapsed=<s>s status=...` with status `OK`, `LOG-OK_POWER-ORDERS-NONCOVARIANT` (the expected
outcome for a hard cutoff, and for proper time with the plain measure), or `MISMATCH`.

**Outcome in one line.** The log-divergent (quadratic-curvature) coefficients come out identical six ways on
the graphs and equal 't Hooft–Veltman (3.40); the covariant recursion on generic backgrounds reproduces them
and adds the Riem² coefficient the graphs cannot see (Euler term −1/360); the ∞⁴ and ∞² assignments of the
abstract are reproduced only by proper time with the covariant measure — the Feynman-parametrised hard cutoff
does not give them with either cutoff placement (the placement changes only the Λ² term, by (1/576)p²(h² + 2h_ijh_ij)),
and the plain measure misses ∞⁴ by exactly ½(ln√g)⁽²⁾ (units Λ⁴/16π²). The finite p⁴ constant (1 Nov): transverse in dim reg, where
it equals −46/15 on Ric²/120 [v 2 Nov: the massless limit of the Gorbar–Shapiro form factor, hep-ph/0210388 eq. 3.7], and, up to Euler's γ, in proper time
with DeWitt's measure; not transverse with the plain measure or with a hard cutoff. Resummed, the spin-2 inverse propagator is
(M_P²/4)p² + (1/3840π²)p⁴ log(p²/Λ²) + …: the abstract's 1/p⁴ with a running coefficient (results note §8.3).
