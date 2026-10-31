# Route 1, rung 1: the fakeon scalar bubble against the published numbers (30 Oct 2026, night)

Working notes for `loops/scalar/`. Status: rung 1 of brief 11 §5b done — the one-loop bubble of scalar fields with
the fakeon (average continuation), Lee–Wick and Feynman–Wheeler prescriptions computed by two independent methods
in D = 2, 3, 4, and the closed forms in 1703.04584, 1801.00915 and 2001.01942 that apply to it reproduced ((5.12)
through route (b) only); the published values certified in ball arithmetic. Reviewed once by a reader who had not
seen the work; findings folded in (History).
Nothing new is claimed here; this is the check the ladder stands on.

## Question

The perturbative route needs two things before the h/χ/φ multiplet: a way to compute amplitudes with the average
continuation that is not a transcription of the papers, and a second way that agrees with it. Does the bubble of
two scalars — the integral every later rung is built from — come out right both ways, to what precision, and where
do the two methods have to be handled with care?

## Answer

**Both methods reproduce every published number, to 10⁻¹² or better where the papers give closed forms, and
agree with each other to 10⁻¹⁴ or better everywhere (10⁻¹⁶–10⁻²⁰ away from the LW threshold; 2·10⁻¹⁴ at D = 2,
t = 2.1).** Which comparisons are independent evidence and which are consistency of one reading is said below. The numbers (M = 1 or m = 1; `results/*.log`):

| published | what | hit to |
|---|---|---|
| 1801.00915 (5.4) | fakeon of equal masses, D = 4: f_AV = ½∫ln[(1 − zx(1−x))²] | 10⁻²⁰ (exact closed form, certified) |
| 2001.01942 (6.6) | Σ_f − Σ = √(u₊u₋)/(16πp²) θ(−u₊), D = 4, m₁ = m₂ (three points above threshold) and m₁ ≠ m₂ (two) | 10⁻²⁰ |
| 2001.01942 (4.1) | 2 Im(−iΣ) = √(u₊u₋)/(8πp²) θ(−u₊) | 10⁻²⁰ |
| 1801.00915 (6.3) | fakeon propagator = lim_{E→0} ½[1/(p²−m²+iE²) + 1/(p²−m²−iE²)] | converges as E² (D = 4), E⁴ (D = 3), coefficient stable to 3 digits |
| 1801.00915 (5.6)–(5.8) | Lee–Wick pair, D = 4: f_AV = f₀ + (2π/t)(√(4−t²) − 2) for 0 < t < 2, f₀ − 4π/t for t > 2; g | 10⁻²¹ below the LW threshold, 10⁻¹⁴ across and above it (the 10⁻¹⁴ end-point offset of the continuation, amplified by the √(t − 2) behaviour) |
| 1801.00915 §5, D = 2 | f_AV = f₀ + 4π/√(4−t²) for 0 < t < 2, f₀ for t > 2 | 10⁻²¹ below the threshold, 10⁻¹² at t = 1.9 and 2.1 (the (t − 2)^{−1/2} behaviour amplifies the end-point offset more than the √ in D = 4) |
| 1801.00915 (5.12) | Lee–Wick pair, D = 3: M_AV ∝ Re[f₀ + g] with the Euclidean closed forms continued | 10⁻¹⁴ (route (b) against the continued elementary functions); the Euclidean forms themselves to 10⁻²⁰ |
| 2001.01942 (2.6) | Feynman–Wheeler nonlocal divergence i lnΛ² (m₁² − m₂²)/((4π)²p²) | coefficient to 4·10⁻⁷ relative at Λ = 500–5000; the deviation falls by 100.0 per decade of Λ (1/Λ²) |
| 2001.01942 (4.3) | 2 Im(−iΣ′) = −√(u₊u₋)/(8πp²) θ(u₋): the pseudothreshold made a threshold, with the wrong sign | 10⁻²⁰ |
| 2001.01942 (3.1) | Σ′ for m₂ = 0 with a cutoff, including its iπ for p² < m² | 10/Λ² (1.2·10⁻⁵ at Λ = 200, 1.2·10⁻⁷ at 2000) — which also fixes their cutoff as |k| ≤ Λ |
| 2001.01942 §7 | D = 3, 2 imaginary parts of Σ′ | 10⁻²⁰ |
| 1703.04584 §3, fig. 13 | propagator (2.1), m = 3, M = 1: Im J continuous with a cusp at p² = 2M² (left slope ∝ (2 − t)^{−1/2}, right slope finite, smooth elsewhere to 10⁻⁷); J⁺ − J⁻ above the threshold equals their rigid-contour residue (1/16π) M⁴/(m⁴+M⁴) √(1 − 4M⁴/(p²)²) | ratio of the half-difference to that residue 0.500000000 at p⁰ = 1.6, 2.0, 2.5 — an algebraic identity once the D = 4 f_d is known, so this certifies the transcription of their residue, not a new fact; Re J = 0 on the grid is a reflection identity of the averaged terms (p⁰ < 2m), listed as a consistency check |

Independent evidence versus consistency: route (a) against route (b) is independent (different integrals, no
shared branch bookkeeping) and holds in every row; route (a) against a closed form is a check of one reading of the
paper, coded once in the generator and once in the checker (a shared mis-transcription would pass both — but would
then break the (a) = (b) agreement, since (b) never sees the closed form). The f_d values above the threshold are
compared with *my continuation* of the published sub-threshold terms, not with a printed formula; they are labelled
so in the logs. Checks that cannot fail (Schwarz reflection of a transcribed expression, Im(B_PV − B_F) = 0 by
construction) are labelled as such.

What f_d is, by dimension, above the LW threshold t > 2 (f₊ − f₋ = 2 f_d; both methods; the continuation of the
published sub-threshold terms): D = 4: f_d = ±(2πi/t)√(t² − 4); D = 3: Im f₀ = ±2π/√t; D = 2: f_d = ±4πi/√(t² − 4). In every
dimension the average kills it and f_AV is real — that is the fakeon: no discontinuity, no imaginary part, the
threshold still visible as a cusp. Below the threshold f₊ = f₋ to 10⁻¹⁴ in all three dimensions (the LW
threshold is a branch point only on the second sheet; the two ways around t = 0 land on the same function).

## What had to be learned to make the two methods agree

1. **The continuation path must stay in the main region.** The (++) and (−−) Lee–Wick terms, with Δ = iM² −
   tx(1−x), have branch points at t = ±4iM² (1703.04584: the cuts xy = 2M², x² ≥ y² in the p⁰ plane). A
   continuation from t < 0 to t > 0 along |t| = T crosses the imaginary axis above 4i for T > 4 and lands on the
   wrong sheet of g (seen: g_AV wrong by 1.5 at T = 8, with a spurious g_d). The path used is the arc |t| = 1
   around t = 0, then the real axis. With it g_AV = g to 10⁻¹⁷ at all T. The (+−) terms have branch points only at
   t = 0 and, on the second sheet, t = 2: either arc gives the same f.
2. **The LW threshold is stepped around, not through.** On the second sheet f ∝ √(t − 2) (D = 4) or
   (t − 2)^{−1/2} (D = 2) at t = 2. Walking along Im t = 10⁻¹⁴ across it, the continuity tracker that picks the
   monodromy integer by linear extrapolation faces a near-tie (best/second-best distance ratio 0.7–1.0) and in
   D = 2 it chose the wrong integer in an early run — giving f_d = 0 above the threshold while f_AV came out
   right, because the two wrong choices on the two paths average to the right thing. The momentum-space route,
   which only sees f_AV, could not catch it. A semicircle around t = 2 (radius min(0.1, half the distance to the
   nearest requested T) = 0.05 in these runs) makes the worst ratio 0.04 (D = 4) and 0.10 (D = 2) and gives the f_d
   above. In D = 2 the two residues are equal and opposite, so the integer labels are degenerate up to (n, n) and
   the logged jump counts there mean nothing; the logs say so. Lesson for later rungs: f_AV agreeing between two methods
   does not check f_d; f_d must be checked on its own, against Schwarz reflection (f₋ = f₊*) at least.
3. **Principal-value windows.** Tanh-sinh quadrature puts nodes within 10⁻ᵈᵖˢ of a breakpoint; a subtracted
   integrand h − R/(k − k*) then cancels catastrophically and the PV came out wrong by a factor 2 before the
   window around k* was integrated by Gauss–Legendre on its two halves (nodes O(width/n²) from the ends).
4. **Factored Δ next to its zeros.** Above threshold Δ = T(x − x₁)(x − x₂); the unfactored form loses all digits
   within 10⁻¹⁰ of a root and spoiled the D = 3 integrals at 10⁻¹². For ordinary masses everything above threshold
   is elementary anyway (ln T + Σ[(1 − x_j)ln(1 − x_j) + x_j ln x_j − 1] ∓ iπ(x₂ − x₁) in D = 4, asinh's in D = 3,
   a log and residues in D = 2), and that is what both the generator and the checker now use.
5. **Non-conjugate complex masses have no real pinch** but a complex pole near the real k axis (the (++) term of a
   nearly real LW pair, distance ∼E²): without a breakpoint there, the E → 0 check broke below E = 0.1.

## Numerics

mpmath at 20 digits; the generator's end-point offset for the continuation is δ = 10⁻¹⁴ (error O(δ), O(√δ) at
the threshold). Route (a) for the Lee–Wick pair: four sweeps (two mass orderings × two half planes) per dimension
over 342 path points each, 82 s total for D = 4 and 2. Route (b): seconds. The checker (`bubble_cert.py`,
python-flint 0.9, 40 digits) encloses the published closed forms to 10⁻³⁹–10⁻⁴⁰ by rigorous integration (the
integrands are holomorphic near [0, 1] for real t ≠ 0) and compares 152 numbers; all pass, worst gap 1.25·10⁻⁵ at
the declared 10/Λ² tolerance of (3.1), everything exact to 10⁻¹¹ or better. `--dps 30` is available and was used
once to confirm that the 10⁻¹² residuals are the generator's and not the physics'.

## Open

- Route (b) is in the rest frame only. The region structure of 1703.04584 (the point P′, Lorentz invariance
  recovered region by region) needs p⃗ ≠ 0 and a two-dimensional pinch surface; the PV then lives on a curve in
  (|k|, cos θ). Not needed for rung 2's self-energies at p⃗ = 0, needed if the paper wants the frame-dependence
  statement checked. Estimated as a two-dimensional version of `B_momentum` with the ±iδ averaging instead of an
  explicit PV.
- 1703.05563 (the cutting equations) is reproduced here only through its consequences — the cut propagator of
  the fakeon vanishes (Im B_PV = 0 identically) and the Feynman–Wheeler one has the wrong sign — not through the
  algebraic cutting identities themselves.
- The D = 3 Lee–Wick continuation is done on route (b) only; the square-root monodromy in route (a) was not
  implemented because route (b) plus the continued closed forms already give two independent checks.

## For rung 2

The two routes carry over as they are: route (a) with the tensor numerators reduced to scalar integrals, route (b)
with the same residue-plus-PV structure. The hdym `ym/hdym_derive.py` Lagrangian and conventions are to be
referenced for the Yang–Mills step, not copied. The toolchain stands: mpmath generator, python-flint checker, no
Mathematica. Reading before rung 2: Anselmi–Piva 1703.05563 in full for the cutting equations, Anselmi 2207.11271
for the purely virtual Yang–Mills extension.

## History

- 30 Oct 2026, night: built and run. Three numerical defects found and fixed in the process (items 2–4 above),
  one path error (item 1); all left in the code as comments where they bit.
- 30 Oct 2026, later: independent review (a reader who had not seen the work; it re-derived the residue formula,
  the Feynman sign, the prefactors and the monodromies by hand and recomputed a dozen of the numbers). No blockers.
  Findings folded in: the headline precisions were overstated (10⁻¹⁴ at the D = 2 threshold, not 10⁻¹⁶; 10⁻¹² for
  the D = 2 closed form across the threshold), the detour radius was misreported (0.05, not 0.1), the (6.6)
  m₁ ≠ m₂ comparison rested on one point (a second added), the fig. 13 "ratio 0.5" and "Re J = 0" were counted as
  reproductions when they are identities given the other results, the f_d comparisons were labelled "published"
  when they are continuations of published terms, the continuation path had the LW-pair geometry hard-coded (now
  derived from the masses, with the arc radius checked against the (a, a) branch points), the D = 2 jump counts
  were undocumented, the 1/Λ² approach to (2.6) was claimed but not logged (now printed and checked), the D = 3
  (5.12) tolerance was 100× looser than the data (now 10×). Code and logs regenerated; all PASS; 163 certified.
