# Rung 2 results, part 2 (3 Nov 2026): (v) the control — the real-axis k⁰ integral of the complex-pair bubble

The control for rung 2's forward pass: the bubble of the complex-conjugate pair computed the way anyone would have written it before
Lee–Wick's contour existed — loop energy along the real k⁰ axis, real spatial momenta, a hard cutoff on |k⃗| in the frame of the
computation. In 2017 this is the "Minkowski formulation" of Aglietti–Anselmi 1612.06510 (EPJC 77, 84). Code and logs: `bubble/bubble_minkowski.py`,
`bubble/results/minkowski.{log,json}` (the author's re-run of 3 Nov evening with P7 and the review's window fix, 277 s, dps 15, commit 3c0ba46; identical to the
sandbox full run on every printed digit; sandbox logs are not committed, rule 15).
Setup and the two techniques before it: `notes/2026-11-02_rung2_1971_crack_results.md`.

**Rule 28 (brief 16 amendment, 2 Nov 23:46), first use.** Nothing from 1612.06510 or later was fetched before the run. The header's
predictions are marked [derived] (pole positions and the large-k expansion, sympy, 2 Nov 23:50–23:58), [obvious to a 1969 reader], or
[memory, tainted] (the abstract, read 31 Oct). The paper was read after the author's log had been compared with the sandbox smoke
(3 Nov 17:45) and the comparison is §4 below and in the header's AFTER THE READ block.

Model: a₁ = i, a₂ = −i (M = 1), D = 4, metric (+,−,−,−), B(p) = ∫d⁴k/(2π)⁴ 1/((k² − a₁)((k − p)² − a₂)); rest frame p = (p₀, 0) and
boosted frame |p⃗| = 0.8 at the same s = p₀² − p⃗²; s ∈ {1.0, 2.5}; cutoff Λ ∈ {10, 30, 100, 300} on the sphere |k⃗| < Λ in the frame;
subtraction at the real point s_ref = 0.5 (§2, the Euclidean point is not available).

## 1. The k⁰ integral: what the real axis selects

For real p₀ the four poles k⁰ = ±ω₁, p₀ ± ω₂ (ω_i = √(k⃗² + a_i), Re > 0) are off the real axis for every real k⃗: Im ω₁ > 0, Im ω₂ < 0.
Closing in the upper half plane selects +ω₁ and p₀ − ω₂ — selected by sign(Im) in the script, generically — and the residue sum is

    h₊₋ = −(i/2) (ω₁ − ω₂) / (ω₁ ω₂ [p₀² − (ω₁ − ω₂)²]),

with **no positive-energy term** (ω₁ + ω₂)/(…) at all: the real-axis contour treats a = +i as an anti-Feynman (−iε) propagator and a = −i
as a Feynman one, and h₊₋ is the mixed term of one +iε and one −iε propagator of 2001.01942 eq. 2.3 (`bubble_lib.h_fw_mixed`, written
there for the (−i, +i) ordering: the sign of ω₁ − ω₂ flips with the ordering) with complex masses [derived; P5: residue sum = direct real-axis quadrature at four sample k⃗, both frames, to 2·10⁻¹⁶]. For (−i, −i) the upper half plane holds −ω₁ and p₀ − ω₂ and the
real-axis value is the Lee–Wick (contour) value; for (+i, +i) it holds +ω₁ and p₀ + ω₂ and the value is **minus** that [derived; P4: both to 0 at Λ = 30, s = 2.5,
against `bubble_lib.B_momentum(kmax = Λ)`]. The denominator p₀² − (ω₁ − ω₂)² stays O(1) at large k instead of O(k²): this is the
"falls off half less rapidly" of 1612.06510 §5, seen before reading it.

## 2. The divergence and its coefficient

Large-k tail in the rest frame (u = 1/k): h₊₋ = u³(p₀² − u²)/(2p₀⁴) + O(u⁷), so

    B_Λ(s) = ln Λ / (4π² s) + C(s, |p⃗|) + O(Λ⁻²)        [derived: the 1/s; the boosted-frame form corrected after the smoke, §3]

— a logarithmic divergence with a **nonlocal** coefficient 1/s, which no counterterm absorbs, **real in B** (the ordinary local log of a bubble is
imaginary in B, +i/(8π²) ln Λ: the two sit in different parts; in 1612.06510's language the ordinary one is anti-Hermitian and this one is not). The Euclidean subtraction of `bubble_frame.py` cannot
remove it, and in fact cannot be made: with p₀ = i the sign(Im) selection is not the continuation of the real-p₀ one (three poles move to
the upper half plane) and B_ref grew like Λ in the first smoke — **the real-axis integral has no Euclidean region to subtract at**
(script error caught by the smoke; 1612.06510 §5 says the same: "the analytic continuation to Euclidean space is not possible").
Subtraction at a real s_ref = 0.5 leaves (1/4π²)(1/s − 1/s_ref) ln Λ.

| | rest, s = 1.0 | boost, s = 1.0 | rest, s = 2.5 | boost, s = 2.5 |
|---|---|---|---|---|
| slope dB/d ln Λ, last pair (100→300) | 0.0253293 | 0.0253287 + 6.6·10⁻⁷ i | 0.0101320 | 0.0101320 + 1.0·10⁻⁷ i |
| predicted 1/(4π² s) | 0.0253303 | same | 0.0101321 | same |
| rel. error | 4.0·10⁻⁵ | 6.7·10⁻⁵ | 1.6·10⁻⁵ | 1.4·10⁻⁵ |
| slope of B(s) − B(s_ref), rel. error vs (1/4π²)(1/s − 1/s_ref) | 1.2·10⁻⁴ | 3.6·10⁻⁴ | 9.7·10⁻⁵ | 2.6·10⁻⁴ |
| C = lim [B_Λ − ln Λ/(4π²s)] | −0.0068343 | −0.0112123 + 0.0021890 i | −0.00047910 (the author's first log had −0.00047983: the window bug, review of 3 Nov) | −0.0012663 + 0.00098393 i |
| convergence of C (last two Λ) | 1.1·10⁻⁶ | 1.9·10⁻⁶ | 1.8·10⁻⁷ | 1.5·10⁻⁷ |

(`minkowski.log`, P1–P3; the O(Λ⁻²) next term is −u⁵/(2p₀⁴), and the slopes approach the prediction like Λ⁻².) P1 and P2 **hold** as re-scored.

## 3. Frames: where the Lorentz-invariance failure is, and where it is not

**Prediction D3 was wrong, and the smoke said so before the paper did.** The header predicted the ln Λ coefficient 1/(4π²(s + p⃗²)) in the
boosted frame — a frame-dependent divergence. The boosted slope is 1/(4π² s), the same as the rest frame (table above). The tail
expansion had dropped the O(1) term k̂·p⃗ in ω₁ − ω₂: with x = k̂·p⃗, ω₁ − ω₂ = x + (a₁ − a₂ − p⃗² + x²)/(2k) + …, the leading term is odd
in k̂ and averages out, and the 1/k term's angular average is (1/2)∫₋₁¹ dc (A + p²c²)/(A − p²c²)² = 1/(A − p²) = 1/s with A = p₀²
(the integrand is d/dc of c/(A − p²c²)). **The coefficient of the nonlocal divergence is frame independent** — under boosts, with the sphere cutoff, by the identity for any
|p⃗| < p₀ [exact] and by the numbers at two s and one |p⃗| [grid]; "Lorentz invariant" with a non-invariant cutoff domain is more than
was shown. The frame dependence is in the finite part C,
in both its real and imaginary parts, and larger at lower s (table). Kept as written in the header with the correction beside it (rule 11).

This is a different crack from Nakanishi's (part 1): that one was finite, linear in |p⃗|, in the imaginary part only, and supplied by the
geometry of a pinch circle; this one is in the finite part of a divergent integral whose divergence is covariant, because the integral
is only conditionally convergent (tail 1/k³ against k² dk) and the sphere |k⃗| < Λ is not a Lorentz-invariant domain [obvious to 1969: O1].

**Smooth in s** (P6, reported): on s ∈ {1.5, 1.7, 1.9, 2.1} at Λ = 30, |p⃗| = 0.8, the second differences are 1.07·10⁻³ and 7.9·10⁻⁴ (real),
3.0·10⁻⁵ and 2.2·10⁻⁵ (imaginary), one sign: nothing happens at s_min(0.8) = 1.707 or at 2M² = 2. No pole crosses the real k⁰ axis for
real p₀, so the real-axis formulation has no thresholds at all [O2].

## 4. The four-pair sum, and the partial-fraction cancellation (D4, P7)

Rest frame, s = 2.5, slopes dB/d ln Λ of each real-axis term (`minkowski.log`, D4):

| term | slope (100→300) | reading |
|---|---|---|
| (+i, +i) | −7.7·10⁻⁷ − 0.0126655 i | the **ordinary local** log, −i/(8π²) = −0.0126651 i, with the anti-Feynman sign |
| (−i, −i) | −7.7·10⁻⁷ + 0.0126655 i | the ordinary local log, +i/(8π²): the Lee–Wick value |
| (+i, −i), (−i, +i) | 0.0101320 each | the nonlocal 1/(4π² s), identical for the two orderings (k → p − k) |
| sum | 0.0202624 | predicted 2/(4π² s) = 0.0202642 |

In the symmetric four-pair sum **the local divergence cancels between (++) and (−−), and only the nonlocal one survives, doubled.** Under the
Lee–Wick contour the same four integrands give the opposite: every term carries the same local tail +i/(8π²) (sum 4i/(8π²), absorbable)
and no nonlocal piece (the tail of `bubble_lib.h_bubble` is i(4 − (3a₁ + 3a₂ − p₀²)u²)u³/16, no 1/p₀² term; `results/lw.log` holds only
the subtracted, finite quantities). One contour choice, and the divergence changes character rather than size. 1612.06510's "the usual divergences of local theories
are anti-Hermitian, these are not" is this fact from the counterterm side.

**P7 (added 3 Nov after the read, author's yes; scored).** 1612.06510's propagator is the physical pole times the pair, S = 1/(x − m² + iε) ·
M⁴/(x² + M⁴), x = p², which at m = 0 is 1/x − (1/2)/(x − i) − (1/2)/(x + i); their D = 4 trivial-numerator coefficient is c₀,₀ = 0 — no
nonlocal log. Our tail formula predicts that: for the (0, +i) pair the upper half plane holds −ω₁ and p₀ + ω₂ (the "(−+)" selection, the mirror of §1's),
whose 1/k³ coefficient is i(a₁ − a₂)/(4p₀²) = +1/(4p₀²) against the bare pair's 1/(2p₀²); the weights 2·(−1/2) and 1/2 give
(1/4)(−1) + (1/2)(1/2) = 0. Sandbox (dps 12 in the smoke, dps 15 in the full sandbox run, Λ = 30 → 300, s = 2.5, the massless Feynman pole
as a₁ = −10⁻¹⁰ i, i.e. 1/(k² + iε)): weighted log-slope −3.6·10⁻⁷ + 9.1·10⁻⁷ i (the O(Λ⁻²) remainder) against the bare pair's 0.0101311
— **PASS** (`minkowski.log`, P7). This is the m = 0 entry only; the m ≠ 0 entries of their table have m-dependent
weights and are not checked. So the bare Lee–Wick pair (the 1801.00915 model) keeps the nonlocal log; a Stelle-type propagator with the
physical pole cancels it at trivial numerator and needs derivative vertices (their c₀,₁ = 1, c₁,₁ = 3, …) to show it. Neither paper
says this in one sentence.

## 5. Against 1612.06510 (read 3 Nov 17:45 after the run, fetch summary of §§1–3, 5, app. A; §3 and footnote 1 read in full 19:10)

| ours | theirs | |
|---|---|---|
| formulation: real k⁰ axis, residues in the UHP, hard cutoff on |k⃗|, rest frame | the same, stated as the definition of "Minkowski theories"; cutoff Λ_UV on space momenta | match |
| D1, the mixed structure | "a propagator calculated on the pole of another falls off half less rapidly" (§5, eq. 5.1–5.2) | match |
| D2, ln Λ/(4π² s) for the bare pair | not stated; their D = 4 trivial-numerator entry is c₀,₀ = 0 for S × S, listed in (3.6) without explanation (c₁,₀ = 0 gets one, "an eikonal cancellation") | reconciled by P7 at m = 0; P7 is the explanation of c₀,₀ = 0 |
| D3 (wrong): frame-dependent divergence | they take p = (0, p_s) and "find (3.2) again": the divergent part is a function of p² only | they had it right; our corrected result agrees |
| P3: the finite part is frame dependent | **not in the paper** [v, full text of §3 read 3 Nov 19:10]: "Σ(p) is expected to be a function of p² only"; computed at p = (p₀, 0) with "the Lorentz invariant result retrieved by means of the replacement p₀² → p²"; spacelike p "obtained by means of the analytic continuation in p²"; the explicit p = (0, p_s) check recovers (3.2), whose finite part is written "+ (finite)"; footnote 1: a Lorentz-invariant regularization would be more rigorous, "we are only interested in the logarithmic divergences, which … are independent of this choice" | **the replacement p₀² → p² is justified for the log and not for the finite part**, and the paper applies it to both. What the era did not write |
| D4: (++) sign flip; local log cancels in the four-pair sum | not in that form; "non-Hermitian" counterterms | same fact, other side |
| no Euclidean region to subtract at | §5: "the analytic continuation to Euclidean space is not possible" | match |
| D = 6 headline Σ = −M⁴/(2(4π)³)[M²/(p²)² − i/p²] ln(Λ²/M²) + … | | same mechanism one dimension up, two nonlocal structures |

Reading note (author, 3 Nov 17:50): their large-k_s "pinch", Ω ≈ p_s − iη and Ω̄ ≈ p_s + iη with η = M²/(2p_s) closing on the real k⁰ axis
from opposite sides, is a conjugate pair squeezing the contour, not a pole collision — the same first misreading hdym's Krein pairs invited.

## 6. What the control says for the ladder

The three prescriptions of the era are now side by side on one integrand: the Lee–Wick contour (Lorentz-invariance crack at p⃗ ≠ 0, finite,
part 1), CLOP (no cure of that crack, part 1), and no prescription (a covariant nonlocal divergence and a frame-dependent finite part, this
note). The 2017 region-wise continuation and the 2018 average continuation are the answers; they come on the back sweep from 2017 (ordering
rule). Rung 2's forward pass is complete for the era; review against intent with the author before rung 3 ('t Hooft–Veltman 1974).

Status: `minkowski.log` (3c0ba46) status=PASS (P1, P2, P4, P5, P7). Readings accepted by the author 3 Nov 19:23 (review against intent):
"CLOP's Lorentz invariance is method P's" (part 1); "p₀² → p² holds for the log, not the finite part"; "the Lee–Wick contour and the
real axis swap which divergence is local". Fresh-agent review (rule 9)
of 3 Nov evening: 18 findings, 2 must-fix (the `B_rest` windows applied to the mixed term — fixed, numbers above updated; the claim that the
author's log and the smoke agreed on every number — withdrawn, they differed at 10⁻⁶ in exactly the affected block), all folded in; the
reviewer re-derived the residue sums, the tails and the angular identity and found them correct. Claim tags (hdym's, adopted 3 Nov pending the author's yes): the
residue formulas and the tail [exact] (sympy); every number above [grid] on the committed Λ and s; the angular-average identity [exact].
