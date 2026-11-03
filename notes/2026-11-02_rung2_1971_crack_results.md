# Rung 2 results, part 1 (2 Nov 2026): the 1971 Lorentz-invariance crack, two techniques of the era

Setup and the agreement criterion: `notes/2026-11-02_rung2_1971_crack_summaries.md` (§1, §4). Code and logs: `bubble/` (concept folder;
`bubble_frame.py`, `bubble_clop.py`, machinery appended to `bubble_lib.py` in a dated block). Ladder: chronology note §5, rung 2.
Ordering rule used (author, 2 Nov evening; written into the chronology §5 as a dated correction in the same patch): regulators and computational
methods may be applied at any point of the forward pass; a *prescription* that defines the theory where the integral is undefined enters
only when the timeline reaches it. So this part is 1969–71 only: the Lee–Wick contour with real spatial momenta (Nakanishi's
computation) and CLOP; the 2017 region-wise continuation, the 2018 average continuation and Anselmi–Piva's model wait for the back
sweep from 2017. Method P (Feynman parameters, average continuation) appears here only as the Lorentz-invariant anchor that the rest
frame was already shown to equal on 30 Oct.

All numbers below are **sandbox smoke** (2 cores, dps 15); the author's full runs are pending (rule 15). Model: the complex-conjugate
pair a = ±i (M = 1, the Lee–Wick pair of 1801.00915 eq. 5.5, so p⃗ = 0 is the committed `results/lw.json`); D = 4; UV subtraction at the
Euclidean point s = p² = −1 in the same frame.

## 1. Technique of the day: the Lee–Wick contour with real q⃗ at p⃗ ≠ 0 (`bubble_frame.py`)

**What is computed.** The k⁰ integral along the Lee–Wick contour by residues, as `bubble_lib.B_momentum` does in the rest frame, with
ω₂ = √((k⃗ − p⃗)² + a₂); the remaining 3-dimensional integral over real k⃗ in cylindrical coordinates; compared at equal s = p₀² − p⃗² with
the rest-frame principal-value (Lee–Wick) amplitude and with method P.

**Predictions (header, before the first run) and outcomes:**

| | outcome |
|---|---|
| P1: below s_min(|p⃗|) = 4(Re E_{p⃗/2})² − p⃗² the boosted and rest-frame values agree | **holds**, to 10⁻¹⁷ at s = 1.0, 1.5 for |p⃗| = 0.3, 0.8 |
| P2: above s_min they differ; at s = 1.8 only |p⃗| = 0.8 (s_min = 1.705) has crossed, |p⃗| = 0.3 (s_min = 1.956) has not | **holds**, exactly that pattern; differences 1.5·10⁻³ to 10⁻² |
| P3: continuous across s_min (no jump) | **holds** on the corrected numbers: two-sided difference 4.7·10⁻⁴, 1.6·10⁻⁴, 5.5·10⁻⁵ at h = 0.1, 0.03, 0.01 — roughly linear in h |
| P4: |p⃗| → 0 at s = 2.5 tends to the rest-frame **principal value**, not to either iε value | **holds**: |B − PV| = 2.97·10⁻³, 1.01·10⁻³, 3.0·10⁻⁴ at |p⃗| = 0.3, 0.1, 0.03; |B − (±iε)| stays 1.2·10⁻² |
| P5: method P = rest-frame method M (PV) at every s | **holds** to 10⁻¹² |

**The crack, as numbers** (sweep, Im B; Re B ≈ 10⁻²⁰ everywhere):

| s | rest frame (PV) | |p⃗| = 0.3 | |p⃗| = 0.8 |
|---|---|---|---|
| 1.0 | 0.005330679809 | same to 10⁻¹⁷ | same to 10⁻¹⁷ |
| 1.5 | 0.008162407909 | same to 10⁻¹⁷ | same to 10⁻¹⁷ |
| 2.0 | 0.01835553842 | 0.01422716657 | 0.008560223554 |
| 2.5 | 0.01373245291 | 0.01076563082 | 0.00654714647 |
| 3.0 | 0.01049714443 | 0.008233642335 | — |
| 4.0 | 0.006159520298 | 0.004683950527 | — |

(|p⃗| = 0.8 at s = 3.0, 3.5, 4.0 in `frame_sweep.log`.) Readings, not predicted: the frame dependence sits entirely in the imaginary part
(the would-be discontinuity channel) and is linear in |p⃗| near p⃗ = 0 (the P4 numbers), not quadratic.

**A correction, recorded (rule 11).** The header said there is no pinch at p⃗ ≠ 0 and that the denominator vanishes at one point of R³.
Wrong: the symmetry k⃗ → p⃗ − k⃗ fixes the plane k_∥ = p/2 pointwise, so ω₁ + ω₂ is real on that whole plane and vanishes-with-p₀ on the
**circle** k_∥ = p/2, k_⊥ = k*(p₀) for s > s_min — the rest-frame pinch sphere becomes a circle. Im(ω₁ + ω₂) is odd across the plane,
so the 3-dimensional integral is absolutely convergent (a 1/r singularity on a circle) and well defined without a prescription — the
iterated k_⊥ integral picks up ±iπ times a residue on the two sides of the plane, and their average is the principal value — but a
quadrature that does not resolve the near-pole at k* gives breakpoint-dependent numbers: the first run's above-threshold values were off
by 0.15% to 7% (worst at |p⃗| = 0.3, s = 4.0: 0.0043486 against 0.0046840; at s = 2.5, |p⃗| = 0.8: 0.0065774, and 0.0065237 in an
unlogged exploratory run with other breakpoints, against the converged 0.0065471464695). Found by the GMRS certification below (agreement
only to 1.5% in an unlogged exploratory run), not by any check in the script. Fixed by the pinch-aware quadrature
`bubble_lib.B_momentum_boost` (Gauss–Legendre windows of width ∝ |Im(ω₁ + ω₂)| around the near-pole; converged to 14 digits across window
choices, `results/superseded/frame_window_variation_check.log`); the first logs kept as `results/superseded/{frame,clop}_sandbox_pointwise_quadrature.log`. No prediction was changed; P1, P2,
P4, P5 were unchanged in kind, P3 re-scored. And P4 got its reason: the average of the two iε signs across the plane *is* the principal
value.

**Certification against the GMRS closed form (Nakanishi 2006 eq. 3.7), `--gmrs`.** Two steps, both derived rather than assumed.
(a) The "fish" D is the image of the (k_∥, k_⊥ ≥ 0) half-plane under f = E_k + E*_{p−k}; its boundary lies in the image of the edge
k_⊥ = 0 or of the fold set (Jacobian of (Re f, Im f) zero). The edge is a fold by evenness in k_⊥; the Jacobian keeps one sign on a
13×12 grid over k_⊥ > 0 and at |k| ≈ 50 (min |J| on the grid 10⁻³, 2·10⁻³ for |p⃗| = 0.3, 0.8), so Γ lies in the image of the axis and
touches the real p₀ axis at k = p⃗/2 — Nakanishi's touching point, checked on a grid (a grid cannot exclude a sign change between nodes;
not a proof). (b) With Γ_s the Im s′ < 0 half of the image of Γ under s′ = f² − p², the closed form
F(s) = (1/16π²) Re ∫_{Γ_s} ds′/(s − s′)[ |p⃗|(M² − M*²)/(s′√(s′ + p⃗²)) + √((s′ − (M+M*)²)(s′ − (M−M*)²))/s′ ] with principal square roots
satisfies **B(s) − B(s₀) = −i [F(s) − F(s₀)] with the constant equal to −i to 10⁻¹⁴ over the whole sweep**, both frames — and −i is
exactly the factor the two normalisations predict (F = iB from Nakanishi's 1/(2(2π)³) against our i/(4π²); derived after the first
sweep, now scored). Taking F without the Re, or the other half of Γ, or the conjugate branch, fails by O(1) (probed in an unlogged
exploratory run and discarded; the choice is post hoc and disclosed). So the 1971 computation is reproduced, with the independent closed
form of GMRS as the second implementation — the independent second implementation the author asked for as the real check of a rung
("steam vac"), here for the rung's first number.

## 2. The next technique: CLOP (`bubble_clop.py`)

M → M + δ in one propagator, compute, δ → 0 from both sides (1703.04584 §4's description of Cutkosky–Landshoff–Olive–Polkinghorne
1969; the paper itself unseen). δ ∈ {±0.1, ±0.03, ±0.01}.

| | outcome |
|---|---|
| P1: rest frame, method M — the two limits are the two iε values, not the PV (Im(ω₁+ω₂) ∝ δ of definite sign lifts the pinch to one side) | **holds**: δ → 0⁺ → −iε, δ → 0⁻ → +iε (shift on line 1), monotone; the gap between the limits is 0.989 and 0.993 of the full iε gap at s = 2.5, 4.0 |
| P2: rest frame, method P — upper = lower continuation for δ ≠ 0, methods P = M, limits are f₊, f₋ | **holds**: upper − lower = 2·10⁻¹¹, P − M = 10⁻¹¹ at every δ; the approach to f₋ (δ > 0) / f₊ (δ < 0) monotone |
| P3: the full Lee–Wick propagator (four-pair sum, shift on line 1's pair) — the cross terms' iπ pieces cancel, no ambiguity for equal masses | **holds**: both limits approach the PV sum monotonically |
| P4: boosted frame |p⃗| = 0.8 — the shift is a regular perturbation and δ → 0± returns the frame-dependent integral | **holds** (re-run with the near-pole windowed for the shifted, non-conjugate pair as well — the second review found the first re-run's δ ≠ 0 points still pointwise and moving by 1–3% between breakpoint sets): |B(δ) − B(0)| decreasing from both sides, numbers in `clop.log` |

Reading (P5 of the header, not scored): method P is Lorentz invariant by construction, so "CLOP + method P" is frame independent while
"CLOP + method M" is not. The Lee–Wick reply of 1971 restores Lorentz invariance by replacing the real-momentum (Hamiltonian) integral
with the covariant parametric one plus a rule for its singularities — which is Nakanishi's objection in one sentence (the definition
changed, not the integral; hep-th/0609206 §1). And on a single (M, M*) pair CLOP's ambiguity is the *entire* would-be discontinuity; it
cancels only in the symmetric four-pair sum — consistent with 1703.04584 §4's ±π for unequal masses only, on the 2 Nov reading of that
section; a 2017-model statement, left for the back sweep.

## 3. Where rung 2 stands on the forward pass

Smoke passed for (i) the technique of the day and (ii) CLOP, each with predictions first, P and M both, the closed form as certification;
done when the author's logs are committed (rule 15). Open on this
pass: (v) the Minkowski real-axis k⁰ integral as the era's control (rule 13; expected: the nonlocal divergence of 1612.06510) — next.
Waiting for the back sweep from 2017: (iii) region-wise continuation, (iv) average continuation on the 1971 model, Anselmi–Piva's own
model against their fig. 14, the unequal-mass CLOP ±π.

Ledger entries for `notes/TECHNIQUES.md`: the pinch-circle geometry (the boost turns the pinch sphere into a circle whose iε sign flips
across the fixed plane of the k → p − k symmetry; the average of the signs is the principal value) — reusable for any complex-conjugate
pair at p⃗ ≠ 0; the fold-set argument for the boundary of an image domain; the "conditionally convergent, so resolve the near-pole with
windows scaled to the distance from the symmetry plane" quadrature.

History: written 2 Nov evening from the sandbox logs; to be re-read against the author's logs.
