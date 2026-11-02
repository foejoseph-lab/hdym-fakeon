# Rung 1 results: Utiyama & DeWitt 1962, redone three ways (31 Oct 2026)

Setup and scoring: `notes/2026-10-31_rung1_utiyama_dewitt_setup.md`. Ladder: `notes/2026-10-31_chronology_1939_2014.md` §5.
Code and logs: `matter_loop/` (concept folder — one-loop matter on a classical background; named by
concept, not year or author, per the author 31 Oct; techniques accumulate inside it).

## 1. What was computed

One real scalar (mass m, coupling ξ) on g = δ + h, Euclidean, flat-space perturbation theory to second
order in h, one loop: the bubble (two first-order vertices) and the tadpole (one second-order vertex).
The one-loop effective action Γ₂ = ∫_p ½ h(p)·Π(p)·h(−p), reported as the quadratic form Q = Γ₂/(2·Vol)
whose value equals the coefficient of ∫d⁴x F for each covariant structure F. All 256 tensor components kept
and fitted; nothing projected by hand. Vertices derived from the action (`ud1962_lib.py`, header); the
second-order densities √g⁽²⁾, (√g g^μν)⁽²⁾, (√g R)⁽²⁾ computed from the metric (Christoffels → Ricci) rather
than copied.

**Technique sequence, honestly labelled** (corrected 31 Oct after independent review): (a) 1962 hard momentum
cutoff — done, for one cutoff placement only (on the Feynman-parametrised, shifted loop momentum; the
unshifted placement promised in the setup note was not run); (b1) Schwinger 1951 proper time applied to the
flat-space perturbative graphs, lower cutoff t ≥ 1/Λ² — done; **(b2) DeWitt 1965, the covariant recursion for
the heat-kernel coefficients on an arbitrary background, no expansion in h — NOT done**; its result was taken
from Vassilevich's table and used as the check. On the ladder's rule (b2) is a rung technique in its own
right and is the first item of the next session; it is also the only technique here that sees the Riem²
coefficient of a₄, which Gauss–Bonnet hides at second order in h. (c) 1972–74 dimensional regularization —
done. Signature: Euclidean throughout, including (a); the setup note said Lorentzian for 1962 — a deviation
made for uniformity, which does not affect divergence coefficients.

Conventions needed to read the comparison: our R is +2 on the unit S²; **'t Hooft–Veltman's R has the
opposite sign** (their eq. 3.18, R = −3F⁻²f_γγ − …, is negative on spheres), and their Lagrangian is
L = −½(∂φ)² + … with Euclidean metric, i.e. L_theirs = −L_E. Both flips are used below; with them the
"M → −ξR" step gives ¼(ξ − 1/6)², without them it would give ¼(ξ + 1/6)² and a spurious mismatch.

Three regulators, two choices of path-integral measure:

| | (a) 1962: hard momentum cutoff Λ on the Feynman-parametrised loop momentum | (b1) 1951 Schwinger proper time on the perturbative graphs, t ≥ 1/Λ² | (c) 1972–74: dimensional regularization, 1/ε = 1/(8π²(n−4)) |
|---|---|---|---|
| measure φ (operator K = √g·D) | `rung1_cutoff_plain.log` | `rung1_propertime_plain.log` | `rung1_dimreg_plain.log` |
| measure χ = g^{1/4}φ (operator g^{−1/4}Kg^{−1/4}; DeWitt's covariant measure) | `rung1_cutoff_covariant.log` | `rung1_propertime_covariant.log` | `rung1_dimreg_covariant.log` |

Run times 1–8 min each, single-threaded sympy; exact rationals throughout. The superseded scripts and logs
that show the development (including a wrong prediction, §4.1) are kept in `results/superseded/`.

## 2. Results

Units: coefficient of (1/16π²) log Λ² for the log order; of (1/16π²)Λ⁴, (1/16π²)Λ² for the power orders;
of 1/ε for (c). Structures: √g⁽²⁾; (√g R)⁽²⁾; Ric² ≡ (R⁽¹⁾_μν)²; R² ≡ (R⁽¹⁾)².

**Log order — identical in all six runs:**

    p⁰ :  −¼ m⁴ √g                      (one-loop vacuum energy, from ½ Tr ln)
    p² :  −(6ξ−1)/12 · m² √g R  =  ½(1/6 − ξ) m² R      (vanishes at conformal coupling)
    p⁴ :  −(1/120) Ric²  −  [¼(ξ − 1/6)² − 1/360] R²

In (c) the same numbers appear with the opposite sign as coefficients of 1/ε, which with ΔL_tHV = +Γ_div
(the two sign flips above) are exactly **'t Hooft–Veltman 1974 eq.
(3.40)** for a real scalar: (1/120) R_μνR^μν + ¼(M − N² − D·N − R/6)² with M → −ξR, i.e. R² coefficient
¼(ξ−1/6)² − 1/360, Ric² coefficient 1/120. At ξ = 1/6 the square vanishes and the p⁴ part is
(1/120)(Ric² − R²/3) = (1/240) C² mod Gauss–Bonnet. The full ξ-polynomial matches, not just the endpoints.

**Power orders:**

| | Λ⁴ | Λ² |
|---|---|---|
| (a) cutoff, φ | −(1/192)(h² + 2 h_ijh_ij): **not** ∝ √g⁽²⁾ = ⅛(h² − 2h_ijh_ij) | not a combination of m²√g⁽²⁾ and (√g R)⁽²⁾ |
| (a) cutoff, χ | still not covariant | still not covariant |
| (b1) proper time, φ | −(1/32)(h² + 2h_ijh_ij): not covariant; off from −¼√g⁽²⁾ by exactly −(1/8) h_ijh_ij = ½·(ln√g)⁽²⁾ | not covariant |
| (b1) proper time, χ | **−¼ √g⁽²⁾** | **½ m² √g⁽²⁾ − ½(1/6 − ξ)(√g R)⁽²⁾** |
| (c) dim reg, either | absent | absent |

**Heat-kernel coefficients at second order in h, from the small-t expansion of Tr e^{−tD}** (b1, χ measure) —
derived perturbatively from the same Feynman rules, *not* by the DeWitt recursion (that is (b2), not done):

    a₀⁽²⁾ = √g⁽²⁾
    a₂⁽²⁾ = [√g(−m² + (1/6 − ξ)R)]⁽²⁾
    a₄⁽²⁾ = [√g(½m⁴ + (ξ − 1/6)m²R + (1/60)Ric² + ((60ξ² − 20ξ + 1)/120)R²)]⁽²⁾

each equal, term by term, to Vassilevich hep-th/0306138 eqs. 4.26–4.28 with E = −ξR − m², Ω = 0, expanded
to second order (Riem² → 4Ric² − R² at this order by Gauss–Bonnet; total derivatives drop). This checks the
O(h²) projection of a₄, not a₄ itself: the Riem² coefficient is not separately visible here. With the φ
measure a₄⁽²⁾ still agrees but a₀⁽²⁾ is off by ½ h_ijh_ij and a₂⁽²⁾ is not covariant.

## 3. Scoring against the 1962 abstract and the checks

| claim (abstract) | reproduced? | how |
|---|---|---|
| three orders of divergence ∞⁴, ∞², log∞ | yes | (a) and (b) both show exactly three; (c) shows the log only, as it must |
| log∞ removed by a counterterm quadratic in the Riemann tensor | **yes, with numbers** | p⁴ log coefficients = 't Hooft–Veltman (3.40), six ways |
| ∞⁴ removed by a cosmological term | **yes, with proper time and the covariant measure** | (b1)+χ: −¼√g; (a), shifted-momentum placement, either measure: not covariant; (b1)+φ off by the measure term |
| ∞² removed by a renormalization of G | **same condition** | (b1)+χ: ½m²√g − ½(1/6−ξ)√gR, the a₂ structure |
| semiclassical graviton Green's function ∼ 1/p⁴ | **argued, not computed** | the p⁴ log(Λ²/p²) term of Π, resummed: G⁻¹ ≈ M_P² p² P + (c/16π²) p⁴ log(p²/Λ²) with c = 1/120 on the spin-2 projector, so G ∼ 1/(p⁴ log p²) for p² ≫ M_P²; this is Stelle's 1/k⁴ fifteen years early, with a running coefficient rather than a constant — the bridge to rung 3 and to 1977. The finite p⁴ log p² part is in the code (`'fin'`) but was neither printed nor fitted; to be computed with (b2) |

Checks, earliest first (setup note §1): 't Hooft–Veltman 1974 eq. (3.40) — matched exactly, in their
normalization; Salvio 2018 eqs. 3.15–3.16 — the per-scalar N_S/60 and (1/12)(1+6ξ)² follow from the p⁴
coefficients by d(1/(2f₂²))/d ln μ ⟷ (1/16π²)(1/120)·2 and d(1/(6f₀²)) ⟷ (1/16π²)·2·[¼(ξ−1/6)² − 1/360] up to
the Gauss–Bonnet rearrangement into (W², R²) and Salvio's sign convention for ξ (his (1+6ξ) ⟷ our (1−6ξ),
i.e. ξ_Salvio = −ξ_here); the N_S/60 is immediate, the f₀ term is the (1+6ξ)² piece once the −1/360 is moved
into the W² coefficient. (Written out, not run: a one-line conversion.)

## 4. Findings beyond the reproduction

1. **The measure matters at ∞⁴ and ∞², and only there** [corrected 1 Nov: and in the *finite* p⁴ constant; the log
   order never sees it — §8.2]. With φ as the integration variable the Λ⁴ term is
   off from −¼√g⁽²⁾ by −⅛h_ijh_ij = ½(ln√g)⁽²⁾ (units Λ⁴/16π²); with χ = g^{1/4}φ it is right. The
   identification with ln√g is a hand step on the logged difference, not a separate logged check. **How this
   was found, for the record:** the first proper-time script (`results/superseded/ud1965_propertime.py`)
   predicted in its header that "proper-time regularization respects general covariance, so all three orders
   should now assemble"; its log says `all three orders covariant and as expected: False`. That prediction was
   wrong as written — the regulator was covariant, the operator was not — and the covariant measure was built
   in response. The status label `LOG-OK_POWER-ORDERS-NONCOVARIANT` for that case was introduced afterwards,
   once the reason was understood; it describes the result, not a prediction. The log order and
   't Hooft–Veltman's result never see this — dimensional regularization sets the offending δ(0) to zero
   (their remark below the tadpole figure, "∫dⁿp p_μp_ν/p² = 0 in the continuous dimension method", is this
   statement). The covariant measure is Fradkin & Vilkovisky, Phys. Rev. D 8, 4241 (1973) [m?: reference from
   memory, to verify] / DeWitt 1967; it was found here by the computation disagreeing, not quoted.
2. **The Feynman-parametrised hard cutoff does not produce the 1962 assignments for ∞⁴ and ∞²**, with either
   measure. Only this placement of the cutoff was run; the unshifted placement (setup §6) was not, so this is
   one regulator variant failing, not a theorem about cutoffs. What a covariant regulator gives is (b1)+χ. Whether
   Utiyama–DeWitt had one (DeWitt's proper-time methods existed by 1962; the 1965 book systematizes them) or
   argued past non-covariant pieces is not decidable from the abstract. What is decidable: the 1965 move to
   proper time was necessary for the ∞⁴ and ∞² halves of the 1962 claim to be well-defined.
3. **(√g g^μν)⁽²⁾ is traceless in four dimensions** (known: det of the densitised metric is 1 in n = 4), so
   the seagull's derivative part vanishes identically and the tadpole is pure m² and ξR. Not a finding — a
   textbook fact met in the code; recorded because it is why the Λ⁴ term is bubble-only with the φ measure.
4. **The dead horse, partly flogged:** the three log coefficients came out identical under cutoff, perturbative
   proper time and dimensional regularization, with either measure — six runs. No crack at the log order in
   the techniques run; the 1965 covariant recursion (b2) is still owed.

## 5. Lessons for later rungs

- Keep all tensor components and fit; projecting by hand would have hidden the non-covariant Λ⁴ structure.
- Derive the second-order densities from the metric every time; the tracelessness of (√g g^μν)⁽²⁾ and the
  sign of (√g R)⁽²⁾ are the kind of thing one gets wrong by hand.
- Bubble-only numbers are not checks: the vacuum-energy "check" of the first run (+¼ m⁴) was wrong by the
  missing tadpole (−½); only the sum (−¼) means anything. Recorded so it is not repeated at rung 3, where
  the graviton loop has many more pieces.
- The covariant-measure vertex is the K-dependent correction −(δ_μν/4)[K₀(k)+K₀(q)]; at rung 3 the analogous
  measure/gauge-fixing bookkeeping is where 't Hooft–Veltman spend their §4–5.
- Runtime is dominated by expanding 256 components of products of 4×4 symbolic matrices; it is embarrassingly
  parallel over components. Rung 3 gets a `--jobs` switch and runs on the author's 22 cores.

## 6. Open

- **(b2) DeWitt 1965 covariant recursion** for a₀, a₂, a₄ on an arbitrary background — the technique skipped;
  first item next session. Its O(h²) expansion must reproduce §2; its Riem² coefficient is new information.
- The 1/p⁴ resummation: print and fit the finite p⁴ log p² part of Π, then resum — "argued" above, not done.
- Dependence of the (a) power divergences on cutoff placement (shifted vs unshifted loop momentum): not run;
  needed before §4.2 can say anything stronger than it now does.
- Spinor and vector matter (the N_F/20, N_V/5 of Salvio 3.15): not done; additions to the same code if the
  ladder wants them before rung 3 (the vector would be abelian — Maxwell — not Yang–Mills).
- The 1962 text: unseen; §3 re-scored in one line if it ever is.

## 7. Addendum, 31 Oct evening: (b2) DeWitt 1965 covariant recursion — done

Setup, predictions and the two corrections made before the logged runs: `2026-10-31_rung1_b2_dewitt_recursion_setup.md`.
Code `matter_loop/rung1_dewitt.py`; logs `results/rung1_dewitt_samples.log` (the result) and
`results/rung1_dewitt_quartic0.log` (the control). Sympy, exact rationals, 6 random rational algebraic
curvature tensors in parallel; 6.5 min on 2 cores in the sandbox, `--jobs` over the samples. **The two logs
were produced by the code as it stood before the review fixes of 31 Oct night** (the only later change to
the code is two extra report fields per background: the `dRiemann_at_0` check and the count of nonzero x²
coefficients of R(x); the computation is untouched). Under the smoke-test rule (smoke test in the sandbox, full runs on the author's machine) they are to be
regenerated on the author's machine with the two README commands and committed in their place.

**Result.** The recursion σ^{;μ}∇_μa_k + k a_k = Δ^{−1/2}(□+E)(Δ^{1/2}a_{k−1}), run in normal coordinates on
generic backgrounds (no expansion in h), gives at coincidence

    [a₁] = E + R/6
    [a₂] = ½E² + ⅙RE + (1/72)R² − (1/180)Ric² + (1/180)Riem²        (□R = 0 on every background used)

as an exact fit over 6 backgrounds, 18 equations, 5 unknowns, residual zero; unit S⁴ gives 29/15 and S²×S² gives
11/45 at E = 0, as predicted. Both constraints inherited from the six O(h²) logs hold exactly
(c_{Ric²} + 4c_{Riem²} = 1/60, c_{R²} − c_{Riem²} = 1/120), and the numbers equal Vassilevich 4.28, read after the
run. In the (W², R², Euler) basis: W²/120 + R²/72 − E₄/360 at E = 0, E₄ = Riem² − 4Ric² + R²; with E = −ξR the
R² term is ½(ξ − 1/6)²R². **The Euler coefficient −1/360 is the number the flat-space graphs could not see;
the log order of rung 1 now has four techniques agreeing, one of them not a graph.**

**Control (`--quartic 0`).** With the x⁴ term of the normal-coordinate metric removed, □R(0) ≠ 0 on every
background (S⁴: 32; S²×S²: 64/9), [a₂] shifts (S⁴: 3, S²×S²: 13/27) by exactly (1/30)□R(0), and the fit over
generic samples degenerates to a one-parameter family, c₂ = 1/108 − 4c₄/9, c₃ = 1/36 − 2c₄/3, which contains
the result at c₄ = 1/30. Within this family of metrics □R(0) is itself a combination of R², Ric², Riem², so
the curvature-squared coefficients are unambiguous only on a metric with □R(0) = 0 — which is what the
standard (2/45) R·R x⁴ term delivers. Two things were learned from this control and are in the setup note's
History: the x⁴ term is load-bearing (the code header had claimed otherwise, and was corrected), and a
"check" first planned — ∇∇R symmetrised in the derivative pair alone, vanishing at the origin — cannot hold
for generic curvature and was replaced before the logged runs (the fully symmetrised object does vanish on
this metric; see the setup note). Reviewed by a fresh agent 31 Oct night before the patch was sent; its
findings (text claims exceeding the code and the logs, none touching the numbers) are folded in and listed
in the setup note's History.

**Scoring.** Agreement on all counts; no crack at the log order in four techniques. DeWitt 1965 itself stays
unseen (rule 10); the recursion is attributed to it through Avramidi hep-th/9510140 and Vassilevich hep-th/0306138.

**Still open from §6:** the 1/p⁴ resummation (finite p⁴ log p² part of Π), the unshifted-cutoff placement for
§4.2, spinor and vector matter.

## 8. Addendum, 1 Nov evening: the three owed items — rung 1 closed

Code: `matter_loop/rung1_run.py --placement unshifted` and `matter_loop/rung1_finite.py` (new), machinery appended to
`ud1962_lib.py` (dated block). Smoke-tested in the sandbox at full settings (there is no smaller setting: 256 components
either run or do not); the logs below are to be regenerated on the author's machine with the README commands and
committed (rule 15). Predictions were written in the two headers before the runs; one was corrected before the run and
the correction is in the header (P5, below).

### 8.1 Cutoff placement (`rung1_run.py --regulator cutoff --placement unshifted`)

Hard cutoff |k| < Λ on the propagator momentum itself, no Feynman parameter; the two-propagator product expanded at
large k (every power of p lowers the UV degree by one, so the Λ⁴, Λ², log Λ² coefficients are exact and the finite part
from this expansion is not reported). The script also recomputes the shifted placement and prints the difference.

Predicted (surface-term formula for shifting a cut-off integral, header): Λ⁴ identical, log identical, Λ² different by a
non-covariant p² term. **Found, plain measure:** Λ⁴ identical; log identical; Λ² (unshifted − shifted) =
(1/576) p² (h² + 2 h_ij h_ij) in units Λ²/16π² — the same non-covariant tensor structure as the Λ⁴ term itself,
−(1/192)(h² + 2 h_ij h_ij). Status `LOG-OK_POWER-ORDERS-NONCOVARIANT`, as predicted.

**§4.2 upgraded:** the hard momentum cutoff fails covariance at Λ⁴ with *either* placement, Feynman-shifted or unshifted
(the Λ⁴ term does not depend on the placement at all); the placement dependence is confined to Λ². The coefficient
(1/576) and the tensor structure were not predicted; only which orders move and which do not. What 1962 did with its ∞⁴ and ∞²
remains undecidable from the abstract; what is now decided is that no placement of a hard momentum cutoff gives them.
Covariant-measure variant, run on the author's machine 1 Nov night (`rung1_cutoff_covariant_unshifted.log`): the same (1/576)p²(h² + 2h_ijh_ij) shift at Λ², Λ⁴ and log identical — the placement shift is the same with either measure.

### 8.2 The finite p⁴ part (`rung1_finite.py`)

Massless scalar; the p⁴ part is bubble-only because the tadpole carries no p, but it is the full sum at this order.
Q(p⁴) = c_L log Λ² + c_p log p² + K, each fitted to {Ric², R²} and to the Barnes–Rivers forms p⁴ h·P⁽²⁾·h, p⁴ h·P⁽⁰ˢ⁾·h
(Ric² = ¼P⁽²⁾ + P⁽⁰ˢ⁾, R² = 3P⁽⁰ˢ⁾, checked in the log). Three regulators: (a) Feynman-parametrised hard cutoff;
(b1) proper time on the graphs; (c) dimensional regularization, MS-bar, **with the n-dimensional angular averages and
n-dimensional dot products** (the 31 Oct dim-reg driver took them in d = 4, which is exact for the pole and wrong for
the finite part — a bug of 1 Nov's first draft, caught by c_L coming out 16π² off and by the {Ric², R²} fit of K failing in dim reg where
it cannot; the explicit Ward test was added together with the fix and then found the measure effect; fixed before any
number was recorded). Two partial sandbox runs preceded the final header and are kept, with their timestamps, in
`results/superseded/rung1_finite_{cutoff_only_sandbox_2027,dimreg_only_sandbox_2041}.log` (the full sandbox smoke logs
were superseded by the author's runs the same night and removed in the housekeeping patch); the first three-regulator
run, with the d = 4 averages, was overwritten and survives only in the session transcript.

| prediction (header) | outcome |
|---|---|
| P1: c_p = −c_L in all three (bookkeeping; cannot fail if the split is right) | **holds**, all three, both measures |
| P3: spin-2 log coefficient −1/480 for all ξ; spin-0 −(3/4)(ξ−1/6)², zero at conformal coupling (arithmetic on the 31 Oct c_L and the projector identity) | **holds** |
| P4: tree −(M_P²/2)∫√g R → +M_P²/8 on P⁽²⁾, −M_P²/4 on P⁽⁰ˢ⁾ (Fierz–Pauli 1 : −2, Euclidean conformal-factor sign) | **holds**, derived from (√g R)⁽²⁾ in code |
| P2: K regulator-dependent, differences local and covariant | **holds between the covariant regulators; fails for the cutoff and for the plain measure — see below** |
| P5: see 8.3 | corrected before the first execution; the log tests only a > 0, b > 0 (P4 and P3 again) |

**K, the regulator constant** (units 1/16π²; dim reg: pole dropped, MS-bar, μ² identified with Λ²; a quadratic form in h
is "transverse" when its first-order variation under h → h + pε + εp vanishes identically in h and ε — the linearized
Ward identity as a bilinear statement, tested in the log; the weaker test K(pure gauge) = 0 is not equivalent, as the
1 Nov review showed with (p h p)² − p⁴h², and the first version of the script used the weaker one; the verdicts did not
change because every "yes" row is also an exact {Ric², R²} fit):

| regulator | measure | K transverse? | K |
|---|---|---|---|
| (c) dim reg MS-bar | either (δ(0) = 0) | **yes** | −(23/900) Ric² − (1800ξ² − 500ξ + 1)/3600 R² = (1/120)(−46/15) Ric² + … |
| (b1) proper time on graphs | χ = g^{1/4}φ | **yes** | K_dimreg + γ_E·(1/120 Ric² + (60ξ²−20ξ+1)/240 R²) — i.e. **exactly −γ_E·c_L**, the known proper-time-cutoff ↔ MS-bar offset |
| (b1) proper time on graphs | φ | **no** | non-transverse part −(5ξ−1)/240 (p h p)² + (20ξ−3)/960 p⁴h² (on the complete 5-element basis, log) |
| (a) hard cutoff (shifted placement; the unshifted finite part was not computed) | either | **no** | non-transverse part, φ: −(15ξ−2)/1440 (p h p)² + (1/720) p²(p h h p) + (6ξ−1)/576 p⁴h²; χ: different numbers, log |

The status label `LOG-OK_FINITE-NONTRANSVERSE[...]` and the rule that P2 pairs containing a non-transverse K are reported
but not scored were introduced after the first, cutoff-only run (which ended `status=OK` with the same non-transverse
K); they describe the result, not a prediction — as in §4.1. **P2 as written failed** for the cutoff and for proper time
with φ. Also from the logs: the measure changes K only in its spin-0 and non-transverse parts — the spin-2 (TT)
constant is measure-independent in every regulator (cutoff −37/14400 with φ and with χ; proper time (−46+15γ_E)/7200
with both; the fit's P⁽²⁾ coefficient equals h·P⁽²⁾KP⁽²⁾·h because the three non-transverse basis elements vanish under
P⁽²⁾ on both sides).

Readings. (i) The −46/15 on Ric² in MS-bar (μ² ≡ Λ² in the table) is, from memory, the massless-scalar
graviton-self-energy constant of the 1974 dim-reg literature [m?: the scalar-loop paper is Capper & Ramón Medrano PRD 9,
1641 (1974), not Capper–Duff–Halpern PRD 10, 461, which is the photon loop; the Barvinsky–Vilkovisky second-order form
factor for the massless scalar carries 23/450 against a log coefficient 1/60, ratio 46/15 — both unverified; to be
checked against a listing before it is leaned on]. It was not aimed at. If a listing confirms it, it is an independent
check of the n-dimensional machinery; until then the checks are the bilinear Ward identity and the γ_E relation to
proper time. (ii) **Finding 1 of §4 is corrected:** the measure matters at ∞⁴, ∞² *and in the finite p⁴ constant*; the log order
is the only order that never sees it. With φ the proper-time finite part violates the linearized Ward identity by a
local, ξ-dependent polynomial; with χ it is transverse and differs from MS-bar by γ_E alone. (iii) The hard cutoff
breaks the Ward identity in the finite parts with either measure — the finite-part analogue of its Λ⁴ failure; a non-covariant
regulator leaves a Λ-independent non-transverse local polynomial, removed by non-covariant counterterms exactly like
the Λ⁴ piece.
(iv) The "−1/12 discipline" (the author's phrase, 1 Nov: a regulated finite number means something only with the
subtracted piece shown beside it and the unobservability of the subtraction argued), applied: every finite number above is written with its regulator, its measure, the subtracted
piece (c_L log Λ²) and the test of whether the remainder is a local covariant polynomial; the ones that pass are
(c) and (b1)+χ, and their difference is γ_E, unobservable.

### 8.3 The 1/p⁴ resummation

Per spin sector, coefficient of h·P·h in 2Q, 1/16π² restored, massless scalar, Euclidean, MS-bar (μ² ≡ Λ²) /
proper-time-χ / cutoff differing only in the spin-2 constant k (for the cutoff, k is the TT projection of a
non-transverse K):

    G₂⁻¹(p²) = (M_P²/4) p²  +  (1/3840π²) p⁴ [ log(p²/Λ²) − 46/15 (+ γ_E for proper time; −37/30 in place of −46/15 for the cutoff) ]
    G₀⁻¹(p²) = −(M_P²/2) p²  −  (3/4)(ξ − 1/6)² (2/16π²) p⁴ log(Λ²/p²) + …     (no spin-0 running at ξ = 1/6)

P5 as corrected before the first execution of the script (the header quotes the first wording verbatim and records the
sequence): the induced spin-2 p⁴ coefficient carries the asymptotic-freedom sign — it grows with p² — and is therefore
*negative* for p² < Λ² e^{46/15}. What the log tests is only a > 0 and b > 0, which are P4 and P3 again; the rest is a
formula printed and read:
(i) for M_P² ≪ p² ≪ Λ², G₂ ≈ 1/[(1/3840π²) p⁴ log(p²/Λ²)]: the abstract's "∼ 1/p⁴", with a logarithmically running
coefficient — Stelle's 1/k⁴ fifteen years early, not as a pole but as a running term; (ii) with no bare C² term,
G₂⁻¹ has a **Euclidean** zero at p² = 960π² M_P²/(log(Λ²/p²) + 46/15), i.e. p ≈ 30–50 M_P for Λ = 60–4000 M_P (hand
evaluation of the printed formula; self-consistent in Λ only above that, and at the zero the one-loop term equals the
tree term by definition, so the location sits at the edge of the one-loop truncation unless N_S is large), a spacelike
pole = a Lorentzian tachyon, not Stelle's ghost; Stelle's
ghost (negative-residue pole at m² = (M_P²/4)/B) is what a *positive* bare Euclidean C² coefficient B at the scale p
gives when it dominates the negative log; N_S scalars divide the pole position by N_S. Which of the two the 1962
abstract meant by "semiclassical Green's function ∼ 1/p⁴" cannot be decided from the abstract; both are 1/p⁴ at large p.
Nothing here is a statement about the quantum gravity theory — it is the matter loop dressing a classical propagator,
which is all 1962 had — and it is the bridge to rung 3 (the graviton loop lands in the same two slots) and to 1977.

### 8.4 Still open from §6, after this

Spinor and vector (abelian) matter: not done; a decision for the author (brief §6). The 1962 text: unseen. **Rung 1 closed as to computation, 1 Nov
night** (the 31 Oct criterion also listed the spinor/vector decision; the author left that as a standing option rather than a condition, 1 Nov in conversation): the four logs of §8 were run on the author's machine and committed (`results/rung1_cutoff_{plain,covariant}_unshifted.log`,
`results/rung1_finite_{plain,covariant}.log`); the three with sandbox counterparts are identical to them except for
timings; the covariant unshifted run is in. Spinor/vector matter remains the author's decision; Pauli–Villars was
excluded by the author's vintage line (author, in conversation, 1 Nov: "a cutoff on methods has to exist somewhere"), so what a covariant 1962-era
regulator would have given for the ∞⁴, ∞² assignments stays unanswered by design. The sandbox smoke logs were removed in
the same housekeeping patch as this paragraph.

History: §8 added 1 Nov evening; finding 1 of §4 corrected the same evening (the "only there" was wrong at the finite p⁴
order; found by the {Ric², R²} fit and the Ward test of `rung1_finite.py`, not by rereading). Fresh-agent review of
1 Nov night folded in before the patch: Ward test made bilinear; status-label history stated; P5's first wording
quoted verbatim in the header; −46/15 flagged everywhere; "closed" qualified; the 133/10 clause removed.
