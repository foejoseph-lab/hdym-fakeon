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

1. **The measure matters at ∞⁴ and ∞², and only there.** With φ as the integration variable the Λ⁴ term is
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
