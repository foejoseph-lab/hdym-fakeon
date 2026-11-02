# Rung 2, the 1971 Lorentz-invariance crack: paper summaries for the author (2 Nov 2026)

Rule 5: each paper summarised to the author before anything is built — the setup, the computation pinned down so two people
could do it independently, what counts as agreement. Open sources only (rule 18). Provenance marks: **[read: …]** read on an open
page tonight; **[abstract/TOC]** title, pages, date from the APS table of contents, text unseen; **[m]** recollection, not verified
tonight. The four 1971 PRD papers are paywalled and **unseen**; what follows about their content comes from (i) Nakanishi's own
open 2006 account, hep-th/0609206, which restates the 1971 computation with equations, and (ii) Anselmi–Piva 1703.04584's account.

## 0. The papers, verified as listings

| paper | listing | content source |
|---|---|---|
| Nakanishi, "Lorentz Noninvariance of the Complex-Ghost Relativistic Field Theory", Phys. Rev. D 3, 811, published 15 Feb 1971 | [read: journals.aps.org/prd/issues/3/4] | hep-th/0609206 §1, §3 |
| Lee & Wick, "Questions of Lorentz Invariance in Field Theories with Indefinite Metric", Phys. Rev. D 3, 1046, 15 Feb 1971 (same issue) | [read: same TOC] | hep-th/0609206 §1; 1703.04584 |
| Nakanishi, "Remarks on the Complex-Ghost Relativistic Field Theory", Phys. Rev. D 3, 3235, 15 June 1971 (Comments section) | [read: journals.aps.org/prd/issues/3/12] | hep-th/0609206 §1 |
| Gleeson, Moore, Rechenberg & Sudarshan, Phys. Rev. D 4, 2242 (1971) | [m] (not re-listed tonight) | hep-th/0609206 §3 ("GMRS", their eq. (A3) with a factor ½ corrected) |
| Cutkosky, Landshoff, Olive & Polkinghorne, Nucl. Phys. B 12, 281 (1969) — CLOP | [m] | 1703.04584 §4 |
| Lee & Wick, Nucl. Phys. B 9, 209 (1969); B 10, 1 (1969); Phys. Rev. D 2, 1033 (1970) | [read: JHEP06(2017)066 reference list] | 1703.04584 §2 |
| Anselmi & Piva, JHEP 06 (2017) 066, 1703.04584 | [read: arXiv PDF v2] | — |
| Nakanishi, "Wightman-function approach to the relativistic complex-ghost field theory", hep-th/0609206 (2006) | [read: arXiv PDF] | — |

## 1. The 1971 computation, pinned down (from Nakanishi's own restatement, hep-th/0609206 §2–3 [read])

**Theory.** Two hermitian scalars φ₁, φ₂ with opposite-sign kinetic terms and a mixing term −γ m φ₁φ₂ (his eq. 2.1): the
indefinite-metric Lagrangian whose normal modes are a complex-conjugate pair of masses, M² = m² + iγm and M*² (stated after his
2.13). No complex mass appears in the action; the pair appears in the propagators Δ_F(ξ; M²), Δ_F(ξ; M*²). This is the one-mode
version of the Pais–Uhlenbeck/Lee–Wick structure that `bubble/` already uses (the `bubble_lib.py` propagator pair) — check that
the two conventions map before anything is run.

**Contour.** For the M² propagator the k₀ contour C must pass *above* the pole at k₀ = +√(M² + k⃗²) even though that pole is in the
upper half plane (and symmetrically for the others): the Lee–Wick prescription. Nakanishi justifies the needed time integrals through
a "complex δ-function" δ_c(k₀ − λ) (his 3.1–3.2), defined by a contour integral around the complex λ.

**The diagram.** The complex-ghost-pair bubble, his eq. (3.3):

    F_MM*(p) = [1/((2π)⁴ i)] ∫d³q ∫_C dq₀  1 / { (q² − M²) [(p − q)² − M*²] },   UV divergence subtracted,

external p = (p₀, p⃗), **real spatial loop momentum q⃗**, q₀ along C. Poles in q₀: −E_q and p₀ − E*_{p−q} below C; E_q and p₀ + E*_{p−q}
above C, with E_q = √(M² + q⃗²).

**Result.** After the q₀ integration, his eq. (3.4):

    F_MM*(p) = [1/(2(2π)³)] ∫d³q  (1/E_q + 1/E*_{p−q}) / { p₀² − (E_q + E*_{p−q})² }.

The would-be unitarity cut is the locus p₀ = E_q + E*_{p−q} (3.5); as q⃗ ranges over R³ it sweeps a two-dimensional "fish-shaped"
domain in the complex p₀ plane whose boundary touches the real axis at the single point p₀ = 2 Re E_{p⃗/2}. Hence the threshold
(3.6):

    s_min ≡ 4 (Re E_{p⃗/2})² − p⃗²,    and   s_min < (M + M*)²  for p⃗ ≠ 0,

with s = p₀² − p⃗². **This is the crack:** s_min depends on p⃗ through E_{p⃗/2} = √(M² + p⃗²/4), not on s alone — the threshold of a
Lorentz scalar amplitude moves with the frame. GMRS (1971) gave the closed form, his eq. (3.7): F = (1/16π²) Re ∫_{Γ_s} ds'/(s − s')
× [ |p⃗|(M² − M*²)/(s'√(s' + p⃗²)) + √((s' − (M+M*)²)(s' − (M−M*)²))/s' ], Γ_s the lower half of the image of the fish boundary under
s' = p₀² − p⃗², from s_min to ∞ (a factor ½ in GMRS's (A3) corrected by Nakanishi).

At p⃗ = 0: E₀ = M and s_min = 4(Re M)² = (M + M*)² — the frame-independent point; this is the case `bubble/` has done (the instrument check,
chronology §5 rule 22). The p⃗ ≠ 0 computation is the rung.

## 2. The disagreement, in a few sentences each

- **Nakanishi (Feb 1971, PRD 3, 811):** with real spatial momenta and the Lee–Wick energy contour, the complex-ghost-pair
  amplitude is not Lorentz invariant already at one loop (the threshold above); unitarity was bought at that price. [hep-th/0609206 §1]
- **Lee & Wick (Feb 1971, PRD 3, 1046, same issue):** the violation comes from pairing complex energies with *real* spatial momenta;
  the spatial momenta should be complexified too, and the coalescing-singularity ambiguity handled by the CLOP rule. [hep-th/0609206 §1;
  1703.04584 §1–2]
- **Nakanishi (June 1971, PRD 3, 3235) and his later view:** the complex energy is met only on the mass shell, while the theory is
  built off shell; the breaking is intrinsic to the Hamiltonian/indefinite-metric formulation with complex energies, not to the choice of
  real momenta. In the 2006 paper he argues the fundamental fields' two-point Wightman functions are Lorentz invariant while the
  *pair* (composite) Wightman function is not — the violation sits in composite states. [hep-th/0609206 §1, §2, §4]
- **Gleeson–Moore–Rechenberg–Sudarshan (1971, PRD 4, 2242):** confirmed Nakanishi's result explicitly and gave the closed form (3.7).
  [hep-th/0609206 §3]
- **CLOP (1969):** in 1703.04584's description, pinchings are resolved by a limit: shift M → M + δ in one propagator, compute, let
  δ → 0; this leaves a ±π ambiguity (their 4.15) that appears only for unequal masses m₁ ≠ m₂. [1703.04584 §4]

## 3. How 1703.04584 restates and resolves it [read]

**Model (their 2.1).** A single LW scale M, real: iD(p², m², ε) = iM⁴/[(p² − m² + iε)((p²)² + M⁴)], six poles per propagator: the
standard pair p⁰ = ±(ω_m − iε) and two LW pairs p⁰ = ±Ω_M, ±Ω_{M*}, Ω_M = √(p⃗² + iM²) (their 2.2, fig. 1). Decomposition (3.9):
D = D₀ + D_LW with D₀ = [M⁴/(M⁴+m⁴)]/(p² − m² + iε), D_LW = −[M⁴/(M⁴+m⁴)](p² + m²)/((p²)² + M⁴). Note the difference from
Nakanishi's model: here the LW pair sits at p² = ±iM² (purely imaginary mass-squared), not at m² ± iγm; the physics of the crack is
the same, the numbers are not — any comparison with §1 must map the models first.

**Diagram (their 2.3, fig. 3).** J(p) = ∫d^Dk/(2π)^D D(k², m₁²) D((k−p)², m₂²), both propagators LW-type, same M, masses m₁, m₂;
amplitude −iλ²J/2. The LW-pinching part is J_LW (3.10). Numerics in their §5: M = 1, m₁ = m₂ = 3 (and 3, 5; 4, 4), chosen so the
physical threshold does not overlap the LW one; p⃗ = 0 with their formula (4.8), and p_s = 10⁻³ as the regulator in the new formulation.

**Definition.** Start Euclidean; the Wick rotation to Minkowski is **non-analytic**: it gives the LW theory only in the "main region"
Ã₀ of the complex p⁰ plane containing the imaginary axis; other regions (Ã_P around P at p⁰ = E_P = √(2M² + p⃗²), and its mirror Ã_P')
are reached by a definite, non-analytic procedure (fig. 12 for the final regions A₀, A_P, A_P').

**The rule.** Two routes are described: (§3.1) keep the spatial integration domain rigid and add the residues of the poles that cross
it; (§3.2, preferred) **deform the spatial loop-momentum domain (k_s, u) into complex values**, in complex-conjugate ways for the two LW
pairs, with the k⁰ integral done by residues on the imaginary axis. Lorentz invariance is recovered region by region: the k⃗-domain is
deformed until the regions are squeezed onto Lorentz-invariant lines, then J_LW^> = ½[J_LW^{0+} + J_LW^{0−}] (their 3.11), the two
terms being the analytic continuations of the main-region function from Im p⁰ > 0 and < 0 (fig. 11) — the **average continuation**
in its first appearance.

**Pinching condition (their 3.7/4.2).** p⁰ = √(k⃗² + iM²) + √((k⃗ − p⃗)² − iM²) and its three partners (sign of iM² in the second
root; complex conjugates). Real solutions with real k⃗ (4.3): k⃗² = ((p⁰)⁴ − 4M⁴)/(4(p⁰)²) and p⃗² = 2 p⃗·k⃗. At p⃗ = 0: k_s² = (p⁰)²/4 − iM²
(3.5). The LW threshold is p² = 2M² (general form p² = (√2 nM + m_{i₁} + …)², their §6); complex branch points at p² = ±4iM² at p⃗ = 0.

**Result.** For m₁ = 3, m₂ = 5 CLOP gives two different answers depending on the sign of δ; the new formulation gives the average of
the two (fig. 14). For m₁ = m₂ = 4 they agree closely. The paper does not give a closed form for J^> in every region; it gives the
pinch-level formulas (4.8, 4.12–4.15) and numerics. An equal-mass correction when p⁰ crosses the real axis, §3.1:
J → J − (1/16π)·M⁴/(m⁴+M⁴)·√(1 − 4M⁴/(p²)²)·θ₋(p² − 2M²).

Also [read, abstract-level]: 1703.05563 (perturbative unitarity, same formulation); 1612.06510 (Aglietti–Anselmi): integrating k⁰
along the real axis in a higher-derivative theory produces nonlocal, non-Hermitian divergences — 1703.04584 §1 cites this as the
reason a Minkowski-first definition is not viable. Their §3 details were not read tonight.

## 4. Proposed "what counts as agreement" for the p⃗ ≠ 0 run (for the author to accept or change)

The computation of §1 (Nakanishi's F_MM*) is the historical one and is elementary enough to do in `bubble/` with method P and
method M (rule 3) — it is one 3-dimensional integral after the q₀ residues, with a known closed form (3.7) to check against:

1. **Technique of the day (1971): Lee–Wick contour, real q⃗.** Compute F_MM*(p₀, |p⃗|) for a grid of frames at fixed s = p₀² − p⃗²,
   e.g. s/(M+M*)² ∈ {0.8, 0.95, 1.0, 1.05, 1.5} and |p⃗|/|M| ∈ {0, 0.5, 1, 2}. **Prediction to write before the run:** the threshold
   s_min(|p⃗|) = 4(Re E_{p⃗/2})² − p⃗² falls below (M + M*)² as |p⃗| grows, so at fixed s just below (M+M*)² the imaginary part is zero at
   p⃗ = 0 and nonzero at large |p⃗|; the frame dependence is in the threshold, not only in the value. Agreement with GMRS (3.7) to the
   numerical precision of the integral is the check of the implementation. This reproduces the crack as a number.
2. **Next technique (1969, CLOP):** the same with M → M + δ in one propagator, δ → 0 from both sides; expect the ±π ambiguity for
   unequal masses and none for equal masses (1703.04584 §4).
3. **2017, region-wise continuation (1703.04584 §3.2):** in their model (2.1). Expect J to depend on s only, region by region, with the
   LW threshold at p² = 2M² and the equal-mass correction above.
4. **2018, average continuation (1801.00915):** the same model; expect equality with 3 where both are defined (the paper's 3.11 *is*
   the average), and a definite answer for unequal masses where CLOP had two.
5. **1612.06510, Minkowski real-axis k⁰:** the control (rule 13): expect the nonlocal divergence, logged as a result.

"Agreement" = the 1971 computation reproduces (3.6)–(3.7) and *is* frame dependent; the 2017 and 2018 computations are frame independent
in each region and coincide with each other; CLOP coincides with them for equal masses only. A disagreement anywhere else is a finding
(rule 1). The one open modelling choice is whether to run the 2017/2018 prescriptions on Nakanishi's model (M² = m² + iγm) or
Anselmi–Piva's (p² = ±iM²): the chronology says technique of the day first, so the model should be the 1971 one and the later
prescriptions applied to it, with their own model run as the certification against their fig. 14. Author's call.

## 5. Provenance summary

Read tonight: the two APS tables of contents (titles, pages, dates of PRD 3, 811 / 1046 / 3235); hep-th/0609206 (Nakanishi's own
restatement, eqs. 2.1, 3.1–3.7, §1 history); 1703.04584v2 (model, diagram, regions, rule, pinching conditions, numerics, figures);
the JHEP reference list for the 1969–70 Lee–Wick references and Modesto–Shapiro. Unseen: the four 1971 PRD texts, CLOP 1969, the
1969–70 Lee–Wick papers, 1612.06510 §3 and 1703.05563 beyond their abstracts. From memory [m]: that GMRS confirmed Nakanishi
(stated by Nakanishi 2006, so [read] at one remove); the CLOP reference itself.
