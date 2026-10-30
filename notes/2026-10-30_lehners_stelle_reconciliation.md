# Lehners–Stelle 1909.01169 against `bianchi9/`

Their numerical claim cannot be placed on the H/M axis: the paper gives no value for the Weyl coupling σ, no time range, and no unit system, so the Fig. 1 run (BKL persists, crunch faster) is a qualitative statement with no stated ghost mass. Nothing in it contradicts the `bianchi9/` closure; nothing in it can be checked against it either until the run is redone at known H/M. The 2023 recap ([2312.14048](https://arxiv.org/abs/2312.14048)) adds no numerics.

## Convention map

Their action (eq. 1): S = ∫√−g \[ R/κ² − Λ − C²/(2σ) + ωR²/(3σ) \], κ² = 16πG. Yours: R/2 − C²/(4M²), M̄\_P = 1.

| quantity | Lehners–Stelle | `bianchi9/` | agravity (Salvio review) |
| --- | --- | --- | --- |
| Einstein term | 1/κ² | 1/2, so κ² = 2 | M̄\_Pl²/2 |
| Weyl coefficient | −1/(2σ) | −1/(4M²), so σ = 2M² | −1/(2f\_2²), so f\_2² = 2M² |
| spin-2 mass² | σ/κ² (their §II) | M² | M\_2² = f\_2²/2 |
| R² coefficient | ω/(3σ) | 0 | 1/(6f\_0²), so ω = M²/f\_0² |
| spin-0 mass² | σ/(2κ²ω) | ∞ | M\_0² = f\_0²/2 |

The three conventions agree on the spin-2 mass (σ/κ² = 2M²/2 = M²), so M in your plots is their √(σ/κ²) and agravity's M\_2 with no factors. Their fixed-point value ω\* ≈ −0.0228 with σ\* = 0 is RG context, not a run input.

## Their Bianchi IX setup

Metric ds² = −dt² + Σ (l\_m/2)² σ\_m² with l₁ = a e^{(β₊+√3β₋)/2}, l₂ = a e^{(β₊−√3β₋)/2}, l₃ = a e^{−β₊}: your Misner variables up to the factor ½ in the one-forms and a = e^α. Fig. 1 data: a(0) = 100, ȧ(0) = −1/3, β₊(0) = 1, β₋ = 0; the higher-derivative initial data are taken from the two-derivative solution (your Einstein-branch seeding, in effect). β₋ = 0 with β₊(0) = 1 is the head-on θ = π geometry of your sweeps. H₀ = ȧ/a = −1/300 in whatever units κ carries, and with σ unspecified H₀/M is unknown.

Their only quantitative statements are the linearized small-anisotropy equations (11)–(13), with only 1/κ² and 1/σ appearing (no R²), and a dominant-balance analysis about a ∝ t^{1/3}: β₊ ∝ t^b with b ∈ {0, 2/3, 4/3, 2}, all decaying toward t → 0 except b = 0. Then the finite-action argument: on power-law ansätze the Weyl action converges only for s > 1 (accelerating), so BKL s = 1/3 is excluded by postulate, not by dynamics.

## Where the two analyses touch

1. **Exponents.** Their b ∈ {0, 2} are your pure-Weyl indicial roots s³(s−2)³; their 2/3 and 4/3 come from the curvature terms 8β/a² and 64β/a⁴ they keep. But their base is small β about isotropic a ∝ t^{1/3}, which is not a vacuum Bianchi IX solution, while your integrated exponents (δβ ∝ t^{−1} toward t → 0, physical GR mode t^{−0.3}) are about exact Kasner, where Q ≡ 0 forces the exponents to the next order. Not a contradiction; a different linearization. Worth one sentence in the paper, because a reader who has seen their b-list will ask why yours is not in it.
2. **Action divergence.** On exact Kasner C² ∝ t^{−4} and √−g ∝ t, so ∫C²√−g dt ∝ ∫t^{−3}dt diverges at t → 0: their eq. 17 result, reproduced trivially on your base flow. Their exclusion of Mixmaster is a selection principle on top of the same dynamics you integrate; it is orthogonal to the repeller question and can be cited as such.
3. **Chaos.** They write that the behaviour remains chaotic near the big bang, from one run with unstated σ. Your Answer 4 says no billiard at H ≫ M, resting on Cotsakis et al. 1993 and Querella. If their run had M large compared with every H reached (plausible if they picked σ for a small correction), both are right: their run never left your adiabatic regime. If it did reach H ≳ M, their figure would be the first published evidence against your Answer 4 and against Cotsakis et al. Only a rerun at known H/M settles it.

## Suggested actions (your call)

- [ ] Rerun their Fig. 1 data in `bianchi9_seed.py` geometry (θ = π, a(0) = 100, ȧ(0) = −1/3, β₊(0) = 1, β₋ = 0, Einstein-branch higher data) for M chosen so that H₀/M = 0.01, 0.1, 0.3, 1; record where in the collapse H/M crosses 0.2 and whether the anisotropy panel still shows oscillations past it. One Radau run each, seconds.
- [ ] Check their eq. (11) against the small-β limit of `L` in `bianchi9_system.pkl` with 1/σ = 1/(2M²), κ² = 2, a = e^α. An independent derivation of the same minisuperspace equations is free validation of `bianchi9_derive.py`, and it closes the README's open item on comparing the literature Lagrangian with `L`. Caveat: eqs. (11)–(13) were transcribed from a PDF text layer here and may be garbled; read them from the PDF.
- [ ] Paper, outlook item (v): one clause noting that 1909.01169 evolved this system numerically without stating σ, and that its Mixmaster exclusion is a finite-action postulate orthogonal to the dynamical question.
- [ ] Castardelli dos Reis, Chapiro & Shapiro [1903.01044](https://arxiv.org/abs/1903.01044): their Bianchi I result that the sign of the Weyl coefficient decides growth vs. no explosive instability should be matched to your sign (−C²/(4M²), pole at M, non-tachyonic). If their stable sign is yours, cite it as the Bianchi I precursor of your adiabatic regime.

Sources: [1909.01169](https://arxiv.org/abs/1909.01169) (PDF read in full for §§II–III and Fig. 1), [2312.14048](https://arxiv.org/abs/2312.14048) (abstract and § structure), `bianchi9/README.md` at hdym 20fdcf2.
