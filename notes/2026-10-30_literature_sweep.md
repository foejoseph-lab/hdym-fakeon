# Fakeon spin-2 literature sweep (30 Oct 2026)

Exported from the session doc of 30 Oct 2026 (three parallel sweeps, arXiv pages fetched; INSPIRE and the arXiv API were unreachable from the sandbox, so the negative results rest on 47 web searches plus the fakeon papers' own statements of scope). The reconciliation with Lehners–Stelle and the step-1 toy results are separate notes in this folder.

## Verdict

The target is unaddressed: no one has applied the fakeon classical projection, or the Lee–Wick or Merlin treatments, to any background with a nonvanishing Weyl tensor, let alone to Mixmaster. The spin-2 fakeon χ has only ever been projected as a linear tensor perturbation about de Sitter; on FLRW it drops out identically. The covariant projected equations exist on paper (1809.05037) and have never been solved. Anselmi–Marino and Salvio both say in print that generic backgrounds are open.

Three things to pin down (written 30 Oct under a classical plan since set aside; they remain the first things a referee will ask of any route):

1. **Reconcile `bianchi9/` with Lehners & Stelle 1909.01169**, who evolved EH + C² + R² on Bianchi IX and report BKL behaviour persisting with a faster crunch, where `bianchi9/` reports the Weyl ghost ending the chaos. Likely a matter of where in H/M their runs sit relative to the 0.2 boundary and of their not computing chaos indicators, but it has to be shown, and it is the first thing a referee will raise (section 6; separate note).
2. **The spin-2 projection recipe is only spelled out for a linear mode about de Sitter** (ABP20's U/V split, section 3). Its nonlinear minisuperspace form, with C ≠ 0 and the constraints, does not exist in the literature and is the main technical risk; the one-mode toy (section 7) is the place to catch sign and ordering errors before Bianchi IX.
3. **The ABP20 bound is a reality condition on the kernel, m_χ² > H²/4, not a convergence statement.** It gives a sharp prediction for where the projected dynamics should change character, H/m_χ = 0.5, that sits above the unprojected boundary H/M ≈ 0.2. Whether the projected Mixmaster boundary lands at one, the other, or neither is the result.

The convergence of the projection has never been studied beyond: exact for special FLRW sources, asymptotic with factorially growing coefficients for dust, and a finite initial-condition count in flat space (Anselmi & Calcagni, JHEP 2026). A chaotic base flow is the natural stress test and nobody has run it.

The comparison with agravity is clean: same action, both theories keep φ as a normal scalar in the GSF variant, and the only difference is what is done with χ. Salvio's metastability paper already contains a Bianchi I-type run with the unprojected ghost and already names the fakeon approach's FLRW-only status as a limitation. Donoghue–Menezes hold that there is no classical statement to make about the ghost at high curvature. A simulation of the projected ghost on Mixmaster is evidence in a three-way argument that is currently being had without any.

## 1. The fakeon prescription (Anselmi & Piva 2017–2018)

A fakeon is a degree of freedom that lives only on internal lines and is removed from the asymptotic spectrum by the quantization rule itself, not by assuming it decays. That is the whole difference from a Lee–Wick ghost, which is dropped from external states only because it is taken to be unstable ([2202.10483](https://arxiv.org/abs/2202.10483)).

The rule: amplitudes are defined in Euclidean signature and continued to Minkowski by a nonanalytic "average continuation" ([1801.00915](https://arxiv.org/abs/1801.00915)); at tree level the fakeon propagator is the principal value, the half-sum of the +iε and −iε propagators, which in position space is the half-sum of retarded and advanced Green functions:

    G_F(p) = lim_{ℰ→0} (p² − m²)/[(p² − m²)² + ℰ⁴] = ½(G₊ + G₋) = PV 1/(p² − m²)

Inside loops the principal value is not used naively (that is the Feynman–Wheeler propagator, which breaks renormalizability and unitarity, [2001.01942](https://arxiv.org/abs/2001.01942)); the average continuation of the whole amplitude is what preserves both. Unitarity to all orders and ordinary power-counting renormalizability are the two results of the 2017–2018 papers ([1703.04584](https://arxiv.org/abs/1703.04584), [1703.05563](https://arxiv.org/abs/1703.05563), [1704.07728](https://arxiv.org/abs/1704.07728)).

For gravity ([1806.03605](https://arxiv.org/abs/1806.03605)) the multiplet is the graviton h, a massive scalar φ from R² (mass m_φ² = ζ/ξ) and a spin-2 field χ from C² (mass m_χ² = ζ/α). χ is always the fakeon; φ may be physical (the "GSF" theory) or a fakeon too ("GFF"). The one-loop self-energy gives χ a negative width, read as microcausality violation above m_χ. Stelle gravity with χ a genuine ghost is the explicit foil in that paper, and agravity with both modes kept is the one for us.

## 2. The classical limit: the fakeon projection

The classical fakeon theory is the "classicized" action: the sum of tree diagrams with physical external legs and fakeon internal lines. The local higher-derivative action is only an "interim" action; the finalized classical equations are nonlocal and of lower order ([1809.05037](https://arxiv.org/abs/1809.05037)).

**The recipe.** Solve the fakeon's own field equation with the averaged Green function, half retarded plus half advanced, and substitute back. In the toy (a higher-derivative nonrelativistic particle with time scale τ) the kernel is G_F(u,τ) = sin(|u|/τ)/(2τ), supported on both past and future and tending to δ(u) as τ → 0; the projected equation is Newton's law with a time-averaged force, m q̈ = −∂V/∂x evaluated at ⟨q⟩. The runaway solutions of the interim equation are not perturbative solutions and drop out; only position and velocity survive as data.

For gravity the average of an operator O is built from 1/(m² + ∂²) with retarded and advanced parts weighted m²/2 each (eq. 5.9 there). The finalized equations are written for general metrics, with the spin-2 sector included:

    GFF:  R_μν − ½ g_μν R = (1/M̄_Pl²) T^GFF_μν,   T^GFF_μν ∋ ⟨T_μν − U_μν⟩_χ + (trace term with ⟨T − U⟩_φ)
    GSF:  G_μν + (1/3m_φ²)(g_μν ∇² − ∇_μ∇_ν) G = (1/M̄_Pl²) ⟨T_μν − U_μν⟩_χ

Here U_μν collects terms at least quadratic in h around flat space, g = η + 2κh. The expansion is in κ, not in α or ξ; a nonperturbative treatment of the potential would need a resummation the paper does not give. **No convergence analysis exists** beyond the remark that the kernel's effects average away for |u| ≫ τ.

**FLRW ([1901.09273](https://arxiv.org/abs/1901.09273)).** The Weyl tensor vanishes, so χ drops out identically and only the scalar sector is projected. With Σ = 1 + m_φ⁻²(3ȧ/a + d/dt)d/dt and Υ = Σ + 2m_φ⁻²(k/a² + 3 d/dt(ȧ/a)), the projection multiplies the fourth-order equations by Σ⁻¹, Υ⁻¹ defined by the fakeon prescription, giving Friedmann equations with averaged sources ρ̃ = ¼⟨ρ−3p⟩_Σ + ¾⟨ρ+p⟩_Υ. The projection is exact for vacuum energy, radiation and their sum; for dust it is only asymptotic, the 1/(m_φ t)² series having a coefficient of order 10⁸ at order 14. The classical theory inherits the quantum one's inability to write exact complete field equations. This is the closest thing in the literature to a convergence statement, and it is on the easy background.

**Constant-H FLRW ([1909.12873](https://arxiv.org/abs/1909.12873)).** The fakeon Green function picks up a θ(Ht) factor, so for H > 0 only the past contributes and microcausality is restored. The paper says outright that the prescription is understood around flat space and that generic backgrounds are open.

**Initial-condition count ([2510.05276](https://arxiv.org/abs/2510.05276), Anselmi & Calcagni, JHEP 2026).** The nonlocal classicized equations do not need infinitely many initial data: their solution space is a subspace of the interim local theory's, and the count drops (Dirac radiation reaction 3 → 2; a two-derivative tachyonic fakeon 4 → 2). The fakeon prescription must be reached as a limit of integration bounds a → −∞, b → +∞ to avoid overcounting. Divergent nonlocal integrals on runaway solutions are called an artifact; there is still no convergence theorem and no background beyond flat-space linear toys. Calcagni's diffusion-method follow-up ([2510.15157](https://arxiv.org/abs/2510.15157)) reproduces the count.

**What is and is not available for spin 2.** The covariant ⟨T − U⟩_χ equations exist on paper. They have never been solved on any background with C ≠ 0. The only place χ is actually projected is in tensor perturbations about de Sitter (section 3), where it is linear and decouples at quadratic order.

## 3. Fakeon cosmology and the ABP20 bound

The bound is m_χ > m_φ/4, eq. (4.9) of Anselmi, Bianchi & Piva ([2005.10293](https://arxiv.org/abs/2005.10293), JHEP 2020), equivalently 1 < m_φ²/m_χ² < 16 as restated in the review [2203.02516](https://arxiv.org/abs/2203.02516). It is a no-tachyon condition on the projected tensor fakeon, not a convergence bound.

**Where it comes from.** The tensor modes of R + R² + C² about quasi-de Sitter have a fourth-order Lagrangian. An auxiliary field S and the redefinitions u = U + αV, S = V + βU split a physical mode U from a fakeon V. V's equation is projected with the averaged Green function, which on de Sitter solves

    (d²/dt² + m_χ² − H²/4 + k²/a²) Ĝ = δ,   n_χ = √(m_χ²/H² − 1/4)

and is fixed as half retarded plus half advanced by matching the flat-space limits k/(aH) → ∞ and → 0. Superhorizon, Ĝ ∝ sin(H n_χ |t − t′|)/(2H n_χ): real and oscillatory only if n_χ is real, i.e. m_χ² > H²/4. With the inflationary de Sitter limit H = m_φ/2 this is m_χ > m_φ/4. A stricter optional criterion (positivity of the superhorizon effective mass for all k/aH) would give m_χ ≳ 0.312 m_φ, eq. (4.15), and is not adopted.

**What happens below the bound** is stated only as: n_χ imaginary, the oscillatory kernel is lost, the superhorizon effective mass turns tachyonic, and the nonlocal fakeon contribution is no longer negligible; a full analysis is deferred. No power-law or runaway statement is made. That is the gap a Bianchi IX calculation would speak to directly, since Mixmaster is where H/m_χ is not small and the kernel's reality condition is tested dynamically rather than at a fixed H.

**Predictions at leading order**, for the parameter map: n_R − 1 = −2/N; r = 24 m_χ² / [N²(m_φ² + 2m_χ²)]; n_T = −3 m_χ² / [N²(m_φ² + 2m_χ²)], so r ≃ −8 n_T and 4/3 < N² r < 12. m_χ → ∞ recovers Starobinsky. The RG-flow reformulations ([2007.15023](https://arxiv.org/abs/2007.15023), [2010.04739](https://arxiv.org/abs/2010.04739), [2103.01653](https://arxiv.org/abs/2103.01653), [2105.05864](https://arxiv.org/abs/2105.05864)) stay on FLRW and add nothing on the background-level projection of χ.

**Parameter range for the comparison.** Any run that wants to be inside the fakeon theory's own consistency window should take m_φ/4 < m_χ < m_φ; agravity has no such window, which is itself a point of comparison. In agravity's conventions (section 4) m_χ ↔ M_2 = f_2 M̄_Pl/√2 and m_φ ↔ M_0 = f_0 M̄_Pl/√2, so the window reads f_0/4 < f_2 < f_0.

## 4. Agravity and the Salvio–Strumia line on the Weyl ghost

Agravity ([1403.4226](https://arxiv.org/abs/1403.4226)) keeps both higher-derivative modes as genuine degrees of freedom and treats the ghost's negative norm as a postponed problem. That is the comparison: same action, no projection.

**Action and masses** (Salvio's review [1804.09944](https://arxiv.org/abs/1804.09944), signature + − − −):

    L = R²/(6f_0²) − W_μνρσ W^μνρσ/(2f_2²) − (M̄_Pl²/2) R − Λ,   W²/2 = R_μν² − R²/3 + G/2
    M_2 = f_2 M̄_Pl/√2,   M_0 = f_0 M̄_Pl/√2 + …,   f_0² = 1/(2β + 2γ + 6α),   f_2² = −1/(β + 4γ)

for L = αR² + βR_μν² + γR_μνρσ². No tachyons needs f_0², f_2² > 0. In agravity proper the M̄_Pl² R term is absent and induced by ξ_S⟨s⟩². The `bianchi9/` action S = ∫√−g [R/2 − C²/(4M²)] is the f_0 → ∞, M̄_Pl = 1 corner with M² ↔ f_2²/2 up to signature convention; the R² term is the one addition path B step 3 needs.

**One-loop results** (eq. 34 of agravity): f_2 is asymptotically free, (4π)² df_2²/dlnμ = −f_2⁴(133/10 + N_V/5 + N_f/20 + N_s/60); f_0 is not, (4π)² df_0²/dlnμ = (5/3)f_2⁴ + 5f_2²f_0² + (5/6)f_0⁴ + …. Naturalness wants f_0, f_2 ∼ 10⁻⁸. [1705.03896](https://arxiv.org/abs/1705.03896) argues the f_0 growth is benign: the conformal mode decouples and agravity flows to conformal gravity plus two conformally coupled scalars.

**The ghost, classically: metastability.** Salvio's [1902.09557](https://arxiv.org/abs/1902.09557) is the paper to measure against. Its claims: (i) the Weyl tensor vanishes on FRW, so the ghost is inert there and only anisotropy or inhomogeneity excites it; (ii) runaways are avoided below thresholds E_2 = M_2/√f_2 = √(f_2/2) M̄_Pl and E_m = f_2^{1/4} M̄_Pl, with spin-2 field values below M̄_Pl; (iii) **its section 3 evolves a homogeneous diagonal anisotropic metric with shear variables β±, Bianchi I in all but name, numerically in the Einstein frame**: below threshold the anisotropy decays and ~62 e-folds of inflation follow, above it the scale factor collapses; (iv) linear perturbations about de Sitter are bounded for M_2 > H. No lifetime is computed, and the same paper remarks that the fakeon approach had only been developed on flat and FRW metrics. That remark is the opening.

Gross, Strumia, Teresi & Zirilli ([2007.05541](https://arxiv.org/abs/2007.05541)) make the metastability precise in mechanics: weakly coupled non-resonant ghosts lock down by KAM/Nekhoroshev, with lowest-order Lyapunov time τ_0 = ω_1ω_2/(8λ J_max); cubic resonances can break it; in field theory the infinitely many resonances destroy lockdown and the system drifts, but for Planck-suppressed couplings Γ ∼ T³/M_Pl² ≲ H. The QFT runaway rate is admitted unknown. This is the same island/escape language as hdym's `islands/` and `rabbit_holes/ym/` work, applied to a two-oscillator toy rather than to a flow with a billiard.

Later Salvio papers ([2012.11608](https://arxiv.org/abs/2012.11608), [2404.08034](https://arxiv.org/abs/2404.08034), [2406.12958](https://arxiv.org/abs/2406.12958), [2507.08803](https://arxiv.org/abs/2507.08803)) move to Dirac–Pauli / indefinite-norm quantization and a background-independent Euclidean formulation; none adds a classical anisotropic computation. The only other Bianchi-type quadratic-gravity evolution found in this pass is de Medeiros et al., Bianchi V with tilt ([2403.10377](https://arxiv.org/abs/2403.10377)), treating the theory as healthy and classical.

**Other camps, one line each** (not comparison theories, but reviewers will ask): Holdom & Ren's confinement conjecture, the ghost pole removed nonperturbatively like the gluon's ([1512.05305](https://arxiv.org/abs/1512.05305)); Mannheim's PT-symmetric inner product with Bender ([0706.0207](https://arxiv.org/abs/0706.0207), [1801.03220](https://arxiv.org/abs/1801.03220)), carried to full quadratic gravity by Kuntz ([2410.08278](https://arxiv.org/abs/2410.08278)); Kumar & Marto's sign-flipped Weyl term ([2603.07150](https://arxiv.org/abs/2603.07150)), which lumps Salvio–Strumia, Donoghue–Menezes and Anselmi together as keeping a negative Hamiltonian for the mode.

## 5. Donoghue & Menezes: Merlin modes

The Merlin mode is the heavy pole of the quartic propagator with both residue and width of the wrong sign, iD ≃ iZ/(q² − m² + iZγ) with Z = −1, so it is the time-reverse of an ordinary resonance and carries positive energy backward in time for ∼1/γ ([1908.02416](https://arxiv.org/abs/1908.02416)). Unitarity holds to all orders with cuts through stable states only; the Lee–Wick contour is needed only when the ghost is treated as stable in the narrow-width limit.

**Difference from the fakeon.** Donoghue–Menezes keep the ordinary Feynman iε for every field, let loops give the ghost its wrong-sign width, and drop it from the asymptotic spectrum because it is unstable. Anselmi defines the theory in Euclidean signature and reaches Minkowski by the average continuation, removing χ at every order with no width needed. In "On quadratic gravity" ([2112.01974](https://arxiv.org/abs/2112.01974)) they say Anselmi's prescriptions seem to work and locate the disagreement in which signature is fundamental; the word fakeon is not used, and Anselmi is otherwise only in reference lists.

**Their classical limit is the one that matters here, and it is the opposite of Anselmi's.** In "The Ostrogradsky instability can be overcome by quantum physics" ([2105.00898](https://arxiv.org/abs/2105.00898)) the classical regime is the low-energy, long-wavelength limit at fixed ħ, in which only the light field survives and no runaway appears; taking ħ → 0 in the high-energy propagator is called unphysical because the width is intrinsically quantum. Stability at high curvature is left open in [2112.01974](https://arxiv.org/abs/2112.01974). So for Mixmaster at H ∼ m_χ there are three distinct classical statements on the table: agravity's unprojected ghost (the `bianchi9/` runs), Anselmi's projected nonlocal dynamics, and Donoghue–Menezes's position that there is no classical statement to make there. A simulation that shows what the projected dynamics does at H ≈ 0.2M and above is evidence in that three-way argument.

Supporting papers: gauge-assisted quadratic gravity ([1804.04980](https://arxiv.org/abs/1804.04980)); the single near-physical Lee–Wick pole ([1812.03603](https://arxiv.org/abs/1812.03603)); the arrow-of-causality papers ([1908.04170](https://arxiv.org/abs/1908.04170), [2003.09047](https://arxiv.org/abs/2003.09047), [2106.05912](https://arxiv.org/abs/2106.05912)), where the macroscopic arrow is set by stable states and Merlin modes run backward for about a Planck time.

## 6. Anisotropic backgrounds in quadratic gravity, 2015–2026: the unaddressed check

No paper applies the fakeon, Lee–Wick or Merlin prescription to any anisotropic homogeneous background. The fakeon papers themselves state the FLRW/flat-space restriction (1901.09273, 1909.12873, 2510.05276), and Salvio's metastability paper says so of the fakeon approach from the outside. That part of the goal is unoccupied.

**But Bianchi IX with the Weyl term active is not empty, and one result needs reconciling with `bianchi9/`** (separate note).

| Paper | Year | C² on Bianchi IX | Prescription on anisotropic background | What it finds |
|---|---|---|---|---|
| [Lehners & Stelle 1909.01169](https://arxiv.org/abs/1909.01169) | 2019 | yes, numerical | none | EH + C² + R² on Bianchi IX toward the crunch; the quadratic terms do not drastically alter BKL behaviour but speed up the crunch. The Mixmaster exclusion is a finite-action selection argument, not dynamical. Ghost status left open. |
| [Lehners & Stelle 2312.14048](https://arxiv.org/abs/2312.14048) | 2023 | yes, analytic | none | Kasner then Bianchi IX; on-shell action converges only for accelerating a(t), so BKL a ∝ t^{1/3} is excluded by the finite-action postulate. |
| [Castardelli dos Reis, Chapiro & Shapiro 1903.01044](https://arxiv.org/abs/1903.01044) | 2019 | Bianchi I only | none | Nonlinear numerics of a₁C² + a₂R² on Bianchi I with the anisotropy as a zero-frequency gravitational wave. For a₁ > 0 anisotropies grow; for a₁ < 0 no explosive instability, matching linear theory. The sign of the Weyl coefficient decides. |
| [de Medeiros et al. 2403.10377](https://arxiv.org/abs/2403.10377), [2502.03653](https://arxiv.org/abs/2502.03653) | 2024–25 | Bianchi V | none | Tilted fluid in R² + Ric²; Kasner or isotropic past attractors by density; ghost-induced shear oscillations noted, not analysed. |
| [Salvio 1902.09557](https://arxiv.org/abs/1902.09557) §3 | 2019 | Bianchi I-type, Einstein frame | none | Isotropization and inflation below the energy threshold, collapse above. |
| [Middleton 1007.4669](https://arxiv.org/abs/1007.4669) | 2010 | Kasner in pure C²; IX qualitative | none | Fully anisotropic Kasner-like solutions exist only in pure Weyl²; oscillations may terminate if all exponents can be positive. |
| [Kumar & Marto 2603.07150](https://arxiv.org/abs/2603.07150) v2 §6 (added 31 Oct evening) | 2026 | Kasner/BKL, analytic | sign-flipped Weyl term (an inverted-oscillator mode in place of the ghost; a change of theory, not a prescription on the same theory) | argues BKL/Mixmaster singularities give infinite on-shell action in quadratic gravity, after Lehners–Stelle, and that their mode suppresses shear growth; no numerics. The nearest competitor claim about anisotropic backgrounds in a ghost-free version of quadratic gravity; listed as a one-liner under "other camps" until this row |
| Cotsakis, Demaret, De Rop & Querella, PRD 48, 4595 | 1993 | yes, general quadratic | none | BKL chaos structurally unstable, replaced by a monotonic isotropic power-law approach. Pre-dates every modern ghost prescription. |
| Barrow & Sirousse-Zia, PRD 39, 2187; Barrow & Cotsakis 1989 | 1989 | yes | none | First Mixmaster-in-quadratic-gravity analyses. |
| Demaret, Querella & Scheen [gr-qc/9811086](https://arxiv.org/abs/gr-qc/9811086); Querella thesis [gr-qc/9902044](https://arxiv.org/abs/gr-qc/9902044) | 1998–99 | Bianchi I only for C² | none | Bach-gravity Hamiltonian, Painlevé non-integrability; Bianchi IX in the full fourth-order theory, chaos structurally unstable. |

**The reconciliation item.** The `bianchi9/` finding is that the Weyl ghost ends the Mixmaster chaos (adiabatic below H ≈ 0.2M, power-law repeller above, no exponential runaway). Lehners–Stelle report BKL persisting qualitatively with a faster crunch in their runs; Cotsakis et al. 1993 report chaos destroyed and replaced by monotonic isotropic approach. All three can be consistent, since Lehners–Stelle do not compute chaos indicators or run long oscillation sequences and work at an unstated H/M, but a reviewer will put 1909.01169 next to the claim. The first `bianchi9/`-related task on the new goal is to run their parameters and initial data through that code and state where in H/M their runs sit relative to 0.2. The sign result of 1903.01044 also needs matching against the C²/(4M²) convention: their a₁ < 0 stable case is presumably the ghost-not-tachyon sign, and the paper should say so explicitly.

**Not relevant despite titles:** Gu, Modesto & Bambi [2511.00585](https://arxiv.org/abs/2511.00585) (dilaton-Weyl-invariant GR, not C²); Moriconi, Montani & Capozziello [1411.0441](https://arxiv.org/abs/1411.0441) (R + qR², no Weyl); Beatriz, Bonotto & Lappicy [2603.13104](https://arxiv.org/abs/2603.13104) (Hořava/f(R) attractor theorem); the asymptotic-safety Bianchi I papers [2509.09759](https://arxiv.org/abs/2509.09759), [2604.21237](https://arxiv.org/abs/2604.21237) (Einstein–Hilbert truncation, no ghost).

**Classical ghost-removal by other means, 2021–2026, none anisotropic:** Deffayet, Mukohyama & Vikman's integrable bounded ghosts ([2108.06294](https://arxiv.org/abs/2108.06294), [2305.09631](https://arxiv.org/abs/2305.09631)); Fring & Turner's Toda-lattice benign ghosts ([2301.11317](https://arxiv.org/abs/2301.11317)); Hell, Lüst & Zoupanos's boundary-condition removal of the conformal-gravity ghost on de Sitter ([2306.13714](https://arxiv.org/abs/2306.13714)); Sutton, de Felice & Sakellariadou invoking the fakeon prescription for −αC² on flat FLRW ([2504.15005](https://arxiv.org/abs/2504.15005)). No Pais–Uhlenbeck–Mixmaster paper exists.

**Caveat on the negative result.** INSPIRE and the arXiv API were unreachable in the session; the negatives rest on 47 web searches (query list available on request) plus the fakeon papers' own statements of scope. A ten-minute INSPIRE citation walk from 1901.09273 and 1909.01169 would close that gap.

## 7. Equations carried into the one-mode toy (classical; ħ → 0 check material, not the route)

Written before the toy was built; the toy's own note records what changed (the iteration was replaced by the order-by-order reduction, and the BVP route was found to be the wrong tool). Kept as the plan it was.

**Interim equation.** A fourth-order oscillator, L = ½ m q̇² − V(q) − ½ m τ² q̈² (Anselmi §4, eq. 4.1: L = (m/2)(v² − τ²a²) − V). Auxiliary Q with x = q + Q; Q = −τ²⟨q̈⟩ in the projected theory.

**Kernel.** G_F(u,τ) = sin(|u|/τ)/(2τ) = ½[G_ret(u) + G_adv(u)], G_F → δ(u) as τ → 0, with integration limits a → −∞, b → +∞ taken as a limit (2510.05276).

**Projected equation.** m q̈(t) = −∫ du G_F(u,τ) V′(q(t+u)) =: ⟨F⟩(t). Two initial data, not four.

**Iteration (as planned).** q = q₀ + q₁ + … with q₀ the τ = 0 solution and each q_n obtained by applying the kernel to the previous force; the convergence criterion the ratio of successive terms as a function of τω_eff. The question for a chaotic base flow: does intermittently large ω_eff make it fail, and how does the failure announce itself.

**What would change in minisuperspace (if the ħ → 0 check is ever taken to Bianchi IX).** On Bianchi IX the operator being inverted is the χ equation from ABP's auxiliary-field split, with C ≠ 0 and time-dependent coefficients from the Misner variables; the kernel must be constructed numerically. Two things from the hdym notes enter here and are not in any fakeon paper: (i) the Lee–Wick paper's open item (v) lesson that constraint terms vanishing on an ansatz couple the normal and ghost sectors off it applies to the Hamiltonian and momentum constraints, so the projection has to be defined on the constraint surface and checked to stay there, and the naive projection in the β± variables alone is to be labelled as such; (ii) the adiabatic/repeller boundary at H ≈ 0.2M is a prediction for where the kernel's reality condition m_χ² > H²/4 (section 3) is lost in the de Sitter approximation, H/m_χ = 0.5, versus where the unprojected ghost actually turns on, H/M ≈ 0.2.

**Agravity on the same backgrounds (same status).** Add R²/(6f_0²) to the `bianchi9/` action; φ is then a normal scalar in both theories (GSF), so the only difference between the fakeon and agravity runs is the treatment of χ. Parameters f_2 and f_0 with f_0/4 < f_2 < f_0 for the fakeon-consistent window, plus runs outside it.

## 8. Sources

All arXiv numbers were seen on an arXiv page during the sweep. Journal references come from the arXiv listings. Items marked † were confirmed from Anselmi's arXiv author listing only.

| arXiv | Authors, short title | Year | Journal |
|---|---|---|---|
| [1703.04584](https://arxiv.org/abs/1703.04584) | Anselmi & Piva, new formulation of LW QFT | 2017 | JHEP 06 (2017) 066 |
| [1703.05563](https://arxiv.org/abs/1703.05563) | Anselmi & Piva, perturbative unitarity of LW QFT | 2017 | PRD 96, 045009 |
| [1704.07728](https://arxiv.org/abs/1704.07728) | Anselmi, QFT of the gravitational interactions | 2017 | JHEP 06 (2017) 086 |
| [1801.00915](https://arxiv.org/abs/1801.00915) | Anselmi, Fakeons and Lee-Wick models | 2018 | JHEP 02 (2018) 141 |
| [1806.03605](https://arxiv.org/abs/1806.03605) | Anselmi & Piva, QG, fakeons and microcausality | 2018 | JHEP 11 (2018) 021 |
| [1809.05037](https://arxiv.org/abs/1809.05037) | Anselmi, fakeons, microcausality and the classical limit of QG | 2018 | CQG 36 (2019) 065010 |
| [1901.09273](https://arxiv.org/abs/1901.09273) | Anselmi, classicization of QG: the FLRW metric | 2019 | JHEP 04 (2019) 061 |
| [1909.12873](https://arxiv.org/abs/1909.12873) | Anselmi & Marino, light cones, GWs and H₀ | 2019 | CQG 37 (2020) 095003 |
| [1911.10343](https://arxiv.org/abs/1911.10343) | Anselmi, fakeons, QG and the correspondence principle | 2019 | proceedings |
| [2001.01942](https://arxiv.org/abs/2001.01942) | Anselmi, fakeons vs Feynman–Wheeler particles | 2020 | JHEP 03 (2020) 142 |
| [2005.10293](https://arxiv.org/abs/2005.10293) | Anselmi, Bianchi & Piva, inflation with the Weyl² term (ABP20) | 2020 | JHEP 07 (2020) 211 |
| [2007.15023](https://arxiv.org/abs/2007.15023)† | Anselmi, inflation as RG flow | 2020 | JCAP 01 (2021) 048 |
| [2010.04739](https://arxiv.org/abs/2010.04739)† | Anselmi, high-order corrections to inflationary spectra | 2020 | JCAP 02 (2021) 029 |
| [2103.01653](https://arxiv.org/abs/2103.01653)† | Anselmi, Fruzza & Piva, RG techniques for single-field inflation | 2021 | CQG 38 (2021) 225011 |
| [2105.05864](https://arxiv.org/abs/2105.05864) | Anselmi, double-field inflation and QG cosmology | 2021 | JCAP 07 (2021) 037 |
| [2202.10483](https://arxiv.org/abs/2202.10483) | Anselmi, purely virtual particles vs Lee-Wick ghosts | 2022 | PRD 105, 125017 |
| [2203.02516](https://arxiv.org/abs/2203.02516) | Anselmi, purely virtual particles review | 2022 | Symmetry 14, 521 |
| [2410.21599](https://arxiv.org/abs/2410.21599) | Anselmi, QG with PVPs from asymptotically local QFT | 2024 | — |
| [2503.01841](https://arxiv.org/abs/2503.01841) | Anselmi, Briscese, Calcagni & Modesto, amplitude prescriptions with complex poles | 2025 | JHEP 05 (2025) 145 |
| [2510.05276](https://arxiv.org/abs/2510.05276) | Anselmi & Calcagni, classicized dynamics and initial conditions | 2025 | JHEP 01 (2026) 104 |
| [2510.15157](https://arxiv.org/abs/2510.15157) | Calcagni, diffusion method in fakeon theories | 2025 | PRD 113, 025006 |
| [2601.06346](https://arxiv.org/abs/2601.06346) | Anselmi, on causality and predictivity | 2026 | Found. Phys. 56, 27 |
| [1403.4226](https://arxiv.org/abs/1403.4226) | Salvio & Strumia, Agravity | 2014 | JHEP 06 (2014) 080 |
| [1512.01237](https://arxiv.org/abs/1512.01237) | Salvio & Strumia, QM of 4-derivative theories | 2016 | EPJC 76, 227 |
| [1705.03896](https://arxiv.org/abs/1705.03896) | Salvio & Strumia, Agravity up to infinite energy | 2018 | EPJC 78, 124 |
| [1709.04925](https://arxiv.org/abs/1709.04925) | Strumia, QM with indefinite norm | 2017 | — |
| [1804.09944](https://arxiv.org/abs/1804.09944) | Salvio, Quadratic gravity (review) | 2018 | Front. Phys. 6, 77 |
| [1902.09557](https://arxiv.org/abs/1902.09557) | Salvio, Metastability in quadratic gravity | 2019 | PRD 99, 103507 |
| [2007.05541](https://arxiv.org/abs/2007.05541) | Gross, Strumia, Teresi & Zirilli, is negative kinetic energy metastable? | 2021 | PRD 103, 115025 |
| [2012.11608](https://arxiv.org/abs/2012.11608) | Salvio, dimensional transmutation in gravity and cosmology | 2021 | IJMPA 36 |
| [2404.08034](https://arxiv.org/abs/2404.08034) | Salvio, non-perturbative background-independent quadratic gravity | 2024 | JCAP 07 (2024) 092 |
| [2406.12958](https://arxiv.org/abs/2406.12958) | Salvio, unimodular quadratic gravity | 2024 | PLB 856, 138920 |
| [2507.08803](https://arxiv.org/abs/2507.08803) | Salvio, Strumia & Vitti, IR-enhanced loops in quadratic gravity | 2025 | JHEP, to appear |
| [1804.04980](https://arxiv.org/abs/1804.04980) | Donoghue & Menezes, gauge-assisted quadratic gravity | 2018 | PRD 97, 126005 |
| [1812.03603](https://arxiv.org/abs/1812.03603) | Donoghue & Menezes, massive poles in Lee-Wick QFT | 2019 | PRD 99, 065017 |
| [1908.02416](https://arxiv.org/abs/1908.02416) | Donoghue & Menezes, unitarity, stability and loops of unstable ghosts | 2019 | PRD 100, 105006 |
| [1908.04170](https://arxiv.org/abs/1908.04170) | Donoghue & Menezes, arrow of causality and QG | 2019 | PRL 123, 171601 |
| [2003.09047](https://arxiv.org/abs/2003.09047) | Donoghue & Menezes, quantum causality and the arrows of time | 2020 | PPNP |
| [2105.00898](https://arxiv.org/abs/2105.00898) | Donoghue & Menezes, Ostrogradsky instability overcome by quantum physics | 2021 | — |
| [2106.05912](https://arxiv.org/abs/2106.05912) | Donoghue & Menezes, causality and gravity | 2021 | JHEP 11 (2021) 010 |
| [2112.01974](https://arxiv.org/abs/2112.01974) | Donoghue & Menezes, On quadratic gravity | 2021 | special volume |
| [1909.01169](https://arxiv.org/abs/1909.01169) | Lehners & Stelle, a safe beginning for the universe? | 2019 | PRD 100 |
| [2312.14048](https://arxiv.org/abs/2312.14048) | Lehners & Stelle, higher-order gravity, finite action, safe beginning | 2023 | — |
| [1903.01044](https://arxiv.org/abs/1903.01044) | Castardelli dos Reis, Chapiro & Shapiro, Bianchi I beyond linear stability | 2019 | — |
| [2403.10377](https://arxiv.org/abs/2403.10377) | de Medeiros, Lazo, Müller & Sales, tilt in quadratic gravity | 2024 | — |
| [2502.03653](https://arxiv.org/abs/2502.03653) | de Medeiros et al., tilt in quadratic gravity II | 2025 | — |
| [1603.02851](https://arxiv.org/abs/1603.02851) | Toporensky & Müller, Bianchi I near the singularity in quadratic gravity | 2016 | — |
| [1007.4669](https://arxiv.org/abs/1007.4669) | Middleton, Kasner-like solutions in f(R, Ric², Riem²) | 2010 | — |
| [gr-qc/0610013](https://arxiv.org/abs/gr-qc/0610013) | Barrow & Hervik, Bianchi I and II in quadratic gravity | 2006 | — |
| [gr-qc/0601118](https://arxiv.org/abs/gr-qc/0601118) | Clifton & Barrow, modified Kasner solutions | 2006 | — |
| [gr-qc/9510065](https://arxiv.org/abs/gr-qc/9510065) | Demaret & Querella, Boulware Hamiltonian for R + βR² Bianchi I, IX | 1995 | CQG |
| [gr-qc/9811086](https://arxiv.org/abs/gr-qc/9811086) | Demaret, Querella & Scheen, Bach gravity on Bianchi I | 1998 | — |
| [gr-qc/9902044](https://arxiv.org/abs/gr-qc/9902044) | Querella, thesis | 1999 | — |
| [gr-qc/0011021](https://arxiv.org/abs/gr-qc/0011021) | Cotsakis, Flessas, Leach & Querella, pure R² Bianchi IX | 2000 | — |
| [gr-qc/0602017](https://arxiv.org/abs/gr-qc/0602017) | Schmidt, fourth-order gravity review | 2006 | — |
| [1512.05305](https://arxiv.org/abs/1512.05305) | Holdom & Ren, QCD analogy for QG | 2016 | PRD 93, 124030 |
| [0706.0207](https://arxiv.org/abs/0706.0207) | Bender & Mannheim, no-ghost theorem for Pais–Uhlenbeck | 2008 | PRL 100 |
| [1801.03220](https://arxiv.org/abs/1801.03220) | Mannheim, loop unitarity of the fourth-order propagator | 2018 | — |
| [2410.08278](https://arxiv.org/abs/2410.08278) | Kuntz, PT-symmetric quadratic gravity | 2024 | — |
| [2603.07150](https://arxiv.org/abs/2603.07150) | Kumar & Marto, inverted-oscillator ghost | 2026 | — |
| [2108.06294](https://arxiv.org/abs/2108.06294) | Deffayet, Mukohyama & Vikman, ghosts without runaway | 2021 | — |
| [2305.09631](https://arxiv.org/abs/2305.09631) | Deffayet, Held, Mukohyama & Vikman, global stability of ghost systems | 2023 | — |
| [2301.11317](https://arxiv.org/abs/2301.11317) | Fring & Turner, benign ghosts from Toda lattices | 2023 | — |
| [2306.13714](https://arxiv.org/abs/2306.13714) | Hell, Lüst & Zoupanos, conformal-gravity ghost by boundary conditions | 2023 | — |
| [2504.15005](https://arxiv.org/abs/2504.15005) | Sutton, de Felice & Sakellariadou, constraints on quadratic Weyl models | 2025 | — |

Non-arXiv: Barrow & Sirousse-Zia, PRD 39, 2187 (1989); Barrow & Cotsakis, PLB 232 (1989); Cotsakis, Demaret, De Rop & Querella, PRD 48, 4595 (1993); Berkin, PRD 42, 1016 (1990); Caprasse, Demaret, Gatermann & Melenk, IJMPC 2, 601 (1991). Johansson & Nohle 1707.02965 and the other items already in the hdym notes were not re-checked.
