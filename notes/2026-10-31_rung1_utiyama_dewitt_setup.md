# Rung 1 setup: Utiyama & DeWitt 1962, reconstructed from the abstract (31 Oct 2026)

Ladder: `notes/2026-10-31_chronology_1939_2014.md` §5. This note is the §2.3 reading before any building:
the papers summarised, the computation pinned down so two people can do it independently, and what counts
as agreement. Nothing is computed here.

## 1. Access, and how the rung is scored

Utiyama & DeWitt, J. Math. Phys. 3, 608–618 (1962), is paywalled and we are not buying or borrowing it (the
author's call, 31 Oct). The rung is therefore a **reconstruction from the abstract**: the setup, regulator,
and claims are fixed by the abstract; the paper's own coefficients are unseen. The 1965 DeWitt book is
likewise out of reach; the proper-time technique is taken from Avramidi hep-th/9510140 and Vassilevich
hep-th/0306138 instead.

**Checks, earliest first, none predating the paper** (author's rule, 31 Oct):

1. **'t Hooft & Veltman, Ann. Inst. H. Poincaré A 20, 69 (1974), §3 — read 31 Oct from the open numdam PDF.**
   Conventions: Euclidean ("purely imaginary time") metric; 1/ε ≡ 1/(8π²(n − 4)) (eq. 3.3); R defined in
   their Appendix A **with the opposite sign to ours** (their eq. 3.18 is negative on spheres) — added 31 Oct
   after review; the comparison in the results note uses this flip. The one-loop algorithm (eqs. 3.1–3.3): for L = −½∂φ_i∂φ_i + φ_i N^μ_ij ∂_μφ_j + ½φ_i M_ij φ_j,
   ΔL = (1/ε) Tr(¼X² + (1/24)Y_μνY_μν), X = M − N^μN^μ, Y_μν = ∂_μN^ν − ∂_νN^μ + N^μN^ν − N^νN^μ. Covariant
   version for a real scalar, L = √g(−½∂_μφ g^μν∂_νφ + φN^μ∂_μφ + ½φMφ) (eq. 3.39):

       ΔL = (√g/ε) Tr{ (1/24) Y^μνY_μν + ¼ (M − N^μN_μ − D_μN^μ − R/6)² + (1/120)(R_μνR^μν − R²/3) }   (3.40)

   introduced by the sentence "Eq. (3.40) contains a well known result" — Utiyama–DeWitt's. Hence the
   **1974 targets for rung 1**: minimal scalar (N = M = 0): ΔL = (√g/ε)[(1/120) R_μνR^μν + (1/240) R²];
   with Rφ²/12 added (the improved tensor, conformal coupling): ΔL = (√g/ε)(1/120)(R_μνR^μν − R²/3), which is
   (√g/ε)(1/240) C² modulo Gauss–Bonnet (their eq. 3.21). Also on these pages: eq. (3.29), the pole part of the
   flat-space bubble with two s_μν = h_μν − ½δ_μν h_αα vertices, written out — the direct target of technique
   (c) below; and the remark that the tadpole ∫dⁿp p_μp_ν/p² vanishes in dimensional regularization — which is
   exactly where the 1962 Λ⁴ and Λ² terms live, so the technique difference (a) vs (c) is stated by the 1974
   authors. This is the primary "as written" check: twelve years after, same quantity. The same paper is
   rung 3, read once for both.

   *Paper consistency check, 31 Oct, before any code:* Vassilevich eq. 4.28 with E = −ξR, Ω = 0, W_div =
   a₄/(n − 4) and a₄'s (4π)^{−2} = 2/(8π²), converted to the 1/ε above, gives (1/120)Ric² + (1/240)R² at
   ξ = 0 and (1/240)C² at ξ = 1/6 (using Riem² = GB + 4Ric² − R²). Same numbers as (3.40). So techniques
   (b) and (c) agree with 1974 on paper; the code is the certified record, and (a) is the open computation.

2. Capper & Duff, Nuovo Cim. A 23, 173 (1974) — the scalar-loop graviton self-energy in dimensional
   regularization, the trace-anomaly coefficient 1/120 (paywalled; cited, not used). [m?]
3. **Salvio, 1804.09944 (2018), eqs. (3.15)–(3.16)** — the modern anchor, §2 below.

If the 1962 text is ever in hand, this note gets a one-line re-score.

## 2. Summaries

**Utiyama & DeWitt 1962, from the abstract only.** A set of quantized matter fields propagates on a classical
metric that is itself sourced by ⟨T_μν⟩ (semiclassical gravity, self-consistently). ⟨T_μν⟩ is divergent with
three orders of divergence in a momentum cutoff Λ: Λ⁴, Λ², log Λ. Each is removed by a counterterm in
Einstein's equation: the Λ⁴ term by a cosmological term, the Λ² term by a renormalization of G, the log Λ
term by a counterterm derivable from a Lagrangian quadratic in the Riemann tensor. Second result: the
gravitational Green's function of this semiclassical theory — the graviton propagator with the matter loop
resummed into it — falls as 1/p⁴ rather than 1/p², so is less singular at short distance than the bare
linearized one. (That 1/p⁴ is Stelle's propagator fifteen years early, and Tomboulis's 1/N mechanism; it is
the thing this rung carries forward to 1977.) Which matter fields: unknown from the abstract; "a set" suggests
scalar, spinor and photon. We do the scalar; the others are additions if wanted.

**Vassilevich, hep-th/0306138, Phys. Rep. 388, 279 (2003) — the heat-kernel coefficients.** For a Laplace-type
operator D = −(g^μν∇_μ∇_ν + E) on a manifold without boundary, Tr(f e^{−tD}) ~ Σ_k t^{(k−n)/2} a_k(f, D)
(eq. 2.21), with (eqs. 4.26–4.28)

    a₀ = (4π)^{−n/2} ∫ √g tr f
    a₂ = (4π)^{−n/2} (1/6) ∫ √g tr f (6E + R)
    a₄ = (4π)^{−n/2} (1/360) ∫ √g tr f (60 E_;kk + 60 RE + 180 E² + 12 R_;kk + 5R² − 2 Ric² + 2 Riem² + 30 Ω²)

Sign convention R = +2 on the unit S² (§2.1). For a real scalar with action ½∫√g(∂φ)² + ½ξRφ² + ½m²φ² the
operator has E = −ξR − m² (eq. 3.6 with U = m²φ²), Ω = 0; conformal coupling ξ = 1/6 in n = 4 (eq. 3.5).
The one-loop effective action is W = ½ ln det D (eq. 1.16); in n = 4 its logarithmic divergence is governed
by a₄ (eqs. 1.21, 2.27–2.31). With a proper-time cutoff the power divergences are governed by a₀ (quartic)
and a₂ (quadratic) — the three orders of the 1962 abstract, in proper-time form.

**Salvio, 1804.09944, Front. Phys. 6, 77 (2018), §3.2.1 — the modern statement of the matter contribution.**
Signature (+,−,−,−); L_gravity = R²/(6f₀²) − W²/(2f₂²) − G − (M̄_P²/2)R − Λ (eq. 2.8), G the Gauss–Bonnet
term. With τ ≡ ln(μ/μ₀)/(4π)² (eq. 3.21):

    df₂²/dτ = −f₂⁴ (133/10 + N_V/5 + N_F/20 + N_S/60)                                   (3.15)
    df₀²/dτ = (5/3) f₂⁴ + 5 f₂² f₀² + (5/6) f₀⁴ + (f₀⁴/12)(δ_ab + 6ξ_ab)(δ_ab + 6ξ_ab)    (3.16)

N_V, N_F, N_S the numbers of gauge fields, Weyl fermions and real scalars; ξ_ab the scalar non-minimal
couplings. The text says these terms are generated by matter loops "as originally showed in [3]" with [3] =
Utiyama–DeWitt. The 133/10, 5/3, 5, 5/6 are the graviton's own loop and are *not* part of rung 1 (no quantized
gravity in 1962). **The rung-1 targets are the matter pieces:** per real scalar, N_S/60 in (3.15) and
(1/12)(1 + 6ξ)² in (3.16).

**Avramidi, hep-th/9510140 (thesis, 1986; Nucl. Phys. B 355 (1991) 712 for the technique).** The Schwinger–
DeWitt / proper-time expansion written out as an algorithm, with the DeWitt iteration for the coefficients.
The source for step (b) below when Vassilevich's tables are not enough — i.e. when we want to *derive* a₄
rather than quote it.

## 3. The computation, pinned down

Theory: one real scalar, mass m, non-minimal coupling ξ, on a background metric. Flat-space perturbation
theory g_μν = η_μν + h_μν, Lorentzian signature for step (a) as in 1962, Euclidean for (b) and (c).

Quantity: the one-loop matter contribution to the graviton two-point function, Π_μνρσ(p) — equivalently the
part of ⟨T_μν⟩ linear in h — from the single scalar loop with two T_μν insertions. (The 1962 abstract speaks of
⟨T_μν⟩; its divergences at second order in h are exactly this self-energy, and the "1/p⁴ Green's function"
is the h-propagator with Π resummed.)

Technique sequence, each giving the same three numbers or explaining why not:

- **(a) 1962: momentum cutoff.** Feynman-parameter method P on the bubble, hard cutoff Λ on the loop momentum
  (or on the Feynman-parametrised denominator — the 1962 choice is unknown; both are run and the difference
  recorded). Expected: Λ⁴ η_μνη_ρσ-type term (cosmological), Λ² term with the structure of the linearized
  Einstein tensor (renormalizes G), log Λ term with the two transverse-traceless structures p⁴ P⁽²⁾ and
  p⁴ P⁽⁰⁾ (renormalize the C² and R² couplings). The log coefficients are the targets; the Λ², Λ⁴ coefficients
  are recorded (they are regulator-dependent and vanish in (c); that difference is a result of the rung, not a
  failure).
- **(b) 1965: proper time.** Same quantity from W = ½ ln det D with a proper-time cutoff; a₀, a₂, a₄ as above
  with E = −ξR − m². Three orders again, now covariant. Derived once via the DeWitt iteration (Avramidi), then
  checked against Vassilevich eq. 4.28.
- **(c) 1972–74: dimensional regularization.** Same bubble, method P, D = 4 − 2ε; only the log survives as a
  1/ε pole; the coefficient must equal (a)'s log coefficient and (b)'s a₄.

**Agreement means:** the log-Λ (or 1/ε, or a₄) coefficients of the C² and R² structures equal, first, the
coefficients of 't Hooft–Veltman eq. (3.40) — (1/120)Ric² + (1/240)R² at ξ = 0, (1/240)C² at ξ = 1/6, in
units of √g/(8π²(n − 4)) — and, converted to Salvio's
normalization, 1/60 and (1/12)(1 + 6ξ)² respectively (per real scalar; the sign convention of ξ between
the two papers to be fixed first, since Vassilevich's conformal value is +1/6 and Salvio's (1 + 6ξ) vanishes
at ξ = −1/6). Equivalently, in a₄ language, a real scalar contributes (1/(4π)²)[(1/120) W² + ½(ξ − 1/6)² R² +
Euler] — this is the number the generator computes and the checker certifies.

**Then the 1962 second result:** resum Π into the h-propagator and show the 1/p⁴ fall-off at p ≫ m. One line
once Π is in hand; recorded because it is the bridge to Stelle.

## 4. Methods

Method P throughout (a) and (c). Method M is not needed here; there is no threshold prescription to compare,
since nothing is a fakeon yet. (b) is not P or M; it is the technique of its day and is on the ladder as such.

## 5. The split (agreed 31 Oct)

On the code, not the math. Claude does the three computations; the author does the engineering around them
where it helps — language porting, memory, threading, and the hybrid move of simplifying a domain by hand,
handing the simplified problem to the computation, and plopping the result back into the full one (past
technique at the key points, modern machinery in between). No independent paper derivation; the second
method on each number is the next technique in the sequence, and the external check is §1.

## 6. Open before building

- ξ sign convention between Vassilevich and Salvio (fix in the first line of code and of the paper derivation).
- Cutoff placement in (a): on the loop momentum before or after Feynman parametrisation. Both run.
- Which matter fields 1962 included. Scalar only here; spinor and vector as later additions if wanted (the
  N_F/20 and N_V/5 of eq. 3.15 are the targets if so). Vector is abelian (Maxwell) — not Yang–Mills.
- Folder: none created yet. Code lands where the author says when there is code.
