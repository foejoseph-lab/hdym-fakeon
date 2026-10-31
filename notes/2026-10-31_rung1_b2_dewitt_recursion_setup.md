# Rung 1, technique (b2): DeWitt's 1965 covariant recursion, set up before building (31 Oct 2026, evening)

Owed item of `2026-10-31_rung1_utiyama_dewitt_results.md` §6. This is the §2.3 summary: the technique in the
words of the open sources, the computation pinned down, the predictions written before any code runs.

## 1. Sources (open; DeWitt 1965 itself unseen, cited as the paper of record)

- **DeWitt, *Dynamical Theory of Groups and Fields* (1965)** — out of print, not bought (rule 10; the one
  copy found earlier was ~$600). Unseen. The recursion below is his; every formula is taken from the two
  sources that follow, which present it as his.
- **Avramidi, hep-th/9510140 (thesis), ch. 1–2** — the Schwinger–DeWitt ansatz and the recursion written as
  an algorithm, with the covariant Taylor expansions in normal coordinates of σ, Δ and the a_k.
- **Vassilevich, hep-th/0306138, §2.2 and §4.1** — the ansatz (eq. 2.21 region), the integrated coefficients
  a₀, a₂, a₄ (eqs. 4.26–4.28). The check, read only after the run (rule 2 of brief 13: the method first).

## 2. The technique

For a Laplace-type operator D = −(g^{μν}∇_μ∇_ν + E) on scalars (here E = −ξR − m²), the heat kernel is
written (DeWitt 1965)

    K(t; x, x′) = (4πt)^{−n/2} Δ^{1/2}(x, x′) exp(−σ(x, x′)/2t) Σ_{k≥0} t^k a_k(x, x′),

with σ half the squared geodesic distance, Δ the Van Vleck–Morette determinant, a₀(x, x) = 1. Inserting into
(∂_t + D)K = 0 and using σ^{;μ}σ_{;μ} = 2σ and Δ^{−1}(Δσ^{;μ})_{;μ} = n gives the recursion

    σ^{;μ} ∇_μ a_k + k a_k = Δ^{−1/2} (□ + E) (Δ^{1/2} a_{k−1}),     k ≥ 1,  with  σ^{;μ}∇_μ a₀ = 0.

The integrated coefficients of Vassilevich are a_{2k}(f, D) = (4π)^{−n/2} ∫ √g f [a_k], where [·] is the
coincidence limit x → x′. The one-loop divergence in n = 4 is governed by [a₂]; the power divergences by
[a₀] = 1 and [a₁].

**The collapse in Riemann normal coordinates about x′.** With x′ at the origin of normal coordinates,
σ = ½ δ_{μν} x^μ x^ν exactly, σ^{;μ} = x^μ exactly (g_{μν}x^ν = δ_{μν}x^ν in normal coordinates), and
Δ(x, 0) = 1/√g(x). So σ^{;μ}∇_μ is the Euler operator x·∂ on scalars, and the recursion is solved order by
order in x: if the right-hand side has homogeneous parts P_d of degree d, then a_k = Σ_d P_d/(d + k).
Nothing is covariant-looking in the intermediate steps; covariance is recovered because the metric fed in is
the normal-coordinate expansion of a generic metric. This is how the coincidence limits were computed in the
1960s; the modern index-free "covariant Taylor" machinery (Avramidi, Barvinsky–Vilkovisky 1985) is the
1985 technique, applied when the ladder reaches 1985, not here.

**Input: the metric in normal coordinates.** For a generic metric with ∇R = 0 at the origin, the standard
normal-coordinate expansion with its derivative terms dropped,

    g_{μν}(x) = δ_{μν} − ⅓ R_{μανβ} x^α x^β + (2/45) R_{αμβλ} R_{λγνδ} x^α x^β x^γ x^δ + O(x⁵),

R = +2 on the unit S² (our convention), all indices Euclidean at the origin. Dropping the derivative terms of
the general expansion (−⅙∇R x³ and −(1/20)∇∇R x⁴) chooses a particular metric; [a₂] is a local invariant of
dimension 4, c₁R² + c₂Ric² + c₃Riem² + c₄□R + c₅□E with universal c's, so the only question is whether □R at
the origin vanishes on the chosen metric. It is computed, not assumed (first draft of this note assumed it
away; corrected before the run, see below). Within this family □R(0) is itself quadratic in R — a combination
of R², Ric², Riem² — so it cannot be fitted as an independent structure; the curvature-squared coefficients
are unambiguous only on a metric with □R(0) = 0, and the standard x⁴ term is what delivers that. A control run
with the x⁴ term removed (`--quartic 0`) must show □R(0) ≠ 0 and a degenerate fit.

*Corrected before the first logged run.* The first draft asked the code to check that ∇_{(γ}∇_{δ)}R_{αβμν} —
symmetrised in the derivative pair only, every component — vanishes at the origin. That object cannot vanish
for generic curvature: the differentiated second Bianchi identity ties its values in different index
positions to the antisymmetrised parts, which are quadratic in R by the Ricci identity. (The object that the
x⁴ term of the normal-coordinate expansion actually involves, symmetrised over all four of α, β, γ, δ, does
vanish on the 2/45 metric — the independent review of 31 Oct checked all its components on one sample — and
that is why □R(0) = 0 there.) A debug run showed the pair-only check failing on generic samples and passing on
S⁴ and S²×S² (symmetric spaces), which is how the mistake surfaced. The check was removed and □R(0) added to
the output and the fit; the removed block is kept as `matter_loop/results/superseded/rung1_dewitt_sym_nabla2_check.py`
(reconstructed from the session; the debug run's output was not saved). A second slip, caught by the `--quartic 0` control
before anything was written up: the first draft of the code header claimed the fitted coefficients do not
depend on the x⁴ term; they do, through □R(0), as the previous paragraph now says.

**Checks the code makes on its own input before the recursion runs** (check, don't assume):

1. Γ^μ_{αβ} x^α x^β = 0 through O(x³) — the defining property of normal coordinates. (This does not test the
   2/45: that term satisfies the condition identically by the antisymmetry of R.)
2. The Riemann tensor computed from g(x) at the origin equals the R_{αβγδ} fed in, with the sign convention
   R = +2 on S²; its first derivative at the origin vanishes (trivially — the metric is even in x — but
   checked); the scalar curvature of g(x) at the origin equals the contraction of the input; □R at the origin
   and the number of nonzero x² coefficients of R(x) are reported (R(x) − R(0) has a traceless x² part on
   generic samples; only its trace, through □E, could enter [a₂]).
3. Δ^{1/2}(x) = 1 + (1/12) R_{αβ} x^α x^β + O(x⁴) (a standard expansion; a check of the 1/√g identification).
4. [a₁] = E + R/6 — the a₂-coefficient of the six logs, obtained here from the recursion at k = 1.

## 3. The computation

Sympy, exact rationals. The 20 components of R_{αβγδ} at the origin are either symbols (`--symbolic`; the
result is then an exact polynomial identity) or random rational algebraic curvature tensors (`--samples N`,
built as sums of Kulkarni–Nomizu squares of random symmetric matrices, which satisfy all symmetries and the
first Bianchi identity by construction; the fit of the three curvature-squared coefficients needs three
independent samples and uses more as a consistency check — this is the embarrassingly parallel index for
`--jobs`). E is a constant symbol; R(x) − R(0) has a traceless x² part under the ansatz, so □E = −ξ□R(0), which
is zero on the backgrounds used (reported in every log), and nothing else of E's variation can enter [a₂].

Steps: metric to O(x⁴) → inverse metric, √g, Δ^{1/2} = g^{−1/4} to O(x⁴) (series of tr ln) → scalar
Laplacian □ = g^{−1/2}∂_μ(g^{1/2} g^{μν}∂_ν) → a₁(x) to O(x²) by the recursion → [a₂] = ½ [Δ^{−1/2}(□ + E)(Δ^{1/2}a₁)]
at the origin → fit to c_{E²}E² + c_{RE} R E + c_{R²}R² + c_{Ric²}Ric² + c_{Riem²}Riem² + c_□ □R (in `--symbolic`
mode the fit is an exact identity test, residual required to be the zero polynomial). E is a constant symbol
(m² and ξR at the origin); the □E total derivative is not represented and not fitted.

## 4. Predictions, written before the run

From the structure of the recursion (not from the table): [a₂] = ½[a₁]² + (curvature-squared terms with no
E) + total derivatives. So c_{E²} = ½ and c_{RE} = 1/6 follow from [a₁] = E + R/6 once that is confirmed.
From the six logs of the results note §2 (the O(h²) projection, where Riem² = 4Ric² − R² + total derivative;
in the a₄ normalisation, which is twice the Γ-coefficient quoted in §2 — the proper-time logs print it as
`Ric^2/60 + ... R^2/120`):

    c_{Ric²} + 4 c_{Riem²} = 1/60 ,     c_{R²} − c_{Riem²} = 1/120 ,

two equations for three unknowns. The recursion must satisfy both and supply the third number. **I expect**
c_{Riem²} = 1/180, c_{Ric²} = −1/180, c_{R²} = 1/72, c_□ = 1/30 where □R is nonzero, which satisfy the two constraints
(−1/180 + 4/180 = 3/180 = 1/60 ✓; 1/72 − 1/180 = (5 − 2)/360 = 1/120 ✓) and are the textbook values
(Vassilevich 4.28: (5R² − 2Ric² + 2Riem²)/360). If the run gives different numbers that satisfy both
constraints, the ladder has a finding; if it gives numbers that violate a constraint, the code is wrong.
Named checks after the fit: unit S⁴ (R_{abcd} = δ_acδ_bd − δ_adδ_bc, R = 12): [a₂]_{E=0} = 29/15;
S²×S² (unit spheres, R = 4, Ric² = 4, Riem² = 8): [a₂]_{E=0} = 16/72 − 4/180 + 8/180 = 11/45.

## 5. What it adds to the rung

- The Riem² coefficient, which the flat-space graphs cannot see (Gauss–Bonnet), and with it the full
  (W², R², Euler) split of the scalar's a₄: W² = Riem² − 2Ric² + R²/3 and Euler = Riem² − 4Ric² + R² give
  W²/120 + R²·(ξ-dependent) − E₄/360 for the predicted numbers, E₄ = Riem² − 4Ric² + R² — the Euler
  coefficient −1/360 is the one genuinely new number. (This sentence first read "+ Euler/360"; the sign was
  mis-predicted, not a convention — E₄ so defined is positive on S⁴ — and the History records the slip.)
- A fourth technique on the log order of rung 1, this one not a flat-space graph at all.

## 6. Scoring

Agreement: the two constraints from the six logs hold exactly, and the three coefficients equal Vassilevich
4.28 read after the run. The 1965 book stays unseen; the recursion is attributed to it via Avramidi and
Vassilevich, who both do so.

## History

- 31 Oct 2026, evening: written before any code. Then, before the first logged run: the symmetrised-∇∇R check
  removed and □R(0) added (§2). After the `--quartic 0` control: the claim that the x⁴ term does not matter
  withdrawn (§2); the Euler-coefficient sentence in §5 corrected from +1/360 to −1/360 (a mis-prediction of
  the sign, recorded as such).
- 31 Oct 2026, night, after the independent review: two input checks the text claimed but the code did not
  make (∇Riem(0) = 0, "R(x) = R(0) + O(x³)") — the first implemented, the second was false as stated and is
  replaced by what the code reports; the removed-check paragraph made precise (pair-only symmetrisation) and
  the removed block archived; the a₄-vs-Γ normalisation of the constraints stated; `--symbolic` marked not run.
