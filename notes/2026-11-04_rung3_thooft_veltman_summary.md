# Rung 3, the paper summarised (4 Nov 2026): 't Hooft & Veltman, "One-loop divergencies in the theory of gravitation", Ann. Inst. H. Poincaré A 20, 69–94 (1974)

Rule 5: summarised to the author before anything is built. Read in full from the numdam scan (open; the PDF itself is not in git).
Status [v] throughout unless marked. This note is the paper; the computation plan is §2; the author's yes fixes it.

## 1. The paper

**All 26 pages read as images 3 Nov 19:45 (the author attached the numdam PDF); everything below is [v] unless marked.** Where memory
(the 19:35 draft) was wrong, it says so — the errors were all in the bookkeeping, which is where code goes wrong.

**Conventions (App. A).** Imaginary time, g_μν = δ_μν in flat space, κ = 1. **1/ε = 1/(8π²(n − 4))** (eq. 3.3) — not 4 − n. Dimensional
regularization throughout; tadpoles ∫d_n p p_μp_ν/p² = 0 (p. 78).

**The master formula (§3).** For real fields, eq. (3.39)–(3.40):

    L  = √g ( −½ ∂_μφ g^μν ∂_νφ + φ N^μ ∂_μφ + ½ φ M φ ),   N^μ antisymmetric, M symmetric in i, j,
    ΔL = (√g/ε) Tr { (1/24) Y^μν Y_μν + ¼ (M − N^μN_μ − R/6)² + (1/120)(R_μν R^μν − R²/3) · 𝟙 },
    Y_μν = D_μN_ν − D_νN_μ + N_μN_ν − N_νN_μ.

(The complex-field form (3.35) has 1/12, ½, 1/60 — memory had those coefficients on the real-field Lagrangian: wrong by 2.) Derived
by: flat-space result (3.2)–(3.5) with the "little theorem" that φ* → φ*Z leaves ΔL unchanged (3.6–3.7); the general form (3.9) with ten
coefficients; the conformally flat case g = δ(1 − f) fixes all but one (3.11–3.25) — and in that case √g(Ric² − R²/3) is itself a total
derivative (3.22), so the last coefficient needs the self-energy graph with two s_μν = h_μν − ½δh vertices (3.27–3.34). The Riem² term is
never needed (App. B: Riem² − 4Ric² + R² is a total derivative in four dimensions only, (B.6); the general-n recipe is given). Scalar in an
external field, minimal (M = N = 0): **(√g/ε)(R²/144 + Ric²/120 − R²/360) = (√g/ε)(R²/240 + Ric²/120)** (5.23) — the first smoke of rung
3 is this number from rung 1's dim-reg logs through (3.40). Adding (1/12)Rφ² kills the "unrenormalizable" MR and N^μN_μR terms; its
energy-momentum tensor is the improved one of Callan–Coleman–Jackiw (p. 79–80).

**Gravity in the loop (§4–5).** L = √g(−R − ½∂φ̄g∂φ̄), ḡ = g + h, φ̄ = φ̃ + φ (4.1–4.2); the second-order expansion (4.3)–(4.25) written out in
full, with the classical fields on shell so the linear part vanishes. Gauge fixing −½C_μ² with **C_α = g^{1/4}(h^ν_{μ,ν} − ½h^ν_{ν,μ} −
φ ∂_μφ̃) t^{μα}**, t the root of g^{μν} (4.27–4.28) — the scalar is in the gauge condition (memory had the pure-gravity gauge). The quadratic
Lagrangian is cast into the §3 form with X^{βμ}_{αν}, Y^α_β, Z (4.30–4.33), ten h-fields plus φ. **Doubling trick** (§5): add primed copies
to make complex fields, use (3.35), halve at the end; the replacement h* → h P⁻¹ with P^{αβμν} = ½g^{αμ}g^{βν} − ¼g^{αβ}g^{μν}
(5.2–5.4), not covariant but allowed since the h are now just φ_i. D_μD^μh_αβ ≠ (D_μD^μ)h as scalars: the connection on the h-bundle is
N^{μγ}_{αβ} = (−2g^{μπ}Γ^γ_{πα}δ^ν_β)_symm (5.6–5.8) and 𝒯 = D_μN^μ − N_μN^μ cancels in M (5.10–5.11). Pieces: Tr Y² = 6R² − 24Ric²
(5.14), Tr(P⁻¹X̄P⁻¹X̄) = 2R² + 6Ric² + 3((∂φ̃)²)² (5.15), Tr(−⅓P⁻¹X̄R + R²/36) = −31/18 R² (5.16), (5.17), (5.18); **the unit-matrix term of
(3.35) counts 11 times** (p. 87). Ghost Lagrangian **√g η*_μ{η^{μ,α}_{,α} − R_{αμ}η^α − (∂_αφ̃ ∂_μφ̃)η^α}** (5.19), Tr Y²_ghost = R² − 4Ric²
(5.20), ghost contribution **−(√g/ε){(1/6)R(∂φ̃)² + (17/60)R² + (7/30)Ric² + R^{αβ}∂_αφ̃∂_βφ̃ + ½((∂φ̃)²)²}** (5.21), minus sign for the loop.
Total (5.22): (√g/ε){(9/720)R² + (43/120)Ric² + ½((∂φ̃)²)² − (1/12)R(∂φ̃)² + 2(D_μD^μφ̃)²}. Subtracting the scalar loop (5.23) and setting
φ̃ = 0: **pure gravity ΔL = (√g/ε)(R²/120 + (7/20)Ric²)** (5.24).

**Equations of motion (§6).** D_μD^μφ̃ = 0, R_μν = −½ D_μφ̃ D_νφ̃, R = −½(Dφ̃)² (6.3) — *not* R_μν = 0 when the scalar is there (memory
said "on shell" loosely); pure gravity R_μν = 0 (6.4). Then **ΔL_grav,scal = (√g/ε)(203/80)R²** (6.5), **ΔL_grav = 0** (6.6): pure gravity's
divergences "can be transformed away by a field renormalization"; gravity + scalar's cannot — "of the non-renormalizable type".

**§7, the improved tensor.** L = √g(−R − ½∂φg∂φ + aRφ² + bR^{μν}∂φ∂φ); b cannot help by dimension; a = 1/12 helps diagrams without
internal gravitons only. Tool: the Weyl transformation g → gf (7.2–7.7) turns √g(−fR − ½(∂φ)² − ½m²φ²), f = 1 − aφ², into the class
(7.10) with power series in φ̃; done by Schoonschip; "a to 1/12 is of little help, the final result seems not to be of any simple form". The
paper ends: "A certain exhaustion however prevents us from further investigation, for the time being."

**Field choice and the measure** (p. 70–71, [v] from the scans): other basic fields give the same S-matrix, "the Jacobian of the
transformation being one (provided dimensional regularization is applied)" — the same sentence that hid rung 1's measure effect
(δ(0) = 0). Rung 3's ledger rows say so.

**The implemented-not-fitted check** [m]: Deser & van Nieuwenhuizen, PRD 10, 401 (1974), Einstein–Maxwell, (√g/ε)(137/60) Ric² on
shell. PRD not open; the number is textbook-stable and the computation is ours regardless.

**What rung 3 computes, and how** (draft, for the author's yes):
1. **(3.40) as code**: given the bundle, N_μ, M on a generic background, return ΔL's coefficients; smoked on rung 1's scalar → (5.23);
   then a vector with its ghosts (Maxwell on a background; also the `matter_loop/` standing item).
2. **The h-operator**: √g(−R − ½(∂φ̄)²) to second order about (g, φ̃) by sympy (rung 1's second-order machinery), the gauge (4.27), the
   P⁻¹ replacement, the connection (5.8), the ghost (5.19); the pieces (5.14)–(5.18), (5.20)–(5.21) reproduced one by one, then (5.22),
   (5.24). The Gauss–Bonnet identity checked on random algebraic curvature tensors (rung 1's fit).
3. **On shell**: (6.3) → 203/80; (6.4) → 0. Einstein–Maxwell → 137/60.
4. **Later techniques**: DeWitt's recursion (1965) on the tensor bundle with endomorphism as method P against (3.40) as method M; the
   modern a₂ tables (Vassilevich 2003) as the check; the 1985 generalized Schwinger–DeWitt only when the chronology reaches it.
5. **Predictions under rule 25**: the four numbers are now [v] from the paper, so they are checks, not predictions; the predictions
   are about *our* machinery — which pieces of (5.14)–(5.21) rung 1's expansion reproduces without new code, and whether the
   δ(0) = 0 field-choice statement hides a measure term a cutoff would show (rung 1's finding carried to rung 3). Flat-space-first
   is the paper's own order (3.2 → 3.25 → 3.34 → 3.40); ours follows it.
Size: large (chronology §5); several sessions; rule 16 per run.


## 2. The computation plan (draft 4 Nov, for the author's yes)

The paper's own order is flat space → conformally flat → the self-energy graph → (3.40) → gravitons; ours follows it.

| step | what | method M / method P (rule 3) | check |
|---|---|---|---|
| 1 | **(3.40) as code**: given the bundle, N_μ, M on a generic background, return ΔL's three coefficients | M: the formula. P: DeWitt's recursion (1965) on the bundle with endomorphism — rung 1's normal-coordinate collapse extended | the minimal scalar → (5.23) from rung 1's dim-reg logs; the conformal scalar; a vector with its ghosts (Maxwell on a background; the `matter_loop/` standing item) |
| 2 | **the h-operator**: √g(−R − ½(∂φ̄)²) to second order about (g, φ̃) by sympy (rung 1's second-order machinery), the gauge (4.27), −½C², the P⁻¹ replacement (5.3–5.4), the connection (5.8), the ghost (5.19) | M: (3.40) applied. P: the pieces (5.14)–(5.18), (5.20)–(5.21) reproduced one by one before the sum | (5.22); φ̃ = 0 → (5.24); Gauss–Bonnet (B.6) on random algebraic curvature tensors (rung 1's fit) |
| 3a | **on shell**: (6.3) → 203/80; (6.4) → 0 | | the paper's (6.5), (6.6) |
| 3b | **§7, the improved tensor** — rule-25 prediction (from the author's overnight s1, 3–4 Nov, hdym; to be re-derived here): with a ≠ 0 and no R² term, aRφ² is removed by the Weyl map (7.7) plus the field redefinition φ̃ = √6 artanh(φ/√6) at a = 1/12 (K = 1/F², F = 1 − φ²/6); the S-matrix is invariant under field redefinitions (p. 70), so **the on-shell divergence with a = 1/12 is (6.5) rewritten in the new field: (√g/ε)(203/80)·¼K²((∂φ)²)²**, which is "not of any simple form" — the paper computed off shell and stopped one step short of this sentence. Off shell the power series in φ̃ is not reproduced (Schoonschip; "unprintable"); on shell it is one line | M: the Weyl map in rung 1's expansion. P: (3.40) with a ≠ 0 directly, on shell | the two agree; if they do not, the prediction failed (rule 11) |
| 3c | **Einstein–Maxwell** (Deser–van Nieuwenhuizen 1974): (√g/ε)(137/60)Ric² on shell [m] | the same code, different matter | the implemented-not-fitted check |
| 4 | later techniques at their date: the modern a₂ tables (Vassilevich 2003) as the check on step 1; the 1985 generalized Schwinger–DeWitt only when the chronology reaches it (fourth-order operators are Stelle's) | | |

**Predictions under rule 25** (the numbers are [v] from the paper, so they are checks; the predictions are about our machinery):
(i) rung 1's second-order expansion reproduces (4.22)–(4.25) and (5.6)–(5.8) without new tensor code; (ii) the "Jacobian one
(provided dimensional regularization is applied)" statement (p. 70) hides a measure term that a cutoff would show — rung 1's finding
carried to rung 3; what it is, in the h-field, is the question; (iii) 3b above.

**Tools** (author, 3 Nov): sympy → generated C for the traces (hdym's codegen pattern) or FORM for the tensor algebra; the 22 cores
for the random-tensor fits; CUDA only if something becomes a sweep. Smokes sized for the sandbox, full runs for the author's machine,
each header says which (rule 15).

**Size**: large (chronology §5); several sessions; rule 16 per run. Step 1 first, this session.

## 3. Step 1, closed 4 Nov (author's runs 916253a, 739b1af) — rule 12 satisfied for scalar, vector, Maxwell

| | method M, `hv1974_master.py` ((3.40) pointwise) | method P, `rung3_dewitt_bundle.py` (DeWitt's recursion on the bundle, Tr a₂) | agreement |
|---|---|---|---|
| scalar(ξ), pre-GB | (ξ²/4 − ξ/12 + 1/240) R² + Ric²/120 | (6ξ − 1)²/72 R² − Ric²/180 + Riem²/180 (= rung 1's [a₂] with E = −ξR) | ½P − M = (1/360)(1, −4, 1)·(R², Ric², Riem²): the Gauss–Bonnet density, as predicted |
| scalar(ξ), after GB | same | (ξ²/4 − ξ/12 + 1/240) R² + Ric²/120 | exact |
| vector (Feynman gauge), pre-GB | −R²/15 + (17/60)Ric² − Riem²/24 | −R²/9 + (43/90)Ric² − (11/180)Riem² (Gilkey's a₂, from memory, now [grid] on four backgrounds) | ½P − M = (1/90)(1, −4, 1) |
| vector, after GB | −R²/40 + (7/60)Ric² | same | exact |
| Maxwell = vector − 2 ghosts | −R²/30 + Ric²/10 = C²/20 | (1/180)(−25 R² + 88 Ric² − 13 Riem²); after GB /2: the same | exact; conformal, 12 × the conformal scalar |

Four random algebraic curvature tensors (seed 20261104), exact rationals; every input check PASS; on each background Tr a₂ equals the
predicted value to the last digit of 7–9-digit rationals (`results/rung3_dewitt_bundle.log`, 2118 s single-threaded on the author's
machine, 4 Nov). The two methods share no code beyond the random tensors. Also found: the parallel propagator's x² part is
+(1/6)R^a_{μbν}x^μx^ν (the header's hand sign was wrong; the run decided — rule 11) and a₁ = R/6 − ξR, −R/3 for the two bundles.
What these backgrounds cannot see (review of 4 Nov, in the header): the two total-derivative coefficients (□E, □R: zero by the
ansatz's x⁴ term) and, for the vector, the Ω² coefficient separately from E² and RE (confirmed in combination only).

Lesson, 4 Nov evening, not physics: `--jobs 4` on Windows re-imports an unguarded script, forking workers that each truncated the log
until the machine was restarted; the author lost an evening of compute to it. Fixed 8 Nov (main guard, non-empty logs refused, a line
per finished background).

## 4. For the ledger when the rows are written

The field-choice / δ(0) = 0 sentence (rung 1 ⇄ rung 3); the doubling trick and the non-covariant P⁻¹ as a 1974 technique; the
Weyl transformation (7.2–7.7) as the tool §7 reaches for; the "11 times the unit term" bookkeeping; (3.22) — on a conformally flat
background Ric² − R²/3 is a total derivative, so one coefficient needs flat-space graphs.
