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

## 4. Step 2, pure gravity by method M, 8 Nov (the author's runs 4d59af7, identical to the sandbox)

**2a, their operator as I first transcribed it.** (4.31) at φ̃ = 0, transcribed (wrongly, see below) and contracted as
½h^α_β X^{βμ}_{αν} h^ν_μ, with the P⁻¹ replacement (5.3)–(5.4): Tr(P⁻¹X̄) = −10R where (5.16) needs 6R, and no sign or index-order
reading of its five terms reproduces (5.15) and (5.16), even modulo Gauss–Bonnet (`graviton_loop/hv1974_graviton_M.py --transcribed`,
`results/hv1974_graviton_M_v1_transcribed431.log`). Its Ω and ghost pieces are fine: Tr Ω² = −6Riem² (5.14) and the ghost's
−(17/60)R² − (7/30)Ric² (5.21) on every background.

**2b, the operator derived from the action** (`graviton_loop/graviton_operator_derive.py`; the machinery validated first on Maxwell in
Feynman gauge, where it gives −D² + Ric exactly). The Euler–Lagrange operator of √g R − ½C² at the origin of normal coordinates, by
numeric jet probes: minimal (no mixed-index second-derivative terms; the ∂²h coefficient is exactly P δ^{μν}), and its endomorphism is
exactly a combination of the six covariant structures,

    E h = −R h + ½R δ tr h + (Ric h + h Ric) − δ tr(Ric h) − Ric tr h + 2R_{ambn}h_{mn}      [Gilkey's E of −(D² + E); standard sign]

with Tr E = −6R (their +6R_HV: (5.16) exact), Tr E² = 2R² + 6Ric² + 3·GB ((5.15) modulo Gauss–Bonnet — their "apart from total
derivatives" after (5.14)), and, through (3.40) with the ghost, on six random backgrounds:

| | pre-GB (R², Ric², Riem²) | after Gauss–Bonnet | the paper |
|---|---|---|---|
| ghost (vector, E = +Ric, ×(−2)) | (−1/5, −17/30, 1/12) | (−17/60, −7/30) | (5.21) at φ̃ = 0 ✓ |
| h-loop alone | (19/24, −17/12, 1/2) — the 8 Nov version of this row read (1/8, 11/15, −1/2), a hand copy that was wrong and that made the 9 Nov method-P run print FAIL (§4b) | (7/24, 7/12) | not printed; = (5.24) − (5.21) |
| **pure gravity, h + ghost** | (71/120, −119/60, 7/12) | **(1/120, 7/20)** | **(5.24) ✓** |

(`results/hv1974_graviton_M_sandbox.log`, 8 Nov; the (3.40) traces of the derived E on the derivation's background differ from the
fixed predictions by (7/12)·GB and (1/2)·GB, the Gauss–Bonnet multiples the fit then removes.)

**The printed (4.31), and a retraction.** Between 15:10 and 18:30 on 8 Nov this note claimed that (4.31) as printed was not the
operator of the paper's §5. The fresh review derived E by hand independently (the same six coefficients), reconstructed which (4.31)
would match it, and the print re-read at 4× showed my transcription had two misplaced subscripts (the δ placements of the two R
terms; δ^β_α R^μ_ν in the second Ricci term). Correctly transcribed, (4.31) **is** the derived operator, exactly, on two seeds, with
E = P⁻¹X̄ when the paper's Ric and R are read as minus the standard ones and its R^{βμ}_{αν} as the standard Riemann — a convention
statement the matter-loop sign argument could not fix for the Riemann sector. The paper is consistent from operator to result; the
derivation stands as the independent confirmation, and the misreading as the record. (`hv1974_graviton_M.py --transcribed` must
print "E == derived E: True"; `--transcribed --misread` is the 15:10 run.)

**Records (rule 11):** the hand sign of the connection term in normal coordinates was wrong (+⅓(Ric h + h Ric), not −⅓; the code
computes it from Γ); the first derivation stored the off-diagonal Euler–Lagrange components unhalved — E then fitted no covariant
structure and the kinetic normalisation read 2 — caught by the Maxwell check, fixed; (5.15) was predicted exact and holds modulo GB;
the method-M script lacked the adjoint symmetrization of X̄ (Tr E blind to it, Tr E² not); and the claim against the print, above.
Four errors of mine in one afternoon, each caught by a check that was already in the rules: the validation on a known case, the
fresh review, the re-read of the source at the point of disagreement.

## 4b. Step 2 closed under rule 12, 9 Nov (the author's overnight run, 4f8680a: `matter_loop/results/rung3_dewitt_graviton.log`, 4680 s)

DeWitt's recursion on the 10-dimensional h-bundle (the induced connection, the derived E as a field) and on the ghost, four random
backgrounds, exact rationals:

| | method P, Tr a₂ pre-GB | ½Tr a₂ after GB | method M after GB | ½Tr a₂ − M, pre-GB |
|---|---|---|---|---|
| graviton | (59/36, −55/18, 19/18) | **(7/24, 7/12)** | (7/24, 7/12) | (1/36, −1/9, 1/36) = (10/360)(1, −4, 1) |
| ghost (E = +Ric, before the −2) | (2/9, 43/90, −11/180) | (17/120, 7/60) | (17/120, 7/60) | (1/90, −2/45, 1/90) = (4/360)(1, −4, 1) |
| **pure gravity = graviton − 2 ghost** | (43/36, −361/90, 53/45) | **(1/120, 7/20)** | (1/120, 7/20) | |

The two methods agree exactly after Gauss–Bonnet and differ before it by exactly the (d/360)·GB density, as for the scalar and the vector:
**'t Hooft–Veltman's (5.24) by two methods sharing no code beyond the random tensors. Step 2 is closed.**

The run's own status line says FAIL. That is a transcription error of mine, not a result: the script's prediction table held method
M's pre-GB graviton coefficients copied by hand from the M log as (1/8, 11/15, −1/2) where the log says (19/24, −17/12, 1/2) (the
false row propagated into §4's table above, now corrected), so the per-sample "predicted" values and the comparison line were wrong
while every computed number was right — visible in the log itself, whose "after GB" line equals method M's. Third hand-copied number
of the week to be wrong. Fix (patch 0013): method M's (3.40) is computed inside the recursion script on the same backgrounds and every
prediction derives from it; nothing is typed from one log into another any more. The FAIL line stays in the committed log (rule 11).

## 4c. Step 3a, gravity + scalar by method M, 9 Nov (sandbox; the author's run pending) — 203/80 reproduced, and a slip in (5.22)

`graviton_loop/hv1974_gravscalar_M.py`: the 11-dimensional bundle (h, φ) of (4.30) with X (4.31), Y (4.32), Z (4.33) as printed — all
re-read at 4× before typing (the first φ̃ term of (4.31) is −½δ^β_ν D^μφ̃D_αφ̃, first derivatives; at page resolution it reads as a
second derivative) — under the convention map of 8 Nov, the ghost (5.19) with E = Ric − ∂φ̃∂φ̃, (3.40), twelve random backgrounds
carrying a curvature tensor, v = ∂φ̃ and W = DDφ̃ independently (off shell), eight structures fitted, Gauss–Bonnet and the Ricci
identity W·W → (tr W)² − Ric·vv as the reductions (the paper's "apart from total derivatives").

| | in standard sign, after the reductions | the paper |
|---|---|---|
| ghost | −(17/60)R² − (7/30)Ric² + (1/6)R v² + Ric·vv − ½v⁴ | (5.21) with Ric_HV, R_HV = −standard ✓ |
| total, five coefficients | (1/80)R² + (43/120)Ric² + ½v⁴ + (1/12)R v² + 0·Ric·vv | (5.22) ✓ (its −(1/12) is +1/12 in standard sign) |
| total, (□φ̃)² | **½** | **(5.22) prints 2** — see below |
| **on shell** (6.3): Ric = ½vv, R = ½v², □φ̃ = 0 | **203/320 (∂φ̃)⁴ = (203/80)R²** | **(6.5) ✓** |
| v = W = 0 | (1/80)R² + (43/120)Ric² = (5.24) + (5.23) | ✓ |

**The (□φ̃)² coefficient.** The paper's own (5.17), 2Tr(ȲP⁻¹Ȳ) = 2(□φ̃)² + 4R^{αβ}∂φ̃∂φ̃, is what our mixed block gives after the
Ricci identity (2(tr W)² − 4Ric·vv in standard sign); through (3.35)'s ½ on Tr(𝓜 − R/6)² and the halving of the doubling it
contributes ½(□φ̃)² to the total, and its R^{αβ}∂φ̃∂φ̃ cancels the ghost's (5.21) term — which is why (5.22) has no such term, as ours
has none. The 2(D_μD^μφ̃)² printed in (5.22) is therefore inconsistent with the paper's own (5.17): a slip in the final line, with no
consequence — it vanishes on shell (□φ̃ = 0) and in pure gravity, so (5.24) and (6.5) never depended on it. The fresh review of
9 Nov (separate agent, from the page images) assembled (5.22) by hand from (5.14)–(5.18), (5.21) and (3.35) and got ½, and confirmed
by sympy that our X and Y reproduce the printed intermediates (5.15)–(5.17) exactly; the slip is in the final addition, where the
(3.35) factor ½ and the halving were both skipped for this one term. One caveat stays: the ½ rests on Y (4.32) as printed — only
(5.17) tests Y's normalization — so step 3b (deriving (4.32) from the action) is what closes the claim beyond "consistent with the
paper's own pieces".

**Records (rule 11):** P3's Ric·vv sign in the header was my conversion error (their −{… + R^{αβ}∂φ∂φ …} is +Ric_std·vv; the run had
+1 and the page at 4× confirms); P1 was written against the printed 2 and re-scored against the ½ derived from (5.17) after the run.
Single method until the recursion runs the (h, φ) bundle with φ̃(x) as a field (method P for 3a: next).

**Next:** method P for 3a; 3b (the improved tensor as a field redefinition: the rule-25 prediction stands); 3c (Einstein–Maxwell, 137/60).

## 4d. Step 3a, method P, 9 Nov evening (sandbox; the author's run pending) — the recursion on the (h, φ) bundle with φ̃(x) a field

`matter_loop/rung3_dewitt_bundle.py --bundles gravscalar,ghostscalar --backend flint`: DeWitt's recursion on the 11-dimensional bundle
(h_ab, φ) and on the ghost with its ∂φ̃∂φ̃ term, on backgrounds that carry, besides the random curvature tensor, a scalar field
φ̃(x) with random jets through x⁴ (a separate random stream, seed + 1, so the curvature tensors are those of steps 1–2). E as a field,
the covariant form of (4.31)–(4.33) after the P⁻¹ replacement (v_a = ∂_aφ̃, W_ab = D_aD_bφ̃; all contractions with g⁻¹(x)):

    (E h)_ab = [the graviton's, §4] − (v_a (hv)_b + (hv)_a v_b) + ½(g_ab vhv + v_a v_b tr h) − ¼ g_ab v² tr h + ½ v² h_ab − 2 W_ab φ,
    (E φ)   = ½ tr W tr h − W^ab h_ab − v² φ;        ghost: E^a_b = Ric^a_b − v^a v_b.

Predictions written in the header before the first run (rule 25), and what the run gave on four backgrounds (26 s in all):

| | prediction | result |
|---|---|---|
| Q1 | E(0) and Tr Ω² equal, as matrices, what `hv1974_gravscalar_M.py`'s `build_E` makes from the printed (4.31)–(4.33) at the same (Riemann, v, W(0)) — imported and evaluated in place | equal on every background, both bundles |
| Q2 | Tr a₂(0) = 2·(3.40) + (d/180)·GB + (1/6)·□(Tr E)(0) | diff 0 on every background, both bundles; without the □E term the diff is exactly −(1/6)□(Tr E) |
| Q4 | Tr a₁(0) = (d/6)R + Tr E(0) | holds |
| total | ½Tr a₂[(h, φ) − 2 ghost] − method M's (3.40) total = (3/360)·GB + (1/12)·□(Tr E_h − 2 Tr E_gh) | holds |

**The □E term.** This is the first run where E is not built from curvature alone, so the first that sees Gilkey's total-derivative
term (zero in steps 1–2 by the ansatz's □R(0) = 0). Its coefficient was checked two ways before the prediction was written, because
memory was not to be trusted (author, 9 Nov): in the literature Gilkey's Theorem 4.18.16 as quoted in arXiv 1411.7851 (Prop. 14.12's
proof) has +60 tr(∇²E)/360, while arXiv 1507.06308 eq. (85), in the convention U = −E, prints the opposite sign for both box terms
(+(1/6)□U, −(1/30)□R) — the two agree on every magnitude and disagree on that sign; and in code, the recursion of this script, unmodified,
on a one-component bundle with symbolic E(x) = e₀ + e_a x^a + ½e_ab x^a x^b gives a₂(0) = ½e₀² + (1/6)Σe_aa in flat space and the same
+(1/6)□E(0) on a curved background (scratch, 54 s; not committed). +1/6, Gilkey's sign, for −(D² + E).

**The flint backend.** `--backend flint` replaces sympy's polynomial arithmetic by python-flint's `fmpq_mpoly` (same exact rationals;
the recursion, bundles and fits are the same lines). Validated by reproducing every computed number of the two committed sympy logs
(`rung3_dewitt_bundle.log`: 2118 s → 4 s; `rung3_dewitt_graviton.log`: 4680 s → 20 s, sandbox), and the new ghost bundle agrees between
the two backends (sympy 14 min, flint 1 s). sympy stays the default and the reference.

**What this establishes and what it does not (fresh review, 9 Nov evening; separate agent, with mutation tests).** My first draft of
this paragraph and of the script header said Q1 + Q2 close 3a under rule 12. Withdrawn: the review's must-fix.
- *Established:* on an 11-dimensional bundle with mixing blocks and an E not built from curvature, the recursion gives
  2·(3.40) + (d/180)·GB + (1/6)·□Tr E. Two methods for the **map** E → divergence.
- *Not tested: E itself.* Q2 holds for any matrix field E(x) — both sides are built from the same field (a mutated coefficient still
  gives diff 0) — and Q1 compares E(0) with method M's own input, the print. The φ̃ terms of X (4.31), Y (4.32), Z (4.33) and the
  ghost's ∂φ̃∂φ̃ term are single-source until **3b derives them from the action; 3a's numbers ((5.21), (5.22) with ½, 203/80) meet
  rule 12 then, not now.** (Step 2 was different: there E came from `graviton_operator_derive.py`, independently of the print.)
- *Not tested: the x-dependence of E(x).* Tr a₂(0) feels it only through □(Tr E)(0), taken from the same field on both sides, and the
  off-diagonal blocks' x-dependence drops out of the trace. Mutations (g(x) → δ, v^a → v_a, extra x-dependent factors) left status=PASS.
  The physical, non-total-derivative part needs only E(0). A guard was added (Tr E(x) = −6R(x) − v²(x); ghost R(x) − v²(x), to x²);
  the reviewer checked □(Tr E)(0) = −□(v²)(0) from the raw jets independently, and the +1/6 from the Mehler kernel (my code derivation
  used this same recursion, so it was not independent on its own).
- The quartic jets of φ̃ enter no printed number (they sit in W's x² part, off-diagonal); the cubic enter through □(v²) only.
- Review otherwise: no arithmetic error found; flint reproduces all computed values of both committed logs; spawn-mode pool and
  `--resume` behave; the covariant E agrees with `build_E` for symbolic v, W.

## 4e. Step 3b, part one, 9 Nov night (sandbox; the author's run pending) — the (h, φ) operator derived from the action

`graviton_loop/gravscalar_operator_derive.py`: step 2b's jet-probe derivation, rewritten on python-flint and extended to the 11
components (h_ab, φ), on backgrounds carrying a curvature tensor and a scalar field φ̃(x); the Lagrangian is √g(R − ½(∂φ̄)²) in
standard sign plus −½C² with the paper's gauge (4.27), C_μ = D^νh_μν − ½D_μh − c·φ∂_μφ̃. Step 2's script is untouched. Predictions
D0–D5 in the header before the first run (rule 25).

| | prediction | sandbox |
|---|---|---|
| D0 | φ̃ = 0 reproduces step 2's E, the scalar decouples, s = 1 | holds |
| D1 | c = 1: no first-derivative terms at the origin, second-derivative part diag(P, 1)δ^{μν}; c = 0: first-derivative terms | holds; with c = 0, 22 nonzero entries and the mixed column and E_φφ change |
| D2 | on shell, derived E = printed E (4.31)–(4.33), 11 × 11 | holds |
| D3 | mixed block (Eφ)_ab = −2W_ab, (Eh)_φ = ½ tr W tr h − W·h; E_φφ = −v² — Y (4.32) with the print's normalisation | holds |
| D4 | (a question) off shell too? | **yes**: the print is the off-shell operator; no field equations were used in (4.31)–(4.33) |

**What it settles.** (i) E for the (h, φ) loop is now two-source — the print and the action — and the fresh review derived the φ̃ terms
by hand as a third. (ii) **The (□φ̃)² coefficient of (5.22) is ½**: it depends only on the mixed block, whose two halves are both
derived, and on there being no mixed first-derivative terms (D1); the caveat of §4c is closed, and the printed 2 is a slip in the
paper's final addition (an inference about the paper, now supported from the action).

**What it does not.** The ghost operator is not derived here: E = Ric − ∂φ̃∂φ̃ (5.19) is still from the print in code (the reviewer's
hand variation of the gauge condition agrees with it). So the ghost row (5.21), the remaining coefficients of (5.22) and **203/80
wait for the ghost derived in code before 3a meets rule 12 in full.** Next small step, then: vary C_μ under a diffeomorphism with
the scalar's transformation, in the same machinery.

**Fresh review (9 Nov night; 5 runs, mutations on scratch copies).** No error found; the comparison is a real two-source test.
Caught by the run: c = −1; kinetic weight 1 instead of ½ (through s — it crashed instead of printing FAIL, fixed); φ̃'s quadratic jet
W instead of W/2. Neutral, correctly: the metric inverse to x² instead of x³. The pool passes under a forced spawn start method.
Folded in: no crash when the normalisation probes disagree; PASS needs D0 and at least one on-shell and one off-shell background;
exit code 1 on FAIL; v drawn without zero components. Limits kept in the header: the first-derivative statement is at the origin
(everywhere by covariance); the randomized derivative-part probes miss with probability ~1/500 and probe 3 of 44 diagonal
second-derivative pairs — one `--full` background closes both; E itself comes from the 11 constant probes only.

**The author's runs, 9 Nov night (committed 36bcd82).** Main run (2 on-shell + 2 off-shell): status=PASS, D0–D3 on every
background, D4 True; the `--c 0` control FAIL as intended (22 first-derivative entries, mixed column and E_φφ changed); `--full` on
one background: every check True — which closes, on that background, the two limits kept above for the randomized probes. Timing
note: the randomized mode does not parallelise (about 250 s per background on 16 workers, as on the sandbox's 2), `--full` took 348 s
for nine times the work; on the author's machine `--full` is the way to run it.

## 4f. Step 3b, part two, 10 Nov morning (sandbox; the author's run pending) — the ghost operator derived

`graviton_loop/ghost_operator_derive.py`, on the machinery of §4e. Predictions G0–G3 in the header before the first run.

| | prediction | sandbox |
|---|---|---|
| G0 | δh_ab = D_aη_b + D_bη_a, δφ = k·η^a∂_aφ̃ with k = 1 is the symmetry: the *ungauged* second-variation operator annihilates it on shell (η vanishing at the origin) | 0 of 32 η fields nonzero on shell; off shell 16 of 16 linear nonzero; k = 2: FAIL as intended |
| G1 | the ghost operator M η = C[δ_ηh, δ_ηφ] has no first-derivative terms, unit second-derivative part | holds, all 60 base fields, on and off shell |
| G2 | E_gh = Ric − c·k·∂φ̃∂φ̃, no W term | holds, on and off shell |
| G3 | c = 0: E_gh = Ric | holds |

**Where 3a stands now** (the fresh review's wording, adopted): both operators behind (5.21), (5.22) and 203/80 are two-source (the
print and the action, plus hand derivations by two reviewers), and the map E → divergence is two-method (§4d). **Shared by both
methods and hand-checked only:** the Faddeev–Popov −2 and the complex ghost (theory input, applied identically in both, as in step
2); the two reductions (Gauss–Bonnet; W·W → (tr W)² − Ric·vv), coded once in `hv1974_gravscalar_M.reduce`; the on-shell substitution
(one line there; redone by hand by the reviewer: 1/320 + 43/480 + ½ + 1/24 = 203/320). Whether that is "closed under rule 12" is
the author's call; a bare "meets rule 12" would overstate.

**Fresh review (10 Nov; mutations on scratch copies).** No error in code or physics. It supplied the reason G0 works (the
derivative of the Noether identity: S″[δ_ηΦ̃] = −S′·(δ_η)′, so with η(0) = 0 only the field equations at the origin enter) and
checked the off-shell output against that formula exactly. It found which η test what — linear: the field equations and k;
diagonal quadratic: k only; off-diagonal quadratic and cubic: nothing — hence the script's 16 + 16 and no random sampling.
Recorded against me: my first smoke printed G0 PASS having drawn no linear η on the on-shell background (fixed before any further
run); the times I typed in the header were guesses, not clock readings (corrected from the logs' file times); an unused import of
the print script in a script whose point is independence from it (removed).

**Record (rule 11): patch 0015 deleted §5 of this note.** My edit script for §4d cut the file at the paragraph it was replacing and
dropped everything after it — the four lines of §5 below. Found while writing this section; restored here verbatim from 459e826.

## 5. For the ledger when the rows are written

The field-choice / δ(0) = 0 sentence (rung 1 ⇄ rung 3); the doubling trick and the non-covariant P⁻¹ as a 1974 technique; the
Weyl transformation (7.2–7.7) as the tool §7 reaches for; the "11 times the unit term" bookkeeping; (3.22) — on a conformally flat
background Ric² − R²/3 is a total derivative, so one coefficient needs flat-space graphs.
