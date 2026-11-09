# graviton_loop — one-loop divergences with the graviton in the loop (rung 3 of the ladder, step 2 on)

Concept folder: everything that puts the graviton (and its ghost) in the loop on a classical background, whatever the year or
technique — 't Hooft–Veltman 1974 §4–7, Deser–van Nieuwenhuizen 1974, DeWitt's recursion on the tensor bundle, the modern a₂
tables, and later the Stelle rungs. Step 1 (the master formula and the recursion on matter bundles) lives in `../matter_loop/`.
Summary of the paper and the plan: `../notes/2026-11-04_rung3_thooft_veltman_summary.md`.

**Asks (step 2):** reproduce 't Hooft–Veltman's pure-gravity counterterm (5.24), (√g/ε)(R²/120 + (7/20)Ric²), with the technique of
the day — their master formula (3.40) on the h-bundle and the ghost — and with DeWitt's recursion as the second method; first from
their operator as printed (2a), then from the operator derived from the action (2b); then gravity + scalar (3a, 3b) and
Einstein–Maxwell (3c).

**Conventions.** Standard curvature sign throughout (the paper's R is minus the standard one; `../matter_loop/hv1974_master.py`).
The h-bundle is the 10-dimensional space of symmetric tensors in the basis E_aa, E_ab + E_ba (a < b); operators are 10 × 10 matrices
in that basis; the field-space metric P = ½gg − ¼gg (5.2) and its inverse (5.4) are operators on it.

**Run** (seconds; exact rationals):

    python3 hv1974_graviton_M.py --samples 6 --log results/hv1974_graviton_M.log
    python3 hv1974_gravscalar_M.py --samples 12 --log results/hv1974_gravscalar_M.log

**Run** (one background, ~15 min single-threaded; the derivation) and the 10 s Maxwell validation:

    python3 graviton_operator_derive.py --log results/graviton_operator_derive.log
    python3 maxwell_operator_check.py > results/maxwell_operator_check.log

| file | what |
|---|---|
| `hv1974_graviton_M.py` (8 Nov) | method M: (3.40) on the h-bundle and the ghost, pointwise on random algebraic curvature tensors, fitted and Gauss–Bonnet-reduced; the paper's §5 pieces (5.14), (5.15) mod GB, (5.16), (5.21) scored one by one, then (5.24). Default: the operator derived from the action (closed form in the header). `--transcribed`: P⁻¹X̄ from the printed (4.31), which must agree with the derived E; `--transcribed --misread`: the first, wrong transcription, kept as the record (`results/hv1974_graviton_M_v1_transcribed431.log`: P2, P3, P5, P6 FAIL) |
| `hv1974_gravscalar_M.py` (9 Nov; step 3a, method M) | gravity + a massless scalar: the 11-dimensional (h, φ) bundle of (4.30)–(4.33) as printed (re-read at 4×), the ghost with its ∂φ̃∂φ̃ term, (3.40) on twelve random (curvature, ∂φ̃, DDφ̃) backgrounds; eight structures fitted, Gauss–Bonnet and the Ricci identity as the reductions. (5.21) ✓, (5.22) ✓ in five coefficients, **(6.5) on shell: (203/80)R² ✓**; the (□φ̃)² coefficient is ½ where (5.22) prints 2 — the paper's own (5.17) gives ½ (rung 3 note §4c) | `hv1974_gravscalar_M.log` |
| `graviton_operator_derive.py` (8 Nov) | step 2b: the gauge-fixed graviton operator derived from √g R − ½C² (C of (4.27)) in Riemann normal coordinates as the Euler–Lagrange operator at the origin, by numeric jet probes on a sparse polynomial ring; minimality and the kinetic normalisation measured; E fitted to the six covariant structures; the paper's traces and (5.24) checked; the printed (4.31) compared. Header: predictions, the record of the first runs (the hand sign of the connection term; the off-diagonal weighting bug caught by the Maxwell check) |
| `maxwell_operator_check.py` | the validation of the derivation machinery on Maxwell in Feynman gauge: minimal, E = −Ric exactly, connection term +⅓Ric |
| `results/` | the author's logs once run: `hv1974_graviton_M.log`, `graviton_operator_derive.log`, `maxwell_operator_check.log`; the sandbox's first-transcription run is committed as the record, `hv1974_graviton_M_v1_transcribed431.log`; the unhalved first derivation as `graviton_operator_derive_sandbox_v1_unhalved.log`; other sandbox logs are not (rule 15) |

**Found, 8 Nov.** The operator derived from the action,

    E h = −R h + ½R δ tr h + (Ric h + h Ric) − δ tr(Ric h) − Ric tr h + 2R_{ambn}h_{mn}   (Gilkey's E of −(D² + E), standard sign),

reproduces (5.16) exactly (Tr E = −6R = +6R_HV), (5.15) modulo Gauss–Bonnet (Tr E² = 2R² + 6Ric² + 3·GB), the ghost (5.21), and
(5.24) exactly after Gauss–Bonnet on six backgrounds; the h-loop alone is (7/24)R² + (7/12)Ric², which the paper does not print.
The printed (4.31), correctly transcribed, **is** this operator, with the paper's Ric and R read as minus the standard ones and its
R^{βμ}_{αν} as the standard Riemann — a convention statement the matter loop could not fix for the Riemann sector. **Record (rule 11):**
between 15:10 and 18:30 on 8 Nov the notes claimed the printed (4.31) was not the operator of §5; the claim rested on my transcription,
which had two misplaced subscripts, found when the fresh review derived E by hand, reconstructed the (4.31) that would match, and the
print was re-read at 4× magnification. Withdrawn; the first-transcription run is kept as `_v1_` and `--transcribed --misread`.
Method P on these bundles (DeWitt's recursion, 10-dimensional) is next, for the author's machine.
