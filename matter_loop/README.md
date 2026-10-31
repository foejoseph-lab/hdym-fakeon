# matter_loop — one-loop matter on a classical background (rung 1 of the ladder)

Concept folder: everything that computes the one-loop effective action of quantized matter on a classical
metric lives here, whatever the year or technique — Utiyama–DeWitt 1962 (cutoff), Schwinger proper time,
DeWitt's 1965 recursion, 't Hooft–Veltman 1974 §3 (dimensional regularization), Salvio 2018's matter
terms. Techniques accumulate inside; no sibling folders per paper.

**Asks:** redo Utiyama & DeWitt, J. Math. Phys. 3, 608 (1962) — the one-loop divergences of a quantized
scalar on a classical gravitational background, and the three counterterms they force (cosmological, Newton,
quadratic-in-Riemann) — with the technique of 1962 (momentum cutoff), then Schwinger proper time (1951, on the
perturbative graphs), then 1972–74 (dimensional regularization), then DeWitt's 1965 covariant recursion on a generic
background, and compare. The paper itself is paywalled and unseen; the computation is
reconstructed from its abstract and scored against 't Hooft & Veltman 1974 eq. (3.40) and Salvio 2018
eqs. 3.15–3.16. Setup: `../notes/2026-10-31_rung1_utiyama_dewitt_setup.md`. Results and findings:
`../notes/2026-10-31_rung1_utiyama_dewitt_results.md`.

**Run** (sympy; 1–8 min each, single-threaded, exact rationals):

    python3 rung1_run.py --regulator cutoff     --measure plain                  --log results/rung1_cutoff_plain.log
    python3 rung1_run.py --regulator cutoff     --measure covariant              --log results/rung1_cutoff_covariant.log
    python3 rung1_run.py --regulator propertime --measure plain     --heatkernel --log results/rung1_propertime_plain.log
    python3 rung1_run.py --regulator propertime --measure covariant --heatkernel --log results/rung1_propertime_covariant.log
    python3 rung1_dimreg.py --measure plain      --log results/rung1_dimreg_plain.log
    python3 rung1_dimreg.py --measure covariant  --log results/rung1_dimreg_covariant.log
    python3 rung1_dewitt.py --samples 6 --jobs 6             --log results/rung1_dewitt_samples.log
    python3 rung1_dewitt.py --samples 6 --jobs 6 --quartic 0 --log results/rung1_dewitt_quartic0.log

`--measure plain` integrates over φ (operator K = √g·D); `--measure covariant` over χ = g^{1/4}φ (DeWitt's
covariant measure). `--heatkernel` derives a₀, a₂, a₄ at second order in h from Tr e^{−tD} and compares with
Vassilevich hep-th/0306138 eqs. 4.26–4.28.

`rung1_dewitt.py` is technique (b2): DeWitt's recursion for the heat-kernel coefficients, run in Riemann normal
coordinates on random generic backgrounds (exact rationals, no expansion in h; `--jobs` over the samples; the
`--quartic 0` run is a control, see the results note §7). Setup and predictions:
`../notes/2026-10-31_rung1_b2_dewitt_recursion_setup.md`.

**Files.** `ud1962_lib.py` — conventions (header), vertices, angular averages, the three regulators'
integrals, metric expansion to second order, the tadpole as an operator symbol, fitting. `rung1_run.py` —
cutoff and proper time. `rung1_dimreg.py` — dimensional regularization. `rung1_dewitt.py` — the covariant
recursion. `results/` — the six graph logs and the two recursion logs (and `results/superseded/`, the development scripts and logs, kept because one of
them records a prediction that failed — see the results note §4.1); each ends
`# END <name> elapsed=<s>s status=...` with status `OK`, `LOG-OK_POWER-ORDERS-NONCOVARIANT` (the expected
outcome for a hard cutoff, and for proper time with the plain measure), or `MISMATCH`.

**Outcome in one line.** The log-divergent (quadratic-curvature) coefficients come out identical six ways on
the graphs and equal 't Hooft–Veltman (3.40); the covariant recursion on generic backgrounds reproduces them
and adds the Riem² coefficient the graphs cannot see (Euler term −1/360); the ∞⁴ and ∞² assignments of the
abstract are reproduced only by proper time with the covariant measure — the Feynman-parametrised hard cutoff
does not give them, and the plain measure misses ∞⁴ by exactly ½(ln√g)⁽²⁾ (units Λ⁴/16π²).
