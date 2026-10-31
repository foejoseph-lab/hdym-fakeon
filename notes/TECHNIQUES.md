# Techniques ledger

The product of the ladder (README, "The rule of the project"): one entry per technique per rung — what it was
built for, what it sees, what it misses, where on the path it is expected to pay off again, and whether it is
standard, forgotten, or ours. Every rung's results note ends by adding here. Entries cite the log they come from.

## Rung 1 — one-loop matter on a classical background (`matter_loop/`)

| technique (year) | sees | misses | pays off again | status |
|---|---|---|---|---|
| hard momentum cutoff on Feynman-parametrised graphs (1962) | the three orders of divergence as three powers of Λ | covariance of the Λ⁴, Λ² terms (`rung1_cutoff_*.log`): the shifted-momentum cutoff is not a covariant regulator, with either measure; one placement run, the unshifted one owed | nowhere as a regulator; as a diagnostic of what a non-covariant regulator does to power divergences, at rung 3 | standard, superseded |
| Schwinger proper time on the graphs (1951/65) | all three orders covariantly — **only with DeWitt's measure χ = g^{1/4}φ**; with φ the Λ⁴ term is off by exactly ½(ln√g)⁽²⁾ (`rung1_propertime_{plain,covariant}.log`) | the Riem² coefficient (Gauss–Bonnet at O(h²)) | rung 3: the measure/gauge-fixing bookkeeping of 't Hooft–Veltman §4–5 is the same issue with more pieces | standard; the measure point is half-forgotten (dim reg hides it by setting δ(0) = 0 — information lost, not simplified) |
| dimensional regularization on the graphs (1972–74) | the log order, cleanly; the full ξ-polynomial of eq. 3.40 (`rung1_dimreg_*.log`) | the power divergences entirely, by construction | every later rung; it is the common technique from 1974 on | standard |
| DeWitt's covariant recursion in normal coordinates (1965) | the Riem² coefficient (Euler term −1/360) on a generic background, with no expansion in h (`rung1_dewitt_samples.log`) | nothing at this order; needs the x⁴ term of the normal-coordinate metric to make □R(0) = 0, else the fit degenerates (`rung1_dewitt_quartic0.log`) | rung 3 (fourth-order operators, where the 1985 index-free method is heavy and the coordinate collapse is light); any a₄-type object | standard in 1965, now mostly replaced by the 1985 method; **the normal-coordinate collapse (σ = x²/2, Δ = 1/√g, recursion → polynomial division) and the fit over random exact algebraic curvature tensors are the two reusable moves** |

Cross-technique: keep all tensor components and fit (projecting by hand would have hidden the non-covariant Λ⁴
structure); derive the second-order densities from the metric every time; bubble-only numbers are not checks.
