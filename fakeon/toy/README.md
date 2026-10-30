# fakeon/toy/ — the fakeon projection in one degree of freedom

The smallest system on which Anselmi's classical fakeon projection can be run and its convergence measured:
a particle with L = (m/2)(q̇² − τ² q̈²) − V(q), whose projected equation of motion is m q̈ = −V′(⟨q⟩) with
⟨q⟩ = (1 + τ² d²/dt²)⁻¹ q. The projection is done the way it is done on gravitational backgrounds — order by
order in τ², each order an ordinary second-order ODE — and the question asked is whether that expansion
converges on (a) a single-frequency base flow and (b) a bounce off a wall. Results and their reading:
`notes/2026-10-30_step1_toy.md`.

## Scripts

| script | what it does | outputs in `results/` |
|---|---|---|
| `fakeon_series.py --harmonic --kmax 6` | order-by-order reduction for V = mω²q²/2; checks the reduced frequency against Anselmi's eq. (4.13) (Catalan series) | `series_harmonic.log`, `series_harmonic.npz` |
| `fakeon_series.py --wall --kmax 7` | one bounce off V = V₀e^{(q−q_w)/w}; trajectory convergence in the order K for τv/w = 0.05 … 0.6 | `series_wall.log`, `series_wall.npz` |
| `fakeon_toy.py --check` | the nonlocal kernel itself: regularised Fourier multiplier → principal value, τ → 0 limit | `check.log` |
| `fakeon_toy.py --picard` | the naive fixed-initial-data iteration of the nonlocal equation, kept because it does not contract | `picard.log` |

**Spyder:**
```
import sys
%cd C:\GitHub\hdym-fakeon\fakeon\toy
!"{sys.executable}" -B -u fakeon_series.py --harmonic --kmax 6   > results\series_harmonic.log
!"{sys.executable}" -B -u fakeon_series.py --wall --kmax 7       > results\series_wall.log
!"{sys.executable}" -B -u fakeon_toy.py --check                  > results\check.log
!"{sys.executable}" -B -u fakeon_toy.py --picard                 > results\picard.log
```

Wall clock on one core: the K = 7 wall sweep about 2 min, the Picard diagnostic about 1 min, the rest seconds; each log ends with its
elapsed seconds. Dependencies: numpy, scipy, sympy. Nothing here reads from hdym.

## Method in one paragraph

With q̈ = f_{K−1}(q, q̇) known to O(τ^{2K−2}), the derivatives q^(2k) are computed along that flow, substituted
into ⟨q⟩ = Σ_k (−τ²)^k q^(2k), and V′(⟨q⟩) is Taylor-expanded and truncated at τ^{2K}; the result is f_K.
The wall exponential is carried as a polynomial generator E with dE/dq = E/w, so every f_K is a polynomial in
(q, q̇, E, τ²) and the reduction is fast when the wall width is substituted numerically. Each f_K is integrated
with DOP853 at rtol 1e-13 from the same two initial data; a runaway or stiff iterate is cut off by an evaluation
budget and recorded as a failure. The convergence observable on the wall is the bounce time delay,
dq_K = |q_K − q_{K−1}| at the end of the run.
