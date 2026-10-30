"""fakeon_toy.py -- the one-mode fakeon kernel (Anselmi 1809.05037 sec. 4) and the Picard diagnostic.

Interim (unprojected) theory, eq. (4.1):   L = (m/2)(q'^2 - tau^2 q''^2) - V(q)
Fakeon kernel, eq. (4.5):                   G_F(u, tau) = sin(|u|/tau) / (2 tau)   (half retarded + half advanced)
Projected ("finalized") equation, (4.8):    m q'' = -V'(<q>),   <q>(t) = int du G_F(u, tau) q(t - u)

What this script does.
--check   The kernel <.> on a uniform grid, computed with the regularised multiplier
          (1/2tau)[(1/tau+nu)/(eps^2+(1/tau+nu)^2) + (1/tau-nu)/(eps^2+(1/tau-nu)^2)]  ->  PV 1/(1 - tau^2 nu^2)
          (e^{-eps|u|} G_F, then eps -> 0 by Richardson; this is the a -> -inf, b -> +inf limit of Anselmi &
          Calcagni 2510.05276 taken as a limit), checked against the analytic multiplier on cosines and against
          the tau -> 0 limit on a Gaussian.
--picard  The naive fixed-initial-data iteration q_{n+1}'' = -V'(<q_n>)/m on the harmonic oscillator, on a
          window of several periods.  It does not contract: each iterate adds polynomial-in-t terms (the
          Taylor-in-t expansion of cos(Omega t) - cos(w t), the Lindstedt problem), which on a window of length
          T grow like (w T)^2 per step and dominate the tau^2 corrections; the undamped kernel also carries
          grid-edge content inward.  Kept as a record that the nonlocal projected equation is not to be solved
          this way; the window dependence means the numbers are diagnostic, not a convergence measurement.  The
          order-by-order local reduction in fakeon_series.py is the working method, and the one used on
          gravitational backgrounds.

Note on the nonlocal equation itself.  For V = m w^2 q^2/2 the line equation m q'' = -m w^2 <q> has four exponential
solutions, Omega^2 = (1 -/+ sqrt(1 - 4 tau^2 w^2))/(2 tau^2): the physical pair (regular as tau -> 0) and the
"fakeon-minus" pair at Omega' ~ 1/tau.  Two initial data do not fix a solution of the nonlocal equation; the
projection is the restriction to the physical pair (Anselmi 4.12-4.13), which the perturbative reduction performs
automatically and a direct solve of the nonlocal equation does not.  At w tau -> 1/2 the two pairs collide and go
complex and the projected solution space is empty.

Runs (Spyder):
    import sys
    %cd C:\\GitHub\\hdym\\fakeon\\toy
    !"{sys.executable}" -B -u fakeon_toy.py --check    > results\\check.log
    !"{sys.executable}" -B -u fakeon_toy.py --picard   > results\\picard.log
Every log ends '# END <name> status=...'.
"""
import argparse
import os
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)


# ----------------------------------------------------------------------------------------------------------
# kernel
# ----------------------------------------------------------------------------------------------------------
def kernel_multiplier(nu, tau, eps):
    """Fourier transform of e^{-eps|u|} sin(|u|/tau)/(2 tau) at angular frequency nu (real, even in nu)."""
    a = 1.0 / tau + nu
    b = 1.0 / tau - nu
    return (a / (eps**2 + a**2) + b / (eps**2 + b**2)) / (2.0 * tau)


def project_grid(f, h, tau, eps):
    """<f> on a uniform grid with step h, via zero-padded FFT convolution with the regularised kernel.

    f is assumed to be known on the whole grid; the kernel is truncated at |u| = n*h where the grid ends,
    so the grid must extend several 1/eps beyond the region where <f> is wanted.
    """
    n = f.size
    npad = 2 * n
    # exact Fourier multiplier of the regularised kernel at the DFT frequencies; the periodic image of the
    # kernel sits at distance npad*h and is suppressed by e^{-eps npad h}
    nu = 2.0 * np.pi * np.fft.rfftfreq(npad, d=h)
    F = np.fft.rfft(np.concatenate([f, np.zeros(npad - n)]))
    out = np.fft.irfft(F * kernel_multiplier(nu, tau, eps), npad)[:n]
    return out


def project(f, h, tau, eps_list):
    """<f> extrapolated eps -> 0 by Richardson on the eps_list (geometric, ratio 2).  Away from nu = 1/tau the
    regularised multiplier differs from the PV one at O(eps^2), so the pairwise step is (4 P(e/2) - P(e))/3."""
    vals = [project_grid(f, h, tau, e) for e in eps_list]
    rich = [(4.0 * vals[i + 1] - vals[i]) / 3.0 for i in range(len(vals) - 1)]
    return rich[-1], vals, rich


# ----------------------------------------------------------------------------------------------------------
# base flow and iteration
# ----------------------------------------------------------------------------------------------------------
def integrate_both_ways(force, t0, q0, v0, t_grid, m=1.0, rtol=1e-10, atol=1e-12):
    """m q'' = force(t) with q(t0) = q0, q'(t0) = v0, on t_grid (which contains t0), Radau both ways."""

    def rhs(t, y):
        return [y[1], force(t) / m]

    tl = t_grid[t_grid <= t0]
    tr = t_grid[t_grid >= t0]
    q = np.empty_like(t_grid)
    if tl.size:
        sl = solve_ivp(rhs, (t0, tl[0]), [q0, v0], t_eval=tl[::-1], method="Radau", rtol=rtol, atol=atol)
        q[: tl.size] = sl.y[0][::-1]
    if tr.size:
        sr = solve_ivp(rhs, (t0, tr[-1]), [q0, v0], t_eval=tr, method="Radau", rtol=rtol, atol=atol)
        q[t_grid.size - tr.size :] = sr.y[0]
    return q


def iterate_projection(Vprime, t_grid, t0, q0, v0, tau, eps_list, m=1.0, nmax=30, tol=1e-9, log=print):
    """q_{n+1}'' = -V'(<q_n>)/m, fixed (q0, v0).  Returns the iterates, the contraction ratios and status."""
    h = t_grid[1] - t_grid[0]
    # n = 0: tau = 0 solution (the force evaluated at the unprojected q)
    def force_self(q_fun):
        return lambda t: -Vprime(q_fun(t))

    q_prev = None
    q_cur = None
    rhos = []
    # start from the two-derivative solution
    tl, tr = t_grid[t_grid <= t0], t_grid[t_grid >= t0]
    q_cur = np.empty_like(t_grid)
    sol = solve_ivp(lambda t, y: [y[1], -Vprime(y[0]) / m], (t0, tr[-1]), [q0, v0],
                    t_eval=tr, method="Radau", rtol=1e-10, atol=1e-12)
    q_cur[t_grid.size - tr.size:] = sol.y[0]
    if tl.size:
        solL = solve_ivp(lambda t, y: [y[1], -Vprime(y[0]) / m], (t0, tl[0]), [q0, v0],
                         t_eval=tl[::-1], method="Radau", rtol=1e-10, atol=1e-12)
        q_cur[: tl.size] = solL.y[0][::-1]
    iterates = [q_cur.copy()]
    status = "max_iter"
    # norm on the central half of the grid (edges are kernel-truncation polluted)
    n = t_grid.size
    core = slice(3 * n // 8, 5 * n // 8)
    for it in range(1, nmax + 1):
        qav, _, _ = project(q_cur, h, tau, eps_list)
        spl = CubicSpline(t_grid, qav)
        q_new = integrate_both_ways(lambda t: -Vprime(spl(t)), t0, q0, v0, t_grid, m=m)
        d_new = np.linalg.norm((q_new - q_cur)[core])
        if q_prev is not None:
            d_old = np.linalg.norm((q_cur - q_prev)[core])
            rho = d_new / d_old if d_old > 0 else np.nan
            rhos.append(rho)
            log(f"  iter {it:2d}  |dq| = {d_new:.3e}  rho = {rho:.4f}")
        else:
            log(f"  iter {it:2d}  |dq| = {d_new:.3e}")
        iterates.append(q_new.copy())
        q_prev, q_cur = q_cur, q_new
        if not np.isfinite(d_new) or d_new > 1e6:
            status = "diverged"
            break
        if d_new < tol * max(1.0, np.linalg.norm(q_cur[core])):
            status = "converged"
            break
    return np.array(iterates), np.array(rhos), status


def dominant_frequency(t, q, core):
    """Angular frequency from zero crossings of q - mean on the core window."""
    tt, qq = t[core], q[core] - q[core].mean()
    s = np.sign(qq)
    idx = np.where(s[:-1] * s[1:] < 0)[0]
    if idx.size < 3:
        return np.nan
    # linear interpolation of crossing times
    tc = tt[idx] - qq[idx] * (tt[idx + 1] - tt[idx]) / (qq[idx + 1] - qq[idx])
    half_periods = np.diff(tc)
    return np.pi / half_periods.mean()


# ----------------------------------------------------------------------------------------------------------
# runs
# ----------------------------------------------------------------------------------------------------------
def run_check(log):
    name = "check"
    t_run = time.time()
    log(f"# {name}: kernel vs analytic multiplier, tau->0, eps extrapolation")
    tau = 0.3
    T, h = 400.0, 0.01
    t = (np.arange(int(round(2 * T / h))) - int(round(T / h))) * h
    ok = True
    # (a) multiplier on a cosine
    for nu in [0.5, 1.0, 2.0, 3.0]:
        f = np.cos(nu * t)
        eps_list = [0.16, 0.08, 0.04]
        fav, vals, rich = project(f, h, tau, eps_list)
        core = slice(3 * t.size // 8, 5 * t.size // 8)
        exact = f / (1.0 - tau**2 * nu**2)
        errs = [np.max(np.abs(v[core] - exact[core])) for v in vals]
        err_r = np.max(np.abs(fav[core] - exact[core]))
        log(f"  cos({nu} t): PV multiplier {1/(1-tau**2*nu**2):+.6f}; "
            f"|err| eps={eps_list[0]},{eps_list[1]},{eps_list[2]}: {errs[0]:.2e} {errs[1]:.2e} {errs[2]:.2e}; "
            f"Richardson {err_r:.2e}")
        ok &= err_r < 5e-3
    # (b) tau -> 0 on a Gaussian
    f = np.exp(-t**2 / 2.0)
    for tau_s in [0.3, 0.1, 0.03]:
        fav, _, _ = project(f, h, tau_s, [0.16, 0.08, 0.04])
        core = slice(3 * t.size // 8, 5 * t.size // 8)
        err = np.max(np.abs(fav[core] - f[core]))
        # exact: Gaussian spectrum e^{-nu^2/2} times 1/(1 - tau^2 nu^2): PV integral; the deviation is O(tau^2)
        log(f"  Gaussian, tau = {tau_s}: max|<f> - f| = {err:.3e}  (expect ~ tau^2 <nu^2> = {tau_s**2:.1e})")
        ok &= err < 3.0 * tau_s**2 + 5e-3
    status = "ok" if ok else "FAIL"
    log(f"# END {name} elapsed={time.time() - t_run:.1f}s status={status}")
    return status


def run_picard(log):
    name = "picard"
    t_run = time.time()
    log(f"# {name}: naive fixed-initial-data iteration on V = m w^2 q^2/2 -- recorded because it diverges")
    m, tau = 1.0, 1.0
    out = {}
    T, h = 200.0, 0.02
    t = (np.arange(int(round(2 * T / h))) - int(round(T / h))) * h
    core = slice(3 * t.size // 8, 5 * t.size // 8)
    for wtau in [0.10, 0.30]:
        w = wtau / tau
        disc = 1.0 - 4.0 * tau**2 * w**2
        Omega = np.sqrt((1.0 - np.sqrt(complex(disc))) / (2.0 * tau**2))  # complex if disc < 0
        log(f" w tau = {wtau:.2f}: (4.13) Omega = {Omega:.6f}  (Omega/w = {Omega / w:.6f}), disc = {disc:+.4f}")
        t1 = time.time()
        its, rhos, status = iterate_projection(lambda q: m * w**2 * q, t, 0.0, 1.0, 0.0, tau,
                                               [0.16, 0.08, 0.04], m=m, nmax=8, log=log)
        Om_num = dominant_frequency(t, its[-1], core)
        log(f"  status = {status}, {its.shape[0]-1} iterations, {time.time()-t1:.1f} s; "
            f"measured Omega = {Om_num:.6f}, rel. error vs (4.13) = {abs(Om_num - Omega.real)/abs(Omega):.2e}"
            if status == "converged" else f"  status = {status}, {its.shape[0]-1} iterations, {time.time()-t1:.1f} s")
        if rhos.size:
            log(f"  rho (last three): {rhos[-3:]}")
        out[f"wtau_{wtau:.2f}"] = dict(Omega_exact=Omega, Omega_num=Om_num, rhos=rhos, status=status,
                                       q_final=its[-1], q_0=its[0])
    log("# note: the fixed-initial-data iteration on a window of several periods does not contract (rho is "
        "erratic and ends above 1; secular polynomial terms grow like (w T)^2 per step); see the docstring and "
        "fakeon_series.py for the working method.")
    log(f"# END {name} elapsed={time.time() - t_run:.1f}s status=ok")
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--picard", action="store_true")
    a = ap.parse_args()
    log = lambda s: print(s, flush=True)
    log(f"# fakeon_toy.py  {time.strftime('%Y-%m-%d %H:%M:%S')}  numpy {np.__version__}  python {sys.version.split()[0]}")
    if a.check:
        run_check(log)
    if a.picard:
        run_picard(log)
    if not (a.check or a.picard):
        ap.print_help()


if __name__ == "__main__":
    main()
