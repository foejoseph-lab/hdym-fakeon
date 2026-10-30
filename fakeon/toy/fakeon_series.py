"""fakeon_series.py -- the perturbative fakeon projection order by order (Anselmi 1809.05037 sec. 4,
1901.09273 sec. 5 for the FLRW analogue), path B step 1.

Projected equation:  m q'' = -V'(<q>),   <q> = 1/(1 + tau^2 d^2/dt^2) q = sum_k (-tau^2)^k q^(2k)   (formal).
Order-by-order reduction (what is done on a gravitational background):  with q'' = f_{K-1}(q, q') known to
O(tau^{2K-2}), compute q^(2k) by differentiating along the flow, substitute, truncate at tau^{2K}:
    f_K(q, q') = -(1/m) V'( q + sum_{k=1}^{K} (-tau^2)^k D^{2k-2} f_{K-1} )  mod tau^{2K+2}.
Each f_K is an ordinary second-order ODE with two initial data: the projection by construction.  The
convergence of the projection is the convergence of the trajectories q_K(t) in K.

Harmonic check:  V = m w^2 q^2/2 gives f_K = -Omega_K^2 q with Omega_K^2 = sum_{k<=K} C_k tau^{2k} w^{2k+2}
(C_k the Catalan numbers), the Taylor series of eq. (4.13), Omega^2 = (1 - sqrt(1 - 4 tau^2 w^2))/(2 tau^2),
radius of convergence tau w = 1/2 -- the same threshold at which the projected solution space becomes empty.

Wall:  V = V0 exp((q - qw)/w), nothing else.  The bounce has duration ~ w/v; the observables are the outgoing
velocity v_out,K (conserved after the bounce) and the position at t_end; the convergence of the projection is
read from dq_K = |q_K - q_{K-1}| at t_end (the bounce time delay) versus K at fixed tau v/w.

Runs (Spyder):
    import sys
    %cd C:\\GitHub\\hdym\\fakeon\\toy
    !"{sys.executable}" -B -u fakeon_series.py --harmonic --kmax 6    > results\\series_harmonic.log
    !"{sys.executable}" -B -u fakeon_series.py --wall --kmax 7        > results\\series_wall.log
Long runs (K >= 8; the reduction grows x8 per order, ~2 min per width at K = 8, ~16 min at K = 9) are resumable:
    !"{sys.executable}" -B -u fakeon_series.py --wall --kmax 9 --out results\\series_wall_k9.npz --resume >> results\\series_wall_k9.log
The npz is checkpointed after every wall width; --resume skips the widths already in --out (refused if that file
was written with a different kmax).  Use >> for the log of a resumable run, > otherwise.
Outputs: results/series_<name>.npz (signed), the log; every log ends '# END <name> status=...'.
"""
import argparse
import os
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

q, p, s, E = sp.symbols("q p s E", real=True)  # s = tau^2; E stands for the wall exponential exp((q - qw)/w)
WALL_W = sp.Symbol("w", positive=True)             # d/dq E = E / w


def truncate(expr, K):
    """Drop powers of s above s^K (expr polynomial in q, p, E, s)."""
    P = sp.Poly(sp.expand(expr), s)
    return sum(c * s**k for (k,), c in P.terms() if k <= K)


def d_q(expr, w=WALL_W):
    """Total q-derivative, with the generator rule d/dq E = E/w."""
    return sp.diff(expr, q) + sp.diff(expr, E) * E / w


def flow_derivative(expr, f, w=WALL_W):
    """d/dt of expr(q, p, E) along q' = p, p' = f(q, p, E)."""
    return d_q(expr, w) * p + sp.diff(expr, p) * f


def reduce_order(Vprime, m, K, log=print, w=WALL_W):
    """Return [f_0, ..., f_K], the reduced right-hand sides q'' = f_K(q, p, E; s), each truncated at s^K.

    Vprime(q, E) is a polynomial in q and E (E = exp((q - qw)/w), or unused); w may be a symbol or a number
    (a number keeps the coefficients numeric and the reduction several times faster)."""
    fs = [sp.expand(-Vprime(q, E) / m)]
    for Kcur in range(1, K + 1):
        f_prev = fs[-1]
        t1 = time.time()
        derivs = [f_prev]
        for _ in range(1, Kcur):
            d = flow_derivative(flow_derivative(derivs[-1], f_prev, w), f_prev, w)
            derivs.append(truncate(d, Kcur))
        delta = truncate(sum((-s) ** k * derivs[k - 1] for k in range(1, Kcur + 1)), Kcur)
        # V'(q + delta) = sum_j (d_q^j V')(q) delta^j / j!, truncated at s^Kcur (delta = O(s))
        Vp = sp.expand(Vprime(q, E))
        ser = 0
        dj = Vp
        dpow = 1
        for j in range(Kcur + 1):
            ser += dj * dpow / sp.factorial(j)
            dj = sp.expand(d_q(dj, w))
            dpow = truncate(dpow * delta, Kcur)
        fK = truncate(-ser / m, Kcur)
        fs.append(fK)
        log(f"  reduced order K = {Kcur}: {sp.count_ops(fK)} ops, {time.time() - t1:.1f} s")
    return fs


def integrate(fK_num, tau, q0, v0, t_end, t_eval, rtol=1e-11, atol=1e-13, qw=None, wv=None):
    s_val = tau**2

    budget = [int(2e6)]  # rhs evaluations; a stiff or runaway iterate is a failure, not a long integration

    def rhs(t, y):
        budget[0] -= 1
        if budget[0] < 0:
            raise RuntimeError("evaluation budget exhausted")
        Ev = np.exp((y[0] - qw) / wv) if qw is not None else 0.0
        return [y[1], fK_num(y[0], y[1], Ev, s_val)]

    def blowup(t, y):  # terminal: a runaway iterate is a failure, not a long integration
        return 1e3 - max(abs(y[0]), abs(y[1]))

    blowup.terminal = True
    sol = solve_ivp(rhs, (0.0, t_end), [q0, v0], t_eval=t_eval, method="DOP853", rtol=rtol, atol=atol,
                    events=blowup)
    return sol


def run_harmonic(kmax, log):
    name = "series_harmonic"
    t_run = time.time()
    log(f"# {name}: Omega_K^2 against the Catalan series of (4.13); trajectories q_K vs exact cos(Omega t)")
    m = 1.0
    w = sp.Symbol("omega", positive=True)
    fs = reduce_order(lambda x, Ex: m * w**2 * x, m, kmax, log)
    ok = True
    catalan = [sp.catalan(k) for k in range(kmax + 1)]
    for K, fK in enumerate(fs):
        coeffs = [sp.simplify(-sp.expand(fK).coeff(s, k) / q) for k in range(K + 1)]
        expect = [catalan[k] * w ** (2 * k + 2) for k in range(K + 1)]
        match = all(sp.simplify(c - e) == 0 for c, e in zip(coeffs, expect))
        ok &= match
        log(f"  K = {K}: Omega_K^2 = {' + '.join(str(c) for c in coeffs)}   Catalan match: {match}")
    tau = 1.0
    out = {}
    for wtau in [0.1, 0.3, 0.45, 0.49, 0.55]:
        wv = wtau / tau
        disc = 1 - 4 * tau**2 * wv**2
        Om = np.sqrt((1 - np.sqrt(complex(disc))) / (2 * tau**2))
        t_end = 20.0 * 2 * np.pi / wv  # twenty base periods
        t_eval = np.linspace(0, t_end, 4001)
        errs = []
        for K, fK in enumerate(fs):
            fnum = sp.lambdify((q, p, E, s), fK.subs(w, wv), "numpy")
            sol = integrate(fnum, tau, 1.0, 0.0, t_end, t_eval)
            if disc > 0:
                errs.append(np.max(np.abs(sol.y[0] - np.cos(Om.real * t_eval))))
            else:
                errs.append(np.nan)
        OmK = [float(sum(float(catalan[k]) * tau ** (2 * k) * wv ** (2 * k + 2) for k in range(K + 1))) ** 0.5
               for K in range(kmax + 1)]
        log(f" w tau = {wtau}: Omega (4.13) = {Om:.6f}; Omega_K = {np.array2string(np.array(OmK), precision=6)}")
        if disc > 0:
            log(f"   max|q_K - cos(Omega t)| over 20 periods, K = 0..{kmax}: {np.array2string(np.array(errs), precision=2)}")
        else:
            log("   disc < 0: no real projected solution; the Omega_K series still has finite terms "
                "(partial sums of a divergent series)")
        out[f"wtau_{wtau}"] = dict(Omega=Om, OmegaK=np.array(OmK), errs=np.array(errs))
    np.savez(os.path.join(RES, f"{name}.npz"), **{k: np.array(v, dtype=object) for k, v in out.items()})
    status = "ok" if ok else "FAIL"
    log(f"# END {name} elapsed={time.time() - t_run:.1f}s status={status}")


def run_wall(kmax, log, out=None, resume=False):
    name = "series_wall"
    t_run = time.time()
    out = out or os.path.join(RES, f"{name}.npz")
    log(f"# {name}: one bounce off V0 exp((q - qw)/w), no other force; convergence in K vs tau v/w  "
        f"kmax {kmax}  out {out}  resume {resume}")
    # V0 = 1/2 puts the turning point of a unit-speed particle at qw; launched from eight widths out
    # (wall force e^{-8} of its peak), run to eight widths out again; nothing else acts, so the outgoing
    # velocity v_out and the position at t_end are the observables (v_out is conserved after the bounce); the
    # projection mainly delays the bounce, so |dq| at t_end carries the signal and |dv_out| sits at machine
    # precision from K ~ 5; the two agree where both are resolved
    m, tau, qw, V0, v0 = 1.0, 0.2, 1.0, 0.5, 1.0
    done = {}
    if resume and os.path.exists(out):
        with np.load(out, allow_pickle=True) as z:
            if int(z["kmax"]) != kmax:
                log(f"# resume refused: {out} was written with kmax {int(z['kmax'])}, this run is kmax {kmax}")
                log(f"# END {name} elapsed={time.time() - t_run:.1f}s status=FAIL")
                return
            done = {k: z[k].item() for k in z.files if k.startswith("ratio_")}
        log(f"# resumed: {len(done)} wall width(s) already in {out}: {sorted(done)}")
    fmt = {"float_kind": lambda x: f"{x:.3e}"}

    def checkpoint(results):
        np.savez(out, tau=tau, qw=qw, V0=V0, v0=v0, kmax=kmax,
                 **{k: np.array(v, dtype=object) for k, v in results.items()})

    for ratio in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]:
        key = f"ratio_{ratio}"
        if key in done:
            log(f" tau v/w = {ratio:.2f}: in checkpoint, skipped")
            continue
        wv = tau * v0 / ratio
        q0 = qw - 8.0 * wv
        t_end = 20.0 * wv / v0
        t_eval = np.linspace(0, t_end, 2001)
        t1 = time.time()
        fs = reduce_order(lambda x, Ex: (V0 / wv) * Ex, m, kmax, log=lambda _: None, w=sp.Float(wv))
        tred = time.time() - t1
        trajs, vels, status = [], [], []
        for K, fK in enumerate(fs):
            fnum = sp.lambdify((q, p, E, s), fK, "numpy")
            try:
                sol = integrate(fnum, tau, q0, v0, t_end, t_eval, rtol=1e-13, atol=1e-15, qw=qw, wv=wv)
                ok = sol.success and sol.y.shape[1] == t_eval.size and np.all(np.isfinite(sol.y))
            except Exception as exc:
                ok = False
                log(f"   K = {K}: integration failed ({str(exc)[:40]})")
            trajs.append(sol.y[0].copy() if ok else np.full(t_eval.size, np.nan))
            vels.append(sol.y[1][-1] if ok else np.nan)
            status.append("ok" if ok else "fail")
        trajs, vels = np.array(trajs), np.array(vels)
        d_v = np.abs(np.diff(vels))                      # |v_out,K - v_out,K-1|: tiny, at machine precision by K ~ 5
        d_q_end = np.abs(np.diff(trajs[:, -1]))          # |q_K - q_K-1| at t_end = v_out x (bounce time delay): primary
        r = d_q_end[1:] / d_q_end[:-1]
        c = r / ratio**2
        turn = np.nanmax(trajs[0])
        log(f" w = {wv:.3f}  tau v/w = {ratio:.2f}  reduction {tred:.0f} s  turning point q = {turn:.4f} (wall at {qw})"
            f"  status {status}")
        log(f"   v_out,K, K = 0..{kmax}: {np.array2string(vels, formatter={'float_kind': lambda x: f'{x:+.10f}'})}")
        log(f"   d_K = |v_out,K - v_out,K-1|, K = 1..{kmax}: {np.array2string(d_v, formatter=fmt)}")
        log(f"   |dq| at t_end, K = 1..{kmax}:           {np.array2string(d_q_end, formatter=fmt)}")
        log(f"   r_K = dq_K/dq_(K-1), K = 2..{kmax}: {np.array2string(r, formatter=fmt)}"
            f"   c_K = r_K/(tau v/w)^2: {np.array2string(c, formatter=fmt)}"
            f"   increments of c_K: {np.array2string(np.diff(c), formatter=fmt)}")
        done[key] = dict(w=wv, tau_rate=ratio, t=t_eval, q0=q0, trajs=trajs, v_out=vels, d_v=d_v,
                         d_q_end=d_q_end, r=r, c=c, status=status)
        checkpoint(done)  # one width per checkpoint: a killed run loses at most the width in progress
    log("# reading: a convergent series has c_K -> c_inf < inf (Catalan: 1, 2, 2.5, 2.8, ... -> 4, increments "
        "decreasing); the local-harmonic radius at the turning point, V'' = V0 E/w^2 = 1/(2 w^2), would give "
        "c_inf = 2.  Increasing increments of c_K, and c_K beyond 2, point to an asymptotic series; a Domb-Sykes "
        "(r_K vs 1/K) extrapolation needs K >= 7 to be decisive.")
    checkpoint(done)
    log(f"# END {name} elapsed={time.time() - t_run:.1f}s status=ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--harmonic", action="store_true")
    ap.add_argument("--wall", action="store_true")
    ap.add_argument("--kmax", type=int, default=5)
    ap.add_argument("--out", default=None, help="npz name for --wall (default results/series_wall.npz)")
    ap.add_argument("--resume", action="store_true",
                    help="--wall: load --out if it exists and skip the wall widths already in it")
    a = ap.parse_args()
    log = lambda m: print(m, flush=True)
    log(f"# fakeon_series.py  {time.strftime('%Y-%m-%d %H:%M:%S')}  numpy {np.__version__} sympy {sp.__version__} "
        f"python {sys.version.split()[0]}  kmax {a.kmax}")
    if a.harmonic:
        run_harmonic(a.kmax, log)
    if a.wall:
        run_wall(a.kmax, log, out=a.out, resume=a.resume)
    if not (a.harmonic or a.wall):
        ap.print_help()


if __name__ == "__main__":
    main()
