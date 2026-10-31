"""bubble_lib.py -- the one-loop bubble of scalar fields under the Feynman, fakeon (average continuation),
Lee-Wick and Feynman-Wheeler prescriptions, in D = 2, 3, 4, computed two independent ways.  Rung 1 of the
perturbative ladder (brief 11 sec. 5b).  Untrusted generator: mpmath, arbitrary precision, no certification here
(that is bubble_cert.py).

Conventions
-----------
Metric (+,-,-,-); propagator 1/(k^2 - a) with a = m^2 complex allowed; p = (p0, 0) in the rest frame; t = p^2.
Units: the Lee-Wick scale M = 1 (so a = +-i for the LW pair of 1801.00915 eq. 5.5) or the physical mass m = 1.

    B_D(t; a1, a2) = int d^D k/(2 pi)^D  1/((k^2 - a1) ((k-p)^2 - a2))

Route (a), Feynman parameters (valid as written in the Euclidean region t < 0, then continued in complex t):
    B_D = i Gamma(2 - D/2)/(4 pi)^{D/2} * I_D(t),   I_D = int_0^1 dx g_D(Delta),  Delta = a1 (1-x) + a2 x - t x (1-x)
    g_2 = 1/Delta,  g_3 = Delta^{-1/2},  g_4 = ln Delta (D = 4 is log divergent; the finite part is -(i/16 pi^2) I_4
    up to an additive constant that cancels in differences I_4(t) - I_4(t_ref)).
The continuation from t < 0 to t > 0 goes around t = 0 along the arc |t| = rho (rho = 1 < 4, inside the main region
of 1703.04584: the (++) and (--) LW terms have branch points at t = +-4 i M^2, and the path must pass between them
and the real axis) in the upper (f_+) or lower (f_-) half plane, then runs along Im t = +-delta to T, stepping
around the LW threshold t = 2 M^2 (a branch point on the second sheet) on a small semicircle.  Along the path the
roots x_1, x_2 of Delta(x) move; the straight-contour integral with the phase of Delta anchored at x = 0
(Delta(0) = a1, principal argument) is single valued off the cuts {t : Delta(x, t) = 0 for some x in [0, 1]} and
jumps across them.  The analytic continuation adds the monodromy of the integrand around the root that crossed the
segment: 2 pi i (1 - x_j) for ln Delta, 2 pi i Res_j for 1/Delta.  The integers are tracked by continuity along the
path (continuation, continuation_sweep; D = 4 and D = 2).  For ordinary (real) masses the only cut is the real
axis above threshold and f_+- are elementary (I_pm_real_masses, all D).  The LW pair in D = 3 (square-root
branch points) is done on the momentum-space route only, against the continued closed forms of 1801.00915 eq. 5.12.
    Average continuation (1801.00915 eq. 4.1/4.4):  f_AV = (f_+ + f_-)/2,   f_d = (f_+ - f_-)/2.

Route (b), momentum space, the nonanalytic Wick rotation of 1703.04584 in the rest frame: the k0 integral along
the imaginary axis by residues (closing to the right, the two poles k0 = w1 and k0 = p0 + w2 with their sum
-(w1 + w2)/(2 w1 w2 [(w1+w2)^2 - p0^2])), with omega_i = sqrt(k^2 + a_i) (Re > 0),
    B_D(p0) = +i int d^{D-1}k/(2 pi)^{D-1}  (w1 + w2) / (2 w1 w2 [(w1 + w2)^2 - p0^2]),
an analytic function of p0 until the two poles pinch the contour at p0 = w1(k) + w2(k) for a real k = k*; there the
LW/fakeon prescription is the principal value in k (1703.04584 eqs. 3.11, 4.8: the average of the two continuations),
the Feynman prescription adds i pi times the residue (p0 -> p0 + i eps; see B_momentum).
The Feynman-Wheeler propagator P 1/(k^2 - m^2) = (1/2)[1/(k^2 - m^2 + i eps) + 1/(k^2 - m^2 - i eps)] gives mixed
terms with one +i eps and one -i eps propagator (2001.01942 eq. 2.3):
    Sigma'(p0) = (i/2) int d^{D-1}k/(2 pi)^{D-1} (w1 - w2) / (w1 w2 [p0^2 - (w1 - w2)^2])   (+ the i eps parts),
whose large-k tail (a1 - a2)/(4 k^3 p0^2) is the nonlocal divergence of 2001.01942 eq. 2.6.

All functions return mpmath numbers; set mp.mp.dps before calling.
"""
import mpmath as mp

PI = mp.pi
I = mp.mpc(0, 1)

# --------------------------------------------------------------------------------------------------------------
# Route (a): Feynman parameters
# --------------------------------------------------------------------------------------------------------------


def delta_roots(t, a1, a2):
    """Roots x_1, x_2 of Delta(x) = a1 (1-x) + a2 x - t x (1-x) = t x^2 + (a2 - a1 - t) x + a1."""
    t, a1, a2 = mp.mpc(t), mp.mpc(a1), mp.mpc(a2)
    A, B, C = t, a2 - a1 - t, a1
    disc = mp.sqrt(B * B - 4 * A * C)
    # numerically stable pair
    if mp.re(B * mp.conj(disc)) >= 0:
        q = -(B + disc) / 2
    else:
        q = -(B - disc) / 2
    x1 = q / A
    x2 = C / q
    return x1, x2


def delta(x, t, a1, a2):
    return a1 * (1 - x) + a2 * x - t * x * (1 - x)


def _anchor_n(t, a1, x1, x2):
    """Integer n with arg(t) + arg(-x1) + arg(-x2) + 2 pi n = Arg(a1) (principal); Delta = t (x-x1)(x-x2)."""
    s = mp.arg(t) + mp.arg(-x1) + mp.arg(-x2)
    return int(mp.nint((mp.arg(a1) - s) / (2 * PI)))


def log_delta_anchored(x, t, a1, a2, roots=None, n=None):
    """ln Delta(x) with the phase continuous along the straight segment [0, 1] and equal to Ln(a1) at x = 0."""
    if roots is None:
        roots = delta_roots(t, a1, a2)
    x1, x2 = roots
    if n is None:
        n = _anchor_n(t, a1, x1, x2)
    d = delta(x, t, a1, a2)
    phase = mp.arg(t) + mp.arg(x - x1) + mp.arg(x - x2) + 2 * PI * n
    return mp.log(abs(d)) + I * phase


def I_straight(D, t, a1, a2, extra_pts=()):
    """int_0^1 dx g_D(Delta) along the straight segment with the anchored phase (route a, one sheet).

    extra_pts: additional subdivision points in (0, 1)."""
    t, a1, a2 = mp.mpc(t), mp.mpc(a1), mp.mpc(a2)
    roots = delta_roots(t, a1, a2)
    n = _anchor_n(t, a1, *roots)
    pts = [mp.mpf(0)]
    for r in list(roots) + list(extra_pts):
        xr = mp.re(r)
        if 0 < xr < 1:
            pts.append(xr)
    pts.append(mp.mpf(1))
    pts = sorted(set(pts))
    if D == 4:
        f = lambda x: log_delta_anchored(x, t, a1, a2, roots, n)
    elif D == 3:
        f = lambda x: mp.exp(-log_delta_anchored(x, t, a1, a2, roots, n) / 2)
    elif D == 2:
        f = lambda x: 1 / delta(x, t, a1, a2)
    else:
        raise ValueError("D must be 2, 3 or 4")
    return mp.quad(f, pts)


def residues_D2(t, a1, a2, roots):
    """Residues of 1/Delta at its roots: 1/Delta'(x_j)."""
    x1, x2 = roots
    dp = lambda x: 2 * t * x + (a2 - a1 - t)
    return 1 / dp(x1), 1 / dp(x2)


def _candidates(D, t, a1, a2, F_s, nmax=2):
    """All F_s + monodromy corrections with integers in [-nmax, nmax]^2, as (value, (n1, n2))."""
    roots = delta_roots(t, a1, a2)
    if D == 4:
        c1, c2 = 2 * PI * I * (1 - roots[0]), 2 * PI * I * (1 - roots[1])
    elif D == 2:
        r1, r2 = residues_D2(t, a1, a2, roots)
        c1, c2 = 2 * PI * I * r1, 2 * PI * I * r2
    else:
        raise ValueError("monodromy tracking implemented for D = 4 and D = 2")
    out = []
    for n1 in range(-nmax, nmax + 1):
        for n2 in range(-nmax, nmax + 1):
            out.append((F_s + n1 * c1 + n2 * c2, (n1, n2)))
    return out


def _lw_thresholds(a1, a2):
    """Real branch points on the second sheet to step around: the threshold t = (sqrt a1 + sqrt a2)^2 when it is
    real (ordinary masses, or a complex-conjugate pair, where it is the LW threshold 2 M^2 for a = +-i M^2)."""
    thr = (mp.sqrt(mp.mpc(a1)) + mp.sqrt(mp.mpc(a2))) ** 2
    if abs(mp.im(thr)) < mp.mpf(10) ** (-(mp.mp.dps - 8)):
        return (mp.re(thr),)
    return ()


def _check_arc(rho, T, a1, a2):
    """The arc |t| = rho around t = 0 must cross the imaginary axis below the nearest branch point of the (a, a)
    type terms, t = 4 a (for a = +-i M^2: +-4 i M^2, the cuts of 1703.04584 that bound the main region)."""
    lim = 4 * min(abs(mp.mpc(a1)), abs(mp.mpc(a2)))
    if not min(mp.mpf(rho), mp.mpf(T)) < lim:
        raise ValueError(f"arc radius {rho} not below the branch point |t| = 4|a| = {lim}: leaves the main region")


def continuation_path(T, upper=True, rho=1, nsteps=160, delta_end=None, dense_at=(2,), targets=()):
    """The t-path used for the continuation: the arc |t| = rho from arg t = pi to arg t = 0+ in the upper
    (upper=True) or lower half plane, then along Im t = +-delta_end from rho to T.  rho < 4 keeps the path in the
    main region of 1703.04584 (below the cut that starts at p^2 = 4 i M^2 for the (++) LW terms); the arc is the
    only way around the branch point t = 0 of the (+-) terms.  dense_at: real points (the LW threshold t = 2 M^2)
    around which the real-axis run is refined.  targets: real values >= rho that must be visited exactly (as
    T + i delta_end).  Returns the ordered list of complex t."""
    T = mp.mpf(T)
    rho = min(mp.mpf(rho), T)
    if delta_end is None:
        delta_end = mp.mpf(10) ** (-(mp.mp.dps - 6))
    sgn = 1 if upper else -1
    th_end = mp.asin(delta_end / rho)
    pts = [mp.mpc(-rho, 0)]
    pts += [rho * mp.expj(sgn * (PI - (PI - th_end) * k / nsteps)) for k in range(1, nsteps)]
    xs = {rho}
    if T > rho:
        xs |= {rho + (T - rho) * k / nsteps for k in range(1, nsteps + 1)}
        xs |= {mp.mpf(x) for x in targets if rho <= mp.mpf(x) <= T}
    detours = []
    for c in dense_at:
        c = mp.mpf(c)
        if rho < c < T:
            # the LW threshold is a branch point on the second sheet (sqrt(t - c), or 1/sqrt in D = 2): step
            # around it on a semicircle of radius r_c instead of along Im t = delta, where the extrapolation
            # across the singularity is ill conditioned
            tg = [mp.mpf(x) for x in targets if abs(mp.mpf(x) - c) > 0]
            r_c = min(mp.mpf(1) / 10, min(abs(x - c) for x in tg) / 2) if tg else mp.mpf(1) / 10
            xs = {x for x in xs if abs(x - c) > r_c}
            # ease in and out of the detour: the function varies as (t - c)^(+-1/2) there
            xs |= {c - r_c * (1 + mp.mpf(k) / 10) for k in range(1, 20)} | {c + r_c * (1 + mp.mpf(k) / 10) for k in range(1, 20)}
            xs = {x for x in xs if rho <= x <= T}
            detours.append((c, r_c))
    real_pts = [mp.mpc(x, sgn * delta_end) for x in sorted(xs)]
    for (c, r_c) in detours:
        arc_c = [c + r_c * mp.expj(sgn * (PI - PI * k / 60)) for k in range(0, 61)]
        # splice: everything with Re < c - r_c, the detour, everything with Re > c + r_c
        before = [p for p in real_pts if mp.re(p) < c - r_c]
        after = [p for p in real_pts if mp.re(p) > c + r_c]
        real_pts = before + arc_c + after
    pts += real_pts
    return pts


def continuation(D, T, a1, a2, upper=True, rho=1, nsteps=160, delta_end=None, log=None, nmax=2, path=None):
    """Analytic continuation of I_D from the Euclidean region to t = T +- i delta_end, T > 0, along
    continuation_path (or the given list of complex t), tracking the monodromy integers by continuity.

    Returns (value, info) with info = dict(ints=(n1, n2), max_ratio=worst best/second-best ambiguity ratio over
    the path, jumps=number of integer changes).  A max_ratio well below 1 means the tracking was unambiguous."""
    T = mp.mpf(T)
    if path is None:
        path = continuation_path(T, upper, rho, nsteps, delta_end, dense_at=_lw_thresholds(a1, a2))
    _check_arc(rho, T, a1, a2)
    prev = None
    prev2 = None
    ints = (0, 0)
    worst = mp.mpf(0)
    jumps = 0
    val = None
    for k, t in enumerate(path):
        F_s = I_straight(D, t, a1, a2)
        if k == 0:
            # Euclidean start: the anchored straight integral is the function itself
            val = F_s
            prev2, prev = None, val
            continue
        cands = _candidates(D, t, a1, a2, F_s, nmax)
        guess = prev if prev2 is None else 2 * prev - prev2
        cands.sort(key=lambda c: abs(c[0] - guess))
        best = cands[0]
        d1 = abs(best[0] - guess)
        # second-best among candidates whose value differs from the best (coincident roots at the LW threshold
        # make two integer vectors give the same value; that is not an ambiguity)
        tol = mp.mpf(10) ** (-3)   # coincident candidates (roots nearly equal near the LW threshold) are not an ambiguity
        d2 = None
        for c in cands[1:]:
            if abs(c[0] - best[0]) > tol:
                d2 = abs(c[0] - guess)
                break
        ratio = d1 / d2 if d2 else mp.mpf(0)
        if ratio > worst:
            worst = ratio
        if best[1] != ints:
            jumps += 1
            ints = best[1]
        val = best[0]
        prev2, prev = prev, val
    if log:
        deg = " [D=2: integer labels degenerate up to (n, n); jumps not meaningful]" if D == 2 else ""
        log(f"    continuation D={D} T={mp.nstr(T, 6)} {'upper' if upper else 'lower'}: ints={ints} "
            f"jumps={jumps}{deg} worst_ratio={mp.nstr(worst, 3)}")
    return val, dict(ints=ints, max_ratio=worst, jumps=jumps)


def average_continuation(D, T, a1, a2, **kw):
    """f_AV, f_d, f_+, f_- at real T > 0 for I_D (route a)."""
    fp, ip = continuation(D, T, a1, a2, upper=True, **kw)
    fm, im = continuation(D, T, a1, a2, upper=False, **kw)
    return (fp + fm) / 2, (fp - fm) / 2, fp, fm, (ip, im)


def continuation_sweep(D, Ts, a1, a2, upper=True, rho=1, nsteps=120, delta_end=None, log=None, nmax=2):
    """continuation() along one path that visits every T in Ts (all > 0): returns {T: (value, info_at_T)}.
    T values below rho get their own small arc (radius T); the rest are read off one pass along the path."""
    Ts = sorted(mp.mpf(T) for T in Ts)
    if delta_end is None:
        delta_end = mp.mpf(10) ** (-(mp.mp.dps - 6))
    rho = mp.mpf(rho)
    out = {}
    small = [T for T in Ts if T < rho]
    big = [T for T in Ts if T >= rho]
    for T in small:
        out[T] = continuation(D, T, a1, a2, upper, rho=T, nsteps=nsteps, delta_end=delta_end, log=log, nmax=nmax)
    if not big:
        return out
    Tmax = big[-1]
    sgn = 1 if upper else -1
    _check_arc(rho, Tmax, a1, a2)
    path = continuation_path(Tmax, upper, rho, nsteps, delta_end, dense_at=_lw_thresholds(a1, a2), targets=big)
    want = {mp.mpc(T, sgn * delta_end): T for T in big}
    prev = prev2 = None
    ints = (0, 0)
    worst = mp.mpf(0)
    worst_t = None
    worst_sep = None
    jumps = 0
    tol = mp.mpf(10) ** (-3)
    val = None
    for k, t in enumerate(path):
        F_s = I_straight(D, t, a1, a2)
        if k == 0:
            val = F_s
            prev2, prev = None, val
        else:
            cands = _candidates(D, t, a1, a2, F_s, nmax)
            guess = prev if prev2 is None else 2 * prev - prev2
            cands.sort(key=lambda c: abs(c[0] - guess))
            best = cands[0]
            d1 = abs(best[0] - guess)
            d2 = None
            for c in cands[1:]:
                if abs(c[0] - best[0]) > tol:
                    d2 = abs(c[0] - guess)
                    break
            ratio = d1 / d2 if d2 else mp.mpf(0)
            if ratio > worst:
                worst, worst_t, worst_sep = ratio, t, (d2 - d1)
            if best[1] != ints:
                jumps += 1
                ints = best[1]
            val = best[0]
            prev2, prev = prev, val
        if t in want:
            out[want[t]] = (val, dict(ints=ints, max_ratio=worst, jumps=jumps, worst_t=worst_t, worst_sep=worst_sep))
    if log:
        deg = " [D=2: Res_1 = -Res_2, so integer labels are degenerate up to (n, n); jumps not meaningful]" if D == 2 else ""
        log(f"    sweep D={D} {'upper' if upper else 'lower'} a=({mp.nstr(a1, 3)},{mp.nstr(a2, 3)}) Tmax={mp.nstr(Tmax, 4)}: "
            f"{len(path)} path points, final ints={ints}{deg} worst_ratio={mp.nstr(worst, 3)} at t={mp.nstr(worst_t, 6) if worst_t is not None else '-'}"
            f" (competing candidates {mp.nstr(abs(worst_sep), 3) if worst_sep is not None else '-'} apart)")
    return out


def I_pm_real_masses(D, T, a1, a2, sign=+1):
    """f_+ (sign=+1) or f_- (sign=-1) at real T for REAL positive a1, a2 (ordinary masses): the exact limit
    t -> T + sign*i0 of the straight-contour integral, which needs no tracking (the only cut is T >= (m1+m2)^2).
    Delta -> Delta_R - sign*i0*x(1-x), so on the segment (x1, x2) where Delta_R < 0:
        ln Delta = ln|Delta_R| - sign*i*pi,   Delta^{-1/2} = |Delta_R|^{-1/2} * exp(+sign*i*pi/2),
        1/Delta = PV 1/Delta_R + sign*i*pi*delta(Delta_R)  ->  + sign*i*pi*sum_j 1/|Delta_R'(x_j)|."""
    T, a1, a2 = mp.mpf(T), mp.mpf(a1), mp.mpf(a2)
    dR = lambda x: a1 * (1 - x) + a2 * x - T * x * (1 - x)
    roots = sorted([mp.re(r) for r in delta_roots(T, a1, a2)]) if T >= (mp.sqrt(a1) + mp.sqrt(a2)) ** 2 else []
    roots = [r for r in roots if 0 < r < 1]
    if len(roots) != 2:
        # below threshold (or the degenerate threshold point): Delta_R > 0 on [0, 1], principal branch
        if D == 4:
            return mp.quad(lambda x: mp.log(dR(x)), [0, mp.mpf(1) / 2, 1])
        if D == 3:
            return mp.quad(lambda x: dR(x) ** mp.mpf(-0.5), [0, mp.mpf(1) / 2, 1])
        return mp.quad(lambda x: 1 / dR(x), [0, mp.mpf(1) / 2, 1])
    x1, x2 = roots
    # above threshold Delta_R = T (x - x1)(x - x2) with 0 < x1 < x2 < 1: everything is elementary
    if D == 4:
        # int_0^1 ln|x - c| dx = (1-c) ln|1-c| + c ln|c| - 1
        Lc = lambda c: (1 - c) * mp.log(1 - c) + c * mp.log(c) - 1
        return mp.log(T) + Lc(x1) + Lc(x2) - sign * I * PI * (x2 - x1)
    if D == 3:
        # int_0^{x1} dx/sqrt(T(x1-x)(x2-x)) = (2/sqrt T) asinh(sqrt(x1/(x2-x1)));  int_{x1}^{x2} = pi/sqrt T (times i);
        # int_{x2}^1 dx/sqrt(T(x-x1)(x-x2)) = (2/sqrt T) asinh(sqrt((1-x2)/(x2-x1)))
        rT = mp.sqrt(T)
        left = (2 / rT) * mp.asinh(mp.sqrt(x1 / (x2 - x1)))
        right = (2 / rT) * mp.asinh(mp.sqrt((1 - x2) / (x2 - x1)))
        return left + right + mp.expj(sign * PI / 2) * PI / rT
    # D = 2: PV int_0^1 dx/(T (x-x1)(x-x2)) = (1/(T(x1-x2))) [ln|(1-x1)/(1-x2)| - ln|x1/x2|],
    # plus sign*i*pi*sum_j 1/|Delta'(x_j)| with Delta'(x1) = T(x1-x2), Delta'(x2) = T(x2-x1)
    pv = (mp.log((1 - x1) / (1 - x2)) - mp.log(x1 / x2)) / (T * (x1 - x2))
    return pv + sign * I * PI * 2 / (T * (x2 - x1))


def _pv_window(sub, lo, c, hi):
    """Integral over the symmetric window [lo, hi] around c of a subtracted integrand that is regular at c but
    cancels catastrophically next to it: Gauss-Legendre on the two halves (its nodes stay O(width/n^2) away from
    the endpoints; tanh-sinh would place nodes within 10^-dps of c)."""
    return (mp.quad(sub, [lo, c], method="gauss-legendre") + mp.quad(sub, [c, hi], method="gauss-legendre"))


# --------------------------------------------------------------------------------------------------------------
# Route (b): momentum space, rest frame
# --------------------------------------------------------------------------------------------------------------


def omega(k, a):
    return mp.sqrt(k * k + a)


def angular(D):
    """d^{D-1}k/(2 pi)^{D-1} = angular(D) * k^{D-2} dk for k = |k| in [0, inf)."""
    if D == 2:
        return 2 / (2 * PI)                     # two half-lines
    if D == 3:
        return 1 / (2 * PI)                     # 2 pi k dk / (2 pi)^2
    if D == 4:
        return 1 / (2 * PI ** 2)                # 4 pi k^2 dk / (2 pi)^3
    raise ValueError


def h_bubble(k, p0, a1, a2):
    """Integrand of B/i after the k0 residues (without the angular factor and k^{D-2})."""
    w1, w2 = omega(k, a1), omega(k, a2)
    s = w1 + w2
    return s / (2 * w1 * w2 * (s * s - p0 * p0))


def pinch_k(p0, a1, a2):
    """Real k* > 0 with w1(k*) + w2(k*) = p0 (real p0 > 0), or None.  For a1, a2 real: closed form; for a
    complex-conjugate pair a2 = conj(a1): w1 + w2 = 2 Re w1 is real for real k and the equation is solved by
    findroot from the real-mass formula seed.  Returns None if no real solution."""
    p0 = mp.mpf(p0)
    if p0 <= 0:
        return None
    # seed from the real-mass formula k*^2 = (p0^4 - 2 p0^2 (a1+a2) + (a1-a2)^2)/(4 p0^2) with Re parts
    s2 = (p0 ** 4 - 2 * p0 ** 2 * mp.re(a1 + a2) + mp.re((a1 - a2) ** 2)) / (4 * p0 ** 2)
    if abs(mp.im(a1)) == 0 and abs(mp.im(a2)) == 0:
        if s2 <= 0 or p0 <= mp.sqrt(a1) + mp.sqrt(a2):
            return None
        return mp.sqrt(s2)
    # a real pinch needs w1 + w2 real for real k: only for a complex-conjugate pair
    if abs(a2 - mp.conj(a1)) > mp.mpf(10) ** (-(mp.mp.dps - 8)):
        return None
    # complex-conjugate pair: is p0 above the pinch threshold (k = 0 value)?
    g = lambda k: mp.re(omega(k, a1) + omega(k, a2)) - p0
    if g(mp.mpf(0)) >= 0:
        return None
    seed = mp.sqrt(s2) if s2 > 0 else p0 / 2
    try:
        ks = mp.findroot(g, seed)
    except Exception:
        ks = mp.findroot(g, (mp.mpf(0), 10 * p0), solver="bisect")
    ks = mp.re(ks)
    if ks <= 0:
        return None
    # the imaginary part of w1 + w2 must vanish there for a genuine pinch on the real axis
    if abs(mp.im(omega(ks, a1) + omega(ks, a2))) > mp.mpf(10) ** (-(mp.mp.dps - 8)):
        raise ValueError("pinch point not on the real k axis: w1 + w2 not real")
    return ks


def B_momentum(D, p0, a1, a2, prescription="PV", p0_ref=None, kmax=None):
    """-i^{-1}... returns B_D(p0; a1, a2) (route b) with the pinch, if any, treated by:
        'PV'       principal value in k (fakeon / Lee-Wick average continuation),
        'feynman'  p0 -> p0 + i eps (adds -i pi * residue),
        'minus'    p0 -> p0 - i eps.
    p0 may be real (> 0) or purely imaginary (Euclidean).  If p0_ref is given, returns the UV-finite difference
    B(p0) - B(p0_ref) (p0_ref Euclidean, no pinch) in one integral; required in D = 4.
    kmax: finite cutoff on |k| (for the Feynman-Wheeler divergence study); default infinity."""
    p0 = mp.mpc(p0)
    ang = angular(D)
    pw = D - 2
    if p0_ref is not None and abs(mp.re(mp.mpc(p0_ref))) != 0:
        raise ValueError("p0_ref must be Euclidean (purely imaginary): the subtraction term has no pinch handling")
    # with no pinch on the real k axis (below threshold, or non-conjugate complex masses) the three prescriptions
    # coincide: the i eps only matters where the contour is pinched
    ks = pinch_k(mp.re(p0), a1, a2) if abs(mp.im(p0)) == 0 else None
    if D == 4 and p0_ref is None and kmax is None:
        raise ValueError("D = 4 needs p0_ref (log divergent) or kmax")
    hi = mp.inf if kmax is None else mp.mpf(kmax)

    def integrand(k):
        v = h_bubble(k, p0, a1, a2)
        if p0_ref is not None:
            v -= h_bubble(k, p0_ref, a1, a2)
        return v * k ** pw

    if ks is None:
        pts = [mp.mpf(0), mp.mpf(1), hi] if hi == mp.inf else [mp.mpf(0), hi]
        # a complex pole of the integrand close to the real k axis (non-conjugate complex masses, e.g. the (++)
        # term of a nearly real LW pair) needs breakpoints: (w1+w2)^2 = p0^2 at k^2 = (p0^4 - 2 p0^2 (a1+a2) + (a1-a2)^2)/(4 p0^2)
        if abs(mp.im(p0)) == 0:
            kt2 = (p0 ** 4 - 2 * p0 ** 2 * (a1 + a2) + (a1 - a2) ** 2) / (4 * p0 ** 2)
            kt = mp.sqrt(kt2)
            if mp.re(kt) > 0 and abs(mp.im(kt)) < mp.mpf(1) / 2:
                r, im = mp.re(kt), abs(mp.im(kt))
                pts += [r - 4 * im, r, r + 4 * im]
                pts = sorted(set(p for p in pts if 0 <= p <= hi))
        val = mp.quad(integrand, pts)
        return I * ang * val
    # pinch: integrand ~ R/(k - k*) with R = N/Dn'(k*)
    w1, w2 = omega(ks, a1), omega(ks, a2)
    s = w1 + w2
    dsdk = ks / w1 + ks / w2
    N = s / (2 * w1 * w2) * ks ** pw
    Dn_prime = 2 * s * dsdk
    R = N / Dn_prime
    wdt = min(ks / 2, mp.mpf(1)) if hi == mp.inf else min(ks / 2, (hi - ks) / 2, mp.mpf(1))

    def sub(k):
        return integrand(k) - R / (k - ks)

    # one interval, Gauss-Legendre: tanh-sinh would place nodes within 1e-dps of k* where the subtraction cancels
    inner = _pv_window(sub, ks - wdt, ks, ks + wdt)
    left = mp.quad(integrand, [mp.mpf(0), ks - wdt])
    if hi == mp.inf:
        right = mp.quad(integrand, [ks + wdt, 2 * ks + 2, mp.inf])
    else:
        right = mp.quad(integrand, [ks + wdt, hi])
    pv = left + inner + right
    # p0 -> p0 + i eps (Feynman) turns Dn = s^2 - p0^2 into Dn - i eps', so 1/Dn -> PV + i pi delta(Dn)
    # = PV + i pi delta(k - k*)/Dn'(k*) with Dn' > 0: the delta term is + i pi R on the integrand of B/i.
    if prescription == "PV":
        val = pv
    elif prescription == "feynman":
        val = pv + I * PI * R
    elif prescription == "minus":
        val = pv - I * PI * R
    else:
        raise ValueError(prescription)
    return I * ang * val


def h_fw_mixed(k, p0, a1, a2):
    """Integrand of the Feynman-Wheeler mixed term Sigma' (one +i eps, one -i eps propagator) after the k0
    residues: Sigma' = (i/2) int d^{D-1}k/(2 pi)^{D-1} (w1 - w2)/(w1 w2 (p0^2 - (w1 - w2)^2))."""
    w1, w2 = omega(k, a1), omega(k, a2)
    d = w1 - w2
    return d / (w1 * w2 * (p0 * p0 - d * d))


def pseudo_pinch_k(p0, a1, a2):
    """Real k* with |w1 - w2| = p0 (pseudothreshold region 0 < p0 < |m1 - m2|), real masses."""
    p0 = mp.mpf(p0)
    m1, m2 = mp.sqrt(a1), mp.sqrt(a2)
    if not (0 < p0 < abs(m1 - m2)):
        return None
    s2 = (p0 ** 4 - 2 * p0 ** 2 * (a1 + a2) + (a1 - a2) ** 2) / (4 * p0 ** 2)
    return mp.sqrt(s2) if s2 > 0 else None


def Sigma_FW_mixed(D, p0, a1, a2, kmax):
    """Sigma'(p0) with a cutoff |k| <= kmax, real masses, real p0 > 0; returns (PV part, delta part), both
    complex, so that Sigma' = PV + delta.  The delta part is the +-i eps contribution at the pseudothreshold pinch
    (2001.01942 eq. 2.3: 1/(X + i eps) = PV - i pi delta(X) on the first term, 1/(Y - i eps) = PV + i pi delta(Y)
    on the second; both give -i pi delta on the bracket)."""
    p0 = mp.mpf(p0)
    ang = angular(D)
    pw = D - 2
    hi = mp.mpf(kmax)
    ks = pseudo_pinch_k(p0, a1, a2)
    f = lambda k: h_fw_mixed(k, p0, a1, a2) * k ** pw
    if ks is None:
        pv = mp.quad(f, [mp.mpf(0), min(hi, mp.mpf(1)), hi] if hi > 1 else [mp.mpf(0), hi])
        return (I / 2) * ang * pv, mp.mpc(0)
    w1, w2 = omega(ks, a1), omega(ks, a2)
    d = w1 - w2
    dddk = ks / w1 - ks / w2
    N = d / (w1 * w2) * ks ** pw
    Dn_prime = -2 * d * dddk            # d/dk [p0^2 - d^2]
    R = N / Dn_prime
    wdt = min(ks / 2, (hi - ks) / 2, mp.mpf(1))
    sub = lambda k: f(k) - R / (k - ks)
    pv = mp.quad(f, [mp.mpf(0), ks - wdt]) + _pv_window(sub, ks - wdt, ks, ks + wdt) + mp.quad(f, [ks + wdt, hi])
    # bracket in eq. 2.3 is (1/4)(1/(w1 w2))[1/(p0 - d + i eps) - 1/(p0 + d - i eps)] * i; the delta parts:
    # -i pi delta(p0 - d) - i pi delta(p0 + d), i.e. for d = +-p0 a term (i/4)(-i pi) delta(.)/(w1 w2) = (pi/4)/(w1 w2)
    # per unit of the delta's argument.  Written through k: delta(p0 - |d|) = delta(k - k*)/|d'(k*)|.
    dpart = (PI / 4) / (w1 * w2) * ks ** pw / abs(dddk) * ang
    return (I / 2) * ang * pv, dpart


# --------------------------------------------------------------------------------------------------------------
# Published closed forms to check against
# --------------------------------------------------------------------------------------------------------------


def lw_closed_D4(t):
    """1801.00915 eqs. 5.6-5.8, M = 1, LW pair a = +-i: (f0, g, f, f_AV) on the real t axis.
    f0 = int ln K, K = (1-2x)^2 + t^2 x^2 (1-x)^2;  g = int ln H, H = 1 + t^2 x^2 (1-x)^2;
    f = f0 + (2 pi/t) theta(t) theta(2-t) (sqrt(4 - t^2) - 2);  f_AV = f - (4 pi/t) theta(t - 2)."""
    t = mp.mpf(t)
    f0 = mp.quad(lambda x: mp.log((1 - 2 * x) ** 2 + t * t * x * x * (1 - x) ** 2), [0, mp.mpf(1) / 2, 1])
    g = mp.quad(lambda x: mp.log(1 + t * t * x * x * (1 - x) ** 2), [0, 1])
    f = f0
    fav = f0
    if 0 < t < 2:
        corr = (2 * PI / t) * (mp.sqrt(4 - t * t) - 2)
        f += corr
        fav += corr
    elif t >= 2:
        fav += -(4 * PI / t)
    return f0, g, f, fav


def lw_closed_D2(t):
    """1801.00915 sec. 5, D = 2, M = 1: f0 = -2t int x(1-x)/K, g = -2t int x(1-x)/H,
    f = f0 + 4 pi/sqrt(4 - t^2) theta(t) theta(2-t), f_AV = f on the whole real axis."""
    t = mp.mpf(t)
    f0 = -2 * t * mp.quad(lambda x: x * (1 - x) / ((1 - 2 * x) ** 2 + t * t * x * x * (1 - x) ** 2),
                          [0, mp.mpf(1) / 2, 1])
    g = -2 * t * mp.quad(lambda x: x * (1 - x) / (1 + t * t * x * x * (1 - x) ** 2), [0, 1])
    f = f0 + (4 * PI / mp.sqrt(4 - t * t) if 0 < t < 2 else 0)
    return f0, g, f, f


def lw_closed_D3_euclid(t):
    """1801.00915 eq. 5.12 region t < 0, M = 1: closed forms of f0 and g (both real for t < 0).
    f0 = (2i/sqrt(-t)) ln((sqrt2 - i sqrt(-t))/(sqrt2 + i sqrt(-t)))
    g  = (i/sqrt(-t)) [ln((2 sqrt(i) + sqrt(-t))/(2 sqrt(i) - sqrt(-t))) + ln((2 sqrt(-i) - sqrt(-t))/(2 sqrt(-i) + sqrt(-t)))]"""
    t = mp.mpf(t)
    assert t < 0
    r = mp.sqrt(-t)
    f0 = (2 * I / r) * mp.log((mp.sqrt(2) - I * r) / (mp.sqrt(2) + I * r))
    si, smi = mp.sqrt(I), mp.sqrt(-I)
    g = (I / r) * (mp.log((2 * si + r) / (2 * si - r)) + mp.log((2 * smi - r) / (2 * smi + r)))
    return f0, g


def fakeon_D4_diff_closed(t, a1, a2):
    """2001.01942 eq. 6.6 (D = 4): Sigma_f - Sigma = sqrt(u+ u-)/(16 pi p^2) theta(-u+),
    u+- = (m1 +- m2)^2 - p^2.  Returned for Sigma = B (no combinatorial 1/2)."""
    t = mp.mpf(t)
    m1, m2 = mp.sqrt(a1), mp.sqrt(a2)
    up, um = (m1 + m2) ** 2 - t, (m1 - m2) ** 2 - t
    if up < 0:
        return mp.sqrt(up * um) / (16 * PI * t)
    return mp.mpf(0)


def feynman_D4_im_closed(t, a1, a2):
    """2001.01942 eq. 4.1: 2 Im(-i Sigma) = sqrt(u+ u-)/(8 pi p^2) theta(-u+)."""
    t = mp.mpf(t)
    m1, m2 = mp.sqrt(a1), mp.sqrt(a2)
    up, um = (m1 + m2) ** 2 - t, (m1 - m2) ** 2 - t
    return mp.sqrt(up * um) / (8 * PI * t) if up < 0 else mp.mpf(0)


def fw_D4_im_closed(t, a1, a2):
    """2001.01942 eq. 4.3: 2 Im(-i Sigma') = -sqrt(u+ u-)/(8 pi p^2) theta(u-)  (for 0 < p^2 < (m1-m2)^2)."""
    t = mp.mpf(t)
    m1, m2 = mp.sqrt(a1), mp.sqrt(a2)
    up, um = (m1 + m2) ** 2 - t, (m1 - m2) ** 2 - t
    return -mp.sqrt(up * um) / (8 * PI * t) if (um > 0 and t > 0) else mp.mpf(0)
