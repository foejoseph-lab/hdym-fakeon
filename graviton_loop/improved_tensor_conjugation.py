#!/usr/bin/env python3
"""
Rung 3, step 3b: 't Hooft-Veltman section 7 (the improved energy-momentum tensor, a = 1/12) -- ROUTE 1, by conjugation.

The claim (ledger row HV-11): with a R phi^2 added, a = 1/12, the one-loop counter-Lagrangian is "of no simple form"; the authors
stop there ("a certain exhaustion").  The author's prediction (3-4 Nov; constant corrected 10 Nov): ON SHELL it is (6.5) rewritten
in the redefined field,
        (sqrt(g)/eps) (203/80) (1/4) K^2 ((d phi)^2)^2,     K = 1/F^2,  F = 1 - phi^2/12,
because the a = 1/12 theory is gravity + a canonical scalar after the Weyl map g_E = F g and phi-tilde = sqrt(12) artanh(phi/sqrt(12))
(weyl_map_check.py: the two Lagrangians differ by a total derivative), and on-shell one-loop divergences do not depend on the
choice of fields.

The method (ledger column "conjugation", ours, 2026; tested against direct derivation on three matter settings in
matter_loop/conjugation_proving.py, the author's run of 10-11 Nov).  THIS run is the second proving step: the same statement for a
field change that also rescales the metric.  It is old mathematics, not a new result: its use here is as a tested tool, run
forward on a case the paper gave up on (author, 10 Nov: the O's are the work).

Set-up, all in exact rationals.
  Einstein frame: exactly step 3b's background (gravscalar_operator_derive.Background): g_E in normal coordinates on a random
  algebraic curvature tensor, phi-tilde(x) a polynomial with rational jets v, W, cubic; only derivatives of phi-tilde enter.
  Original ("Jordan") frame: choose a rational phi_0 with F(phi_0) != 0 and define phi(x) = Phi(phi-tilde(x) - phi-tilde(0)), Phi the
  Taylor series of  dPhi/ds = F(Phi),  Phi(0) = phi_0  (rational coefficients; this is phi-tilde = integral d phi / F without the
  artanh).  The metric is g = g_E / F(phi(x)); the code works with g-hat = (F_0 / F(phi(x))) g_E, which is delta at the origin,
  using  L_J[g-hat / F_0] = (1/F_0) sqrt(g-hat) (F R[g-hat] - 1/2 (d phi)^2_{g-hat})   (a constant rescaling).
  The fluctuations map, to first order, as
        h_E = (F/F_0) h-hat + (F'/F_0) g-hat phi,        phi-tilde-fluct = phi-fluct / F          (F, F' at the background phi(x)),
  an 11 x 11 position-dependent matrix J(x) acting on (h-hat_ab, phi).  The second-order terms of the map multiply the Einstein
  field equations and so drop out of the second variation where those hold.

PREDICTIONS (rule 25; written 11 Nov 07:17 ET, clock read in the same step, before any code below this docstring exists):
  I0 [validation, the Weyl map again but through this script's own Lagrangian evaluators]: at phi_0 = 0 with phi-tilde's jets
      zero (so F = 1 everywhere, J = 1) the Jordan evaluator reproduces the Einstein one: identical ungauged operator.
  I1 [the second proving step; can fail]: ON SHELL at the origin (Ric_E = 1/2 v v, Box phi-tilde = 0), for random phi_0, v, W, cubic,
      the UNGAUGED second variation of the Jordan Lagrangian, derived from that Lagrangian alone by jet probes, equals
            (1/F_0) J^T (ungauged Einstein operator) J
      at the origin on every base field tried (all 11 constant, all 44 linear, and the quadratic ones of the default set [as coded:
      all 110]), all 11 components.  The Einstein side is gravscalar_operator_derive.el_at_origin(..., gauge_sign = 0) on the base field J(x).base.
  I2 [control, must differ]: OFF SHELL (curvature, v, W independent) the two differ on some base field -- the second-order terms
      of the map times the field equations.  Predicted more exactly: they differ already on constant base fields.
  I3 [control]: with the scalar's fluctuation mapped as phi-tilde-fluct = phi-fluct (the 1/F left out; --mutate nof) I1 fails.
  What follows if I1 holds, and by what: taking the gauge condition and the ghosts over from the Einstein frame, the gauge-fixed
  operator of the a = 1/12 theory is the step 3 operator conjugated by J, the ghost operator is unchanged, and by the theorem of
  conjugation_proving.py's header the divergence is that of step 3:  (203/80) sqrt(g_E) R_E^2 on shell
        = (203/320) sqrt(g_E) ((d phi-tilde)^2_{g_E})^2 = (203/320) K^2 sqrt(g) ((d phi)^2_g)^2.
      The last step is an identity (sqrt(g_E) = F^2 sqrt(g), (d phi-tilde)^2_{g_E} = (d phi)^2_g / F^3); the script checks it
      numerically at the origin on each background (I4).
  What this run does NOT do: it does not compute the a = 1/12 divergence in the paper's own gauge and field variables (ROUTE 2,
  the direct computation; that tests gauge independence on shell and is the next script), and it says nothing off shell, where
  the paper's "no simple form" lives.  The local Jacobian of the field change is dropped, as in dimensional regularization --
  the paper's remark on p. 70, ledger row HV-12, still open.
Status PASS if I0, I1 and I4 hold on every on-shell background and I2 shows a difference on every off-shell one.

CORRECTION TO I1, made while writing the code and before any run (11 Nov 15:41 ET, clock read by the script that inserted this line):
  the factor.  The evaluator below is sqrt(g-hat)(F R[g-hat] - 1/2 (d phi)^2_g-hat), WITHOUT the 1/F_0.  Since the Einstein action
  equals (1/F_0) times its integral, its second variation is  F_0 J^T (ungauged Einstein operator) J,  not (1/F_0) as I1 says.
  The prediction as first written stays above (rule 11); the run tests the F_0 form and prints the ratios it finds if they differ.

SMOKE RECORD (rule 11; sandbox 11 Nov, --subset --jobs 2, 1927 s -- over the 20-minute line for smokes; log file time 16:14 ET).
  I0 FAILED as first coded, on exactly one base field: the constant scalar.  The fault was mine, in the test and in the wording of
  I0, not in the evaluators: I had put I0 on a background with random curvature, which with phi-tilde constant is OFF shell (on
  shell means Ricci-flat there), and "F = 1 everywhere, J = 1" overlooked that F(phi_0 + fluctuation) = 1 - fluct^2/12 even at
  phi_0 = 0 -- the second-order term of the map, which gives the scalar the term -(1/12) R fluct^2 and drops out only where the
  Einstein field equations hold.  That is prediction I2 in miniature, and it is the improved tensor itself: (1/12) R phi^2.
  I0 is now run on a Ricci-flat background (the Weyl part random).  The first wording of I0 stays above.
  I1 held on the on-shell background (phi_0 = 1/2, F_0 = 47/48) on all 18 subset base fields, 16 of them with a nonzero result,
  with the factor F_0 of the correction above; I4 matched; I2 as predicted (see the log line for which constants differ).
  Timing: about 97 s per base field on shell (mean, one sandbox core; the Jordan half 54 s), so all 165 base fields are about
  4.5 core-hours; on the author's machine, by the 10 Nov diagnostic and the real runs since, somewhere between 35 minutes and 2 hours.

FRESH REVIEW (11 Nov evening, separate agent; six harness runs, thirteen mutations on scratch copies; repository untouched).
No must-fix: the mathematics and the code are right, and every mutation that can reach the origin is caught.  Folded in or to be read:
  * The on-shell argument holds, and more sharply than I wrote it.  The second-order terms of the map are ALGEBRAIC in the
    fluctuations, so the extra term in the action is (Einstein field equation)(x) times a bilinear with no derivatives, and its
    Euler-Lagrange expression at the origin is EOM_E(0) x B(0)[., base(0)].  No derivative of the field equations can enter: "on
    shell at the origin" is enough for every polynomial base field.  CONSEQUENCE, which I had not seen: for linear and quadratic
    base fields (base(0) = 0) the identity holds OFF shell too (the reviewer checked four of them on the off-shell background).  So
    I2 should read "they differ on CONSTANT base fields only", and the 154 linear and quadratic base fields of I1 test the
    derivative structure of J, not the shell condition.
  * NO DIVERGENCE IS COMPUTED HERE.  The script compares ungauged second variations at the origin and nothing else: no gauge-fixed
    operator of the a = 1/12 theory, no ghost, no E, Omega or trace.  203/80 is step 3's number.  The step from I1 to the
    divergence is an argument, which the reviewer finds correct: the Einstein gauge condition composed with J is a legitimate
    gauge choice; the field-space metric becomes F_0 J^T Pbig J so the operator is J^-1 Delta J; the ghost is unchanged because J
    maps the diffeomorphism of (g-hat, phi) to that of g_E (L_xi((F/F_0) g-hat) = L_xi g_E); the local Jacobian is dropped.  This
    is less than conjugation_proving.py's header promised for "the next run" (the gauge-fixing term carried along in code): the
    gauge term is carried along by argument only.
  * What the other checks carry.  I0 is a regression test of the copied Lagrangian evaluator at F = 1 and of nothing else (it is
    blind to F', the flow, 1/F, J, J^T and the factor F_0; it always uses the 18 subset base fields; only R(0) = 0 is needed for it).
    I2 shows the sensitivity to the shell condition; it does not check that the difference is EOM x B.  I4 is arithmetic at the
    origin (it tests the first flow coefficient and the powers of F_0) and adds almost nothing to I1; its output line is not a result.
  * Which base fields carry what (the mutation table): the constants carry most of the test -- F', the trace term of J^T, g-hat's
    F_0, the second flow coefficient (through the metric constants only), the 1/12; the linear ones carry the derivatives of F, F'
    and 1/F in J; the factor F_0 against 1/F_0 showed on a quadratic one.  Orders x^3, x^4 of phi and g-hat cannot reach the origin,
    so no run tests the third and fourth flow coefficients: they rest on the reviewer's sympy check against sqrt(12) tanh.
  * The square root of det(g-hat + eps h) lacked two terms needed at x^2 (g-hat - delta is O(x) here); they cannot reach the
    origin and the results are identical with them, but they are now included.
  * Timing: my 4.5 core-hours came from the constant-heavy subset; measured per kind (constants ~110 s, linear ~60 s, quadratic
    ~37 s) the full set is about 2.8 core-hours, and scaling from the author's gravscalar_operator_derive_full run suggests about a
    quarter of an hour on his machine.  Forecasts for that machine have been wrong in both directions; the log's per-task line will say.
  * Smoke subset: two of its three quadratic base fields are nearly or wholly vacuous (all-zero results); it concerns smokes only.

Run (the author's machine; every base field):
      python3 improved_tensor_conjugation.py --jobs 16 --log results/improved_tensor_conjugation.log
      control:  python3 improved_tensor_conjugation.py --jobs 16 --subset --mutate nof --offshell 0 --log results/improved_tensor_conjugation_control_nof.log
      (both from graviton_loop; the control needs only the subset: it is caught on constant and linear base fields)
Sandbox smokes use --subset (11 constant, 4 linear, 3 quadratic base fields): they test the code, not the physics.
"""
import sys, os, time, argparse, random, itertools
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gravscalar_operator_derive as D
from gravscalar_operator_derive import idx, XR, EPS, ZERO, ONE, C, T, d, psum, at0_rational, CTX, COMPS, PAIRS

out = None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def Tx(p, xdeg):
    return T(p, xdeg, 0, 0)


def flow_coefficients(phi0, nmax=4):
    """Taylor coefficients c_0..c_nmax of Phi(s) with dPhi/ds = F(Phi) = 1 - Phi^2/12, Phi(0) = phi0 (exact rationals)"""
    c = [sp.Rational(phi0)]
    for n in range(nmax):
        sq = sum(c[i] * c[n - i] for i in range(n + 1))                    # [s^n] Phi^2
        c.append(((1 if n == 0 else 0) - sq / 12) / (n + 1))
    return c


class Jordan:
    """the a = 1/12 theory on the background that the Weyl map sends to the Einstein-frame background bgE"""
    def __init__(self, bgE, phi0, mutate='none'):
        self.mutate = mutate
        s = bgE.phit                                                       # phi-tilde(x) - phi-tilde(0): a polynomial without constant term
        cs = flow_coefficients(phi0)
        phi = C(cs[0]); pw = ONE
        for n in range(1, len(cs)):
            pw = Tx(pw * s, 4); phi = phi + C(cs[n]) * pw
        self.phi = Tx(phi, 4)
        self.F = Tx(ONE - self.phi * self.phi * C(sp.Rational(1, 12)), 4)
        self.Fp = self.phi * C(-sp.Rational(1, 6))                         # F'(phi)
        self.F0 = at0_rational(Tx(self.F, 0)); self.Fp0 = at0_rational(Tx(self.Fp, 0)); self.phi0 = sp.Rational(phi0)
        assert self.F0 != 0
        u = self.F * C(1 / self.F0) - ONE
        acc = ONE; term = ONE
        for _ in range(4):
            term = Tx(-term * u, 4); acc = acc + term
        self.Finv = acc * C(1 / self.F0)                                   # 1/F(phi(x)) to x^4
        self.ghat = [[Tx(self.Finv * bgE.gb[a][b] * C(self.F0), 4) for b in idx] for a in idx]      # (F_0/F) g_E: delta at the origin
        self.dphi = [d(self.phi, a) for a in idx]

    def L2(self, h, ph, c_gauge=None, gauge_sign=None):
        """the eps^2 coefficient of sqrt(G)(F(Phi) R[G] - 1/2 G^{mn} d Phi d Phi), G = g-hat + eps h, Phi = phi + eps ph, to x^2.
        No gauge term.  (Signature as Background.L2 so that gravscalar_operator_derive.el_at_origin can be used unchanged.)"""
        G = [[self.ghat[a][b] + EPS * h[a][b] for b in idx] for a in idx]
        Ginv = D.series_inverse(G, 3)
        A = [[G[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]; A2 = D.mat_mul(A, A, 2); A3 = D.mat_mul(A2, A, 2)
        tr = lambda M: psum(M[i][i] for i in idx)
        A4 = D.mat_mul(A3, A, 2)             # g-hat - delta is O(x) here (not O(x^2) as in the Einstein frame): one more order (review, 11 Nov)
        logdet = T(tr(A) - tr(A2) * C(sp.Rational(1, 2)) + tr(A3) * C(sp.Rational(1, 3)) - tr(A4) * C(sp.Rational(1, 4)), 2)
        l2 = T(logdet * logdet, 2)
        sqrtg = T(ONE + logdet * C(sp.Rational(1, 2)) + l2 * C(sp.Rational(1, 8)) + l2 * logdet * C(sp.Rational(1, 48)) + l2 * l2 * C(sp.Rational(1, 384)), 2)
        dG = [[[T(d(G[a][b], c), 3) for c in idx] for b in idx] for a in idx]
        Gam = [[[T(psum(Ginv[m][l] * (dG[l][b][a] + dG[l][a][b] - dG[a][b][l]) for l in idx) * C(sp.Rational(1, 2)), 3)
                 for b in idx] for a in idx] for m in idx]
        Ric = [[T(psum(d(Gam[a][b][c], a) - d(Gam[a][b][a], c) for a in idx)
                  + psum(Gam[a][a][e] * Gam[e][b][c] - Gam[a][c][e] * Gam[e][b][a] for a in idx for e in idx), 2) for c in idx] for b in idx]
        Rsc = T(psum(Ginv[b][c] * Ric[b][c] for b in idx for c in idx), 2)
        Phi = T(self.phi + EPS * ph, 2)
        FPhi = T(ONE - Phi * Phi * C(sp.Rational(1, 12)), 2)
        dPhi = [T(self.dphi[a] + EPS * d(ph, a), 2) for a in idx]
        kin = T(psum(Ginv[m][n] * dPhi[m] * dPhi[n] for m in idx for n in idx), 2)
        L = T(sqrtg * (T(FPhi * Rsc, 2) - kin * C(sp.Rational(1, 2))), 2)
        return CTX.from_dict({m[:4] + (0,) + m[5:]: c for m, c in L.to_dict().items() if m[4] == 2})

    def J(self, base):
        """(h-hat, phi) -> (h_E, phi-tilde fluctuation) to first order in the fluctuations, as fields to x^2"""
        h, ph = base
        hE = [[Tx((self.F * h[a][b] + self.Fp * self.ghat[a][b] * ph) * C(1 / self.F0), 2) for b in idx] for a in idx]
        phE = ph if self.mutate == 'nof' else Tx(self.Finv * ph, 2)
        return hE, phE

    def JT(self, r):
        """the transpose at the origin, on the 11 Euler-Lagrange components (order COMPS; pairing sum_ab h_ab r_ab + phi r_s)"""
        r = list(r)
        tr = sum(r[COMPS.index((a, a))] for a in idx)
        rs = r[10] if self.mutate == 'nof' else r[10] / self.F0
        return r[:10] + [self.Fp0 / self.F0 * tr + rs]


_W = {}


def _init(Rd, v, W, cub, phi0, mutate):
    _W['bgE'] = D.Background(Rd, v, W, cub); _W['J'] = Jordan(_W['bgE'], phi0, mutate); _W['probes'] = D.probe_fields()


def _task(spec):
    t = time.time()
    base = D.base_from_spec(spec)
    jo = _W['J']
    lhs = D.el_at_origin(jo, base, 0, _W['probes'], gauge_sign=0)                       # the Jordan Lagrangian, from itself
    tj = time.time() - t
    rE = D.el_at_origin(_W['bgE'], jo.J(base), 0, _W['probes'], gauge_sign=0)           # ungauged Einstein operator on J.base
    rhs = [jo.F0 * x for x in jo.JT(rE)]
    return lhs, rhs, tj, time.time() - t


def base_specs(subset):
    mono = lambda *ms: tuple(sum(1 for m in ms if m == i) for i in idx)
    consts = [('const', [(k, (0, 0, 0, 0), 1)]) for k in range(11)]
    lin = [('linear', [(k, mono(m), 1)]) for k in range(11) for m in idx]
    quad = [('quadratic', [(k, mono(m, n), sp.Rational(1, 2) if m == n else 1)]) for k in range(11) for m in idx for n in idx if m <= n]
    if subset:
        lin = [lin[1], lin[14], lin[27], lin[42]]            # three metric components and the scalar, different directions
        quad = [quad[3], quad[57], quad[105]]
    return consts, lin, quad


def compare(Rd, v, W, cub, phi0, mutate, specs, jobs):
    t0 = time.time()
    if jobs > 1:
        import multiprocessing as mp
        with mp.Pool(min(jobs, len(specs)), initializer=_init, initargs=(Rd, v, W, cub, phi0, mutate)) as pool:
            res = []
            for k_, r_ in enumerate(pool.imap(_task, [s_[1] for s_ in specs], chunksize=1)):
                res.append(r_)
                if (k_ + 1) % 20 == 0:
                    say(f"      ... {k_ + 1} of {len(specs)} base fields done, {sum(1 for (sp_, rr_) in zip(specs, res) if list(rr_[0]) != list(rr_[1]))} mismatches so far, {time.time() - t0:.0f} s")
    else:
        _init(Rd, v, W, cub, phi0, mutate)
        res = [_task(s_[1]) for s_ in specs]
    wall = time.time() - t0
    bad = {}; tried = {}
    ratios = set()
    for (kind, _), (lhs, rhs, tj, tt) in zip(specs, res):
        tried[kind] = tried.get(kind, 0) + 1
        if list(lhs) != list(rhs):
            bad[kind] = bad.get(kind, 0) + 1
        for a_, b_ in zip(lhs, rhs):
            if b_ != 0:
                ratios.add(sp.Rational(a_) / b_)
    nonzero = sum(1 for (lhs, rhs, tj, tt) in res if any(x != 0 for x in lhs))
    tts = [r[3] for r in res]; tjs = [r[2] for r in res]
    timing = f"tasks {len(res)}: per task min/mean/max {min(tts):.0f}/{sum(tts) / len(tts):.0f}/{max(tts):.0f} s (Jordan part mean {sum(tjs) / len(tjs):.0f} s); sum {sum(tts):.0f} s, wall {wall:.0f} s, ratio {sum(tts) / wall:.1f}"
    return bad, tried, nonzero, ratios, timing


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261111)
    ap.add_argument('--onshell', type=int, default=1)
    ap.add_argument('--offshell', type=int, default=1)
    ap.add_argument('--subset', action='store_true', help='sandbox smokes only: 11 constant, 4 linear, 3 quadratic base fields instead of all 165 (never status=PASS)')
    ap.add_argument('--mutate', default='none', choices=['none', 'nof'], help="control: 'nof' leaves the 1/F out of the scalar fluctuation's map (I3; must FAIL)")
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; to start over:  del {args.log.replace('/', chr(92))}"); sys.exit(1)
        out = open(args.log, 'w')
    say(f"# improved_tensor_conjugation: seed {args.seed}, {args.onshell} on-shell + {args.offshell} off-shell backgrounds, base fields {'SUBSET (smoke)' if args.subset else 'all 165'}, mutate {args.mutate}, jobs {args.jobs}")
    say(D.version_line())
    rng = random.Random(args.seed)
    consts, lin, quad = base_specs(args.subset)
    sub_c, sub_l, sub_q = base_specs(True)
    ok = True
    # ---- I0: F = 1 (phi_0 = 0, phi-tilde constant): the Jordan evaluator must reproduce the Einstein one -- on a Ricci-flat background
    R = D.with_ricci(D.riemann_sample(rng), sp.zeros(4, 4))
    bad, tried, nz, ratios, timing = compare(R, sp.zeros(4, 1), sp.zeros(4, 4), {}, 0, 'none', sub_c + sub_l + sub_q, args.jobs)
    i0 = not bad and nz > 0
    say(f"I0 F = 1, Ricci-flat background: Jordan evaluator == Einstein evaluator on {sum(tried.values())} base fields ({nz} with a nonzero result): {'PASS' if i0 else 'FAIL ' + str(bad)}   [{timing}]")
    ok &= i0
    i1_all = True; i2_all = True; i4_all = True; n_on = n_off = 0
    for kind, n in (('on shell', args.onshell), ('off shell', args.offshell)):
        for i in range(n):
            R = D.riemann_sample(rng); v, W, cub = D.scalar_sample(rng)
            phi0 = sp.Rational(rng.choice((-3, -2, -1, 1, 2, 3)), rng.choice((2, 3)))             # |phi_0| <= 3/2 < sqrt(12)
            if kind == 'on shell':
                W = W - W.trace() * sp.eye(4) / 4
                R = D.with_ricci(R, v * v.T / 2)
            Ric, Rs = D.ricci_of(R)
            onsh = (Ric - v * v.T / 2) == sp.zeros(4, 4) and W.trace() == 0
            specs = (consts + lin + quad) if kind == 'on shell' else consts                    # off shell: the constants suffice to show the difference
            bad, tried, nz, ratios, timing = compare(R, v, W, cub, phi0, args.mutate, specs, args.jobs)
            jo = Jordan(D.Background(R, v, W, cub), phi0)
            F0 = jo.F0
            say(f"{kind} {i}: phi_0 = {phi0}, F_0 = {F0}, F'_0 = {jo.Fp0}; R_E = {Rs}, (d phi-tilde)^2 = {(v.T * v)[0]}, Einstein field equations hold at 0: {onsh}")
            say(f"   base fields tried {tried}; with a nonzero result {nz}; Jordan != F_0 J^T (Einstein) J on: {bad if bad else 'none'}"
                + (f";  ratios lhs/rhs seen: {sorted(ratios, key=str)[:6]}" if bad else "") + f"   [{timing}]")
            if kind == 'on shell':
                n_on += 1; i1_all &= (not bad) and onsh and nz > 0
                # I4: (203/80) sqrt(gE) R_E^2 = (203/320) K^2 sqrt(g) ((d phi)^2_g)^2 at the origin, from the Jordan background's own fields
                dphi0 = [at0_rational(Tx(jo.dphi[a], 0)) for a in idx]
                sqrtg = 1 / F0**2                                  # g = g-hat / F_0, g-hat(0) = delta
                dphi2_g = F0 * sum(x * x for x in dphi0)           # g^{-1} = F_0 delta
                lhs4 = sp.Rational(203, 80) * Rs**2                # sqrt(gE)(0) = 1
                rhs4 = sp.Rational(203, 320) * (1 / F0**2)**2 * sqrtg * dphi2_g**2
                i4 = sp.expand(lhs4 - rhs4) == 0
                say(f"   I4: (203/80) R_E^2 = {lhs4};  (203/320) K^2 sqrt(g) ((d phi)^2_g)^2 = {rhs4}: {'MATCH' if i4 else 'MISMATCH'}")
                i4_all &= i4
            else:
                n_off += 1; i2_all &= bool(bad)
    say("")
    if n_on:
        say(f"I1 on shell, Jordan second variation == F_0 J^T (ungauged Einstein operator) J, every base field tried, all 11 components{' (mutated: must FAIL)' if args.mutate != 'none' else ''}: {'PASS' if i1_all else 'FAIL'}")
        say(f"I4 (arithmetic only, not a result: step 3's (6.5) rewritten in the original variables at the origin): {'PASS' if i4_all else 'FAIL'}")
        ok &= i1_all and i4_all
    if n_off:
        say(f"I2 off shell they differ on the constant base fields (the only ones that can: linear and quadratic base fields vanish at the origin): {'as predicted' if i2_all else 'NOT as predicted'}")
        ok &= i2_all
    complete = n_on > 0 and n_off > 0 and not args.subset and args.mutate == 'none'
    if not ok:
        say("\nstatus=FAIL")
    elif not complete:
        say("\nstatus=INCOMPLETE (PASS needs every base field, no mutation, and at least one on-shell and one off-shell background)")
    else:
        say("\nstatus=PASS")
    say("No divergence is computed by this script; what I1 licenses, and by which argument, is in the header.")
    say(f"END improved_tensor_conjugation {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
