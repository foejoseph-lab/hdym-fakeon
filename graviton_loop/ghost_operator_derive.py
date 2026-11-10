#!/usr/bin/env python3
"""
Rung 3, step 3b (part two): the Faddeev-Popov ghost operator of gravity + a massless scalar DERIVED, not transcribed -- the piece
the fresh review of 9 Nov night left open (gravscalar_operator_derive.py derives the (h, phi) operator only; the ghost's
E = Ric - d phi-tilde d phi-tilde (5.19) was still from the print in hv1974_gravscalar_M.py and rung3_dewitt_bundle.py).

Two things are computed, with the machinery of gravscalar_operator_derive.py (imported: same background, same Lagrangian, same
Euler-Lagrange jet formula; python-flint, exact rationals):

 (a) THE SYMMETRY, from the action.  The ghost operator is the variation of the gauge condition under the gauge symmetry, so the
     symmetry's form on the fluctuations -- and in particular the RELATIVE normalisation of the metric's and the scalar's
     transformation, which is what puts the v v term into the ghost -- must come from the action and not from memory:
         delta h_ab = D_a eta_b + D_b eta_a,     delta phi = k eta^a d_a phi-tilde        (k = 1 claimed; background derivatives)
     is a symmetry only if the UNGAUGED second-variation operator (gauge term switched off) annihilates it when the background is on
     shell (a necessary condition; it fixes k within this ansatz).  The backgrounds here are on shell at the origin only, so eta is taken to vanish at the origin (then the statement
     needs the field equations at the origin and nothing about their derivatives): eta = one monomial in one component -- the 16
     linear ones and the 16 diagonal quadratic ones x_m^2 e_k (see the review record below for why these and no others).
 (b) THE GHOST OPERATOR, by direct evaluation:  (M eta)_m = C_m[delta h, delta phi]  with the paper's gauge condition (4.27),
     C_m = D^n h_mn - 1/2 D_m h - c phi d_m phi-tilde, on every polynomial eta of degree <= 2 (4 constant, 16 linear, 40 quadratic
     base fields, each separately -- nothing randomized), at the origin of normal coordinates:  M eta = s_gh (D^2 eta + E_gh eta).

PREDICTIONS (rule 25; in this file as first written, before the first run -- the first log's file time is 10 Nov 07:06 ET; the
"07:15" I first typed here, and "07:25" for the first smoke, were guesses, not clock readings; nothing fetched):
  G0 [the symmetry]: on shell, the ungauged operator gives 0 (all 11 components) on every pure-gauge fluctuation tried, with k = 1.
      Controls: off shell it does not (some component nonzero for some eta); on shell with k = 2 (--k 2) it does not.
  G1 [derived by hand, 10 Nov: D^n(D_m eta_n + D_n eta_m) - D_m D.eta = D^2 eta_m + [D_n, D_m] eta^n = D^2 eta_m + Ric_mn eta^n]:
      M has no first-derivative terms at the origin, second-derivative part delta_mn delta^{ab} (s_gh = 1), on and off shell.
  G2 [the content]: E_gh = Ric - c k v v as a 4 x 4 matrix, on and off shell: with c = k = 1 exactly the Ric - v v that
      hv1974_gravscalar_M.py (its ghost line) and rung3_dewitt_bundle.py ('ghostscalar') use.  No W term.
  G3 [control]: with --c 0, E_gh = Ric (the pure-gravity ghost of step 2).
Status PASS if G0 (with its off-shell control), G1, G2 hold on at least one on-shell and one off-shell background, c = k = 1.
What this does NOT derive: that the ghost is a complex (anticommuting) vector entering with the factor -2 -- that is the
Faddeev-Popov determinant itself, taken as in step 2.

Run:  python3 ghost_operator_derive.py --jobs 16 --log results/ghost_operator_derive.log
      controls:  --k 2 --offshell 0   (G0 must FAIL);   --c 0 --no_symmetry   (E_gh = Ric; seconds)

SMOKE RECORD (rule 11).  First smoke (sandbox 10 Nov 07:06, a random 6 of the 136 monomial eta's of degree 1-3 per background): G1,
G2 held; G0 printed PASS but had drawn NO linear eta on the on-shell background -- and off shell only the linear eta gave a nonzero
result (1 of 1; quadratic and cubic 0 of 5).  So that PASS had not tested what G0 is about.  Changed before any further run: all 16
linear eta always included.  Then (07:27, 16 linear + 4 higher): on shell 0 of 20 nonzero, off shell 16 of 16 linear nonzero,
status=PASS; --k 2: on shell 16 of 16 linear and 1 of 3 quadratic nonzero, FAIL as intended; --c 0: E_gh = Ric.
FRESH REVIEW (10 Nov, separate agent; mutations on scratch copies; repository untouched).  No error found in code or physics.
  * The logic of (a), by an independent formula: differentiating the Noether identity S'.delta_eta Phi = 0 gives
    S''[delta_eta Phi-tilde] = -S'.(delta_eta)', which at the origin is EOM(0) x d eta(0) plus eta(0) x d EOM(0); so with eta(0) = 0
    only the field equations AT the origin enter, and the restriction is legitimate.  The reviewer checked more than "nonzero off
    shell": the 11 off-shell components equal (d.eta) E^ab - d_c eta^a E^cb - E^ac d_c eta^b, phi-row (d.eta) tr W, exactly.
  * Which eta test what: linear eta carry the field equations and k; the diagonal quadratic x_m^2 e_k feel k only (phi-row,
    (k - 1) 2 v_k); off-diagonal quadratic and all cubic eta are annihilated identically, on and off shell, for any k -- they test
    nothing here.  Hence the present choice (16 + 16) and no --full.  k is fixed uniquely, sign included ((k - 1) x nonzero).
  * G2's target is Ric - c k v v, so G2 alone never catches a wrong k (--k -1 --no_symmetry prints True): k is guarded by G0 and by
    the c = k = 1 condition for PASS.  G0 also ties the background's Ric = 1/2 v v to the Lagrangian's 1/2.
  * (b): gauge_condition_at_0 (retyped here) agrees with the gauge term inside L2 (the reviewer's check: L2(gauge on) - L2(gauge off)
    at the origin = -1/2 |C(0)|^2 on random eta, c = 1, -1, 0).  Mutations: factor 2 on the connection term dropped, or its sign
    flipped -> G2 False (and G0 nonzero on shell).  Dropping the order-h terms of the transformation is right at one loop.
  * NOT testable here: g^-1 in delta phi = eta^a d_a phi-tilde (replacing it by delta changes no polynomial at the orders used;
    harmless for E(0)); the constant-eta part of the symmetry (would need the field equations to first order in x); G1 away from
    the origin (by covariance).  "Exactly what those scripts use" for Ric - v v is true by READING them (hv1974_gravscalar_M.py's
    ghost line; rung3_dewitt_bundle.py 'ghostscalar'), the derived M = +(D^2 + E) being minus their -(D^2 + E).
  * Folded in: the unused import of the print script removed; no random sampling of eta any more (it shared the backgrounds'
    random stream, so the backgrounds changed with the options); "iff" -> "only if"; the times above corrected.
  * After the fold-in (sandbox 10 Nov 08:07-08:30, this version, --jobs 2): complete run status=PASS, 816 s -- on shell 0 of 32
    eta nonzero, off shell 16 of 16 linear nonzero and 0 of 16 diagonal quadratic; --k 2 --offshell 0: 32 of 32 nonzero, FAIL as
    intended (the diagonal quadratics do feel k, as the review said); --c 0 --no_symmetry: E_gh = Ric.
  * WHAT 3a STANDS ON after this: both operators (print vs action, plus hand derivations) are two-source, and the map E -> divergence
    is two-method (rung 3 note 4d).  Shared by both methods and hand-checked only: the Faddeev-Popov -2 and the complex ghost; the
    two reductions (Gauss-Bonnet; W.W -> (tr W)^2 - Ric.vv), coded once in hv1974_gravscalar_M.reduce; the on-shell substitution
    (one line there; the reviewer redid it by hand: 1/320 + 43/480 + 1/2 + 1/24 = 203/320).
"""
import sys, os, time, argparse, random
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gravscalar_operator_derive as D
from gravscalar_operator_derive import idx, XR, ZERO, C, T, d, psum, at0_rational

out = None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def eta_field(k, mon, coeff=1):
    """the covector field eta_a(x) = coeff * x^mon in component k, zero elsewhere"""
    e = [ZERO for _ in idx]
    e[k] = C(coeff) * XR[0]**mon[0] * XR[1]**mon[1] * XR[2]**mon[2] * XR[3]**mon[3]
    return e


def pure_gauge(bg, eta, k_scalar, xdeg):
    """delta h_ab = D_a eta_b + D_b eta_a (background connection), delta phi = k eta^a d_a phi-tilde, to x^xdeg"""
    h = [[T(d(eta[b], a) + d(eta[a], b) - psum(bg.Gb[c][a][b] * eta[c] for c in idx) * C(2), xdeg) for b in idx] for a in idx]
    ph = T(C(k_scalar) * psum(bg.gbinv[a][b] * eta[b] * bg.dphit[a] for a in idx for b in idx), xdeg)
    return h, ph


def gauge_condition_at_0(bg, h, ph, c_gauge):
    """C_m = g^{nr} D_r h_mn - 1/2 D_m (g^{ab} h_ab) - c phi d_m phi-tilde at the origin (sympy Rationals)"""
    Gb, gbinv = bg.Gb, bg.gbinv
    Dh = [[[T(d(h[a][b], c) - psum(Gb[l][c][a] * h[l][b] + Gb[l][c][b] * h[a][l] for l in idx), 0) for c in idx] for b in idx] for a in idx]
    Dtrh = [T(psum(gbinv[a][b] * Dh[a][b][c] for a in idx for b in idx), 0) for c in idx]
    Cl = [T(psum(gbinv[n][r] * Dh[m][n][r] for n in idx for r in idx) - Dtrh[m] * C(sp.Rational(1, 2)) - C(c_gauge) * ph * bg.dphit[m], 0) for m in idx]
    return [at0_rational(x) for x in Cl]


def ghost_operator(bg, c_gauge, k_scalar):
    """returns s_gh, flags, E_gh (4x4) from all 60 polynomial eta's of degree <= 2"""
    M = lambda eta: sp.Matrix(gauge_condition_at_0(bg, *pure_gauge(bg, eta, k_scalar, 1), c_gauge))
    mono = lambda *ms: tuple(sum(1 for m in ms if m == i) for i in idx)
    const = sp.Matrix(4, 4, lambda i, j: M(eta_field(j, (0, 0, 0, 0)))[i])
    first_zero = all(M(eta_field(k, mono(m))) == sp.zeros(4, 1) for k in idx for m in idx)
    s_set = set(); second_ok = True
    for k in idx:
        for m in idx:
            for n in idx:
                if m > n:
                    continue
                col = M(eta_field(k, mono(m, n), sp.Rational(1, 2) if m == n else 1))
                if m == n:
                    second_ok &= all(col[i] == 0 for i in idx if i != k)
                    s_set.add(col[k])
                else:
                    second_ok &= col == sp.zeros(4, 1)
    second_ok &= len(s_set) == 1
    s = list(s_set)[0] if len(s_set) == 1 else None
    # (D^2 eta)_a(0) on a constant covector: -(d_m Gamma^c_{m a})(0) eta_c
    conn = sp.Matrix(4, 4, lambda a, c: -sum(at0_rational(T(d(bg.Gb[c][m][a], m), 0)) for m in idx))
    E = (const / s - conn) if s not in (None, 0) else None
    return s, second_ok, first_zero, E, sorted(s_set, key=str)


_W = {}


def _init(Rd, v, W, cub, k_scalar):
    _W['bg'] = D.Background(Rd, v, W, cub); _W['probes'] = D.probe_fields(); _W['k'] = k_scalar


def _sym(spec):
    """the ungauged Euler-Lagrange operator at the origin on the pure-gauge fluctuation of eta = x^mon e_k"""
    k, mon = spec
    bg = _W['bg']
    base = pure_gauge(bg, eta_field(k, mon), _W['k'], 2)
    return D.el_at_origin(bg, base, 0, _W['probes'], gauge_sign=0)


def symmetry_test(Rd, v, W, cub, k_scalar, specs, jobs):
    if jobs > 1:
        import multiprocessing as mp
        with mp.Pool(jobs, initializer=_init, initargs=(Rd, v, W, cub, k_scalar)) as pool:
            vals = pool.map(_sym, specs, chunksize=1)
    else:
        _init(Rd, v, W, cub, k_scalar)
        vals = [_sym(s_) for s_ in specs]
    return [sum(1 for x in val if x != 0) for val in vals]


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261110)
    ap.add_argument('--c', type=int, default=1, help='coefficient of phi d phi-tilde in the gauge condition (the paper: 1)')
    ap.add_argument('--k', type=int, default=1, help='coefficient of eta.d phi-tilde in the scalar transformation (1; --k 2 is the control for G0)')
    ap.add_argument('--onshell', type=int, default=1)
    ap.add_argument('--offshell', type=int, default=1)
    ap.add_argument('--no_symmetry', action='store_true', help='skip (a): the ghost operator only (seconds; for the --c 0 control; never status=PASS)')
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; to start over:  del {args.log.replace('/', chr(92))}"); sys.exit(1)
        out = open(args.log, 'w')
    say(f"# ghost_operator_derive: seed {args.seed}, c = {args.c}, k = {args.k}, {args.onshell} on-shell + {args.offshell} off-shell backgrounds, G0 on the 16 linear + 16 diagonal quadratic monomial eta fields, jobs {args.jobs}")
    say(D.version_line())
    rng = random.Random(args.seed)
    unit = lambda m, w: tuple(w if i == m else 0 for i in idx)
    specs = [(k, unit(m, 1)) for k in idx for m in idx] + [(k, unit(m, 2)) for k in idx for m in idx]
    ok = True; g12 = True; kinds_done = set()
    g0_on = g0_off_control = None
    for kind, n in (('on shell', args.onshell), ('off shell', args.offshell)):
        for i in range(n):
            R = D.riemann_sample(rng); v, W, cub = D.scalar_sample(rng)
            if kind == 'on shell':
                W = W - W.trace() * sp.eye(4) / 4
                R = D.with_ricci(R, v * v.T / 2)
            Ric, Rs = D.ricci_of(R)
            onsh = (Ric - v * v.T / 2) == sp.zeros(4, 4) and W.trace() == 0
            bg = D.Background(R, v, W, cub)
            # (b) the ghost operator
            s, second_ok, first_zero, E, sset = ghost_operator(bg, args.c, args.k)
            if E is None:
                say(f"{kind} {i}: FAIL: second-derivative part does not give one normalisation: {sset}"); say("\nstatus=FAIL"); sys.exit(1)
            target = Ric - args.c * args.k * v * v.T
            g2 = E.applyfunc(sp.expand) == target.applyfunc(sp.expand)
            g1 = s == 1 and second_ok and first_zero
            say(f"{kind} {i}: R = {Rs}, v^2 = {(v.T * v)[0]}, tr W = {W.trace()}, field equations hold at 0: {onsh}")
            say(f"   ghost: s_gh = {s}, second-derivative part delta: {second_ok}, no first-derivative terms: {first_zero} (all 60 base fields);  E_gh == Ric - {args.c * args.k} v v: {g2}"
                + ("" if g2 else f";  E_gh - Ric = {(E - Ric).applyfunc(sp.expand).tolist()},  v v = {(v * v.T).tolist()}") + f"  [{time.time() - t0:.0f} s]")
            if args.no_symmetry:
                ok &= g1 and g2; g12 &= g1 and g2; kinds_done.add(kind); continue
            # (a) the symmetry
            nz = symmetry_test(R, v, W, cub, args.k, specs, args.jobs)
            nbad = sum(1 for x in nz if x)
            by_deg = {dg: (sum(1 for sp_, x in zip(specs, nz) if sum(sp_[1]) == dg and x), sum(1 for sp_ in specs if sum(sp_[1]) == dg)) for dg in (1, 2)}
            say(f"   symmetry (ungauged operator on the pure-gauge fluctuation, k = {args.k}): nonzero for {nbad} of {len(specs)} eta fields  [by degree of eta, nonzero/tried: " + ', '.join(f"{dg}: {a}/{b}" for dg, (a, b) in by_deg.items()) + f"]  [{time.time() - t0:.0f} s]")
            if kind == 'on shell':
                g0_on = (g0_on if g0_on is not None else True) and nbad == 0 and onsh
            else:
                g0_off_control = (g0_off_control or False) or nbad > 0
            ok &= g1 and g2; g12 &= g1 and g2
            kinds_done.add(kind)
    say("")
    if g0_on is not None:
        say(f"G0 on shell the pure-gauge fluctuations are annihilated (k = {args.k}): {'PASS' if g0_on else 'FAIL'}")
        ok &= g0_on
    if g0_off_control is not None:
        say(f"G0 control, off shell they are not: {'as predicted' if g0_off_control else 'NOT as predicted (all zero off shell)'}")
        ok &= g0_off_control
    say(f"G1, G2 (c = {args.c}, k = {args.k}): {'PASS' if g12 else 'FAIL'}")
    complete = kinds_done == {'on shell', 'off shell'} and args.c == 1 and args.k == 1 and not args.no_symmetry
    if not ok:
        say("\nstatus=FAIL")
    elif not complete:
        say("\nstatus=INCOMPLETE (PASS needs c = k = 1, the symmetry test, and at least one on-shell and one off-shell background)")
    else:
        say("\nstatus=PASS")
    say(f"END ghost_operator_derive {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
