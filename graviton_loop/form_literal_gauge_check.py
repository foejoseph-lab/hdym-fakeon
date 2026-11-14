#!/usr/bin/env python3
"""
Rung 3, step 3b, section 7: THE GAUGE WITHOUT K ((4.27) exactly as printed, on the class (7.10)), on shell.
PREDICTIONS WRITTEN 13 Nov 2026 20:25 ET (clock read by the command that wrote them), before any code for this existed and
before FORM's literal output or the jet probes' had been compared with anything.  Shown to the author at once; he started the
jet-probe run for it (20:40 ET).  They stand as first written (rule 11).  The reader is this file.

THE QUESTION (item 1 of the author's order, 13 Nov 20:24 ET).  The gauge condition C_m = D^n h_mn - 1/2 D_m h - c phi d_m phi-bg,
with no K, is the other reading of what section 7 could have used.  Is the on-shell divergence still (203/320) K^2 ((d phi)^2)^2
up to total derivatives?  At K' = 0 this is the c-family already tested.  At K' != 0 it is new.

THE ARGUMENT THE PREDICTIONS COME FROM (hand argument; nothing run).
  1. In the canonical field chi (d chi = sqrt(K) d phi, fluctuation psi = sqrt(K) phi) the literal gauge reads
         C_m = D^n h_mn - 1/2 D_m h - gamma(chi-bg) psi d_m chi-bg,      gamma = c / K.
     So it is step 3's c-family with the coefficient a FUNCTION of the background field.  gamma' = d gamma / d chi = - c K1 / K^(5/2).
  2. The operator contains gamma and gamma' and no gamma'' (the gauge term is quadratic in C; one derivative can fall on gamma).
     So the total is a polynomial in gamma, gamma' and never contains gamma''.  In the original variables: K2 does not enter.
  3. On shell at a point there are four structures (as in R8): (v^2)^2, v.Wc.v, Wc.Wc + 1/2 (v^2)^2, Gauss-Bonnet.  The total
     derivatives with a function of chi in front are, on shell,
         D_a(F v^a v^2)                   = F' (v^2)^2 + 2 F v.Wc.v
         D_a(H (Wc^ab v_b - v^a Box chi)) = H' v.Wc.v + H (Wc.Wc + 1/2 (v^2)^2)
     (a third, D_a(F G^ab v_b), is one quarter of the first on shell.)
  4. At constant gamma the author's runs and FORM gave  b = 0  and  C = (gamma^2 + 10 gamma + 1)/12  at six values.  So
     H = (gamma^2 + 10 gamma + 1)/12 and F = 0 as functions of gamma.  If F or H also depended on gamma', their derivative would
     bring gamma'' into the total, which point 2 forbids.  Hence, if the divergence is gauge independent up to total derivatives:

PREDICTIONS (rule 25).  Fit basis and units as in R8 (v = sqrt(K) u, Wc = sqrt(K) Wt; the reader's 'vWv' row is K u.Wt.u).
  L1  [two methods on the new operator]  FORM's literal operator equals the jet probes' (class710_operator_derive.py --gauge
      literal, from its checkpoint) entry by entry: second-derivative part, first-derivative coefficient, E, Tr Omega^2, ghost
      (the ghost is Box eta + (Ric - c u u) eta: no K).
  L2  The four structures suffice: zero residual on the backgrounds beyond the four fitted, for every K jet set and c tried.
  L3  [THE NUMBER -- can fail]  alpha = 203/320, for every (K, K1, K2, c): no dependence on gamma or on gamma'.
      In particular NO gamma'^2 term.
  L4  C = (gamma^2 + 10 gamma + 1)/12 with gamma = c/K at the point, whatever K1 is.
  L5  [new, sharp]  b = gamma' (gamma + 5)/6 per v.Wc.v;  in the reader's units  b = - (c k1 / K)(c/K + 5)/6,  k1 = K1/K.
      Linear in K1; zero when K1 = 0 or c = 0.
  L6  g = 7/12.
  L7  K2 enters none of alpha, b, C, g.
  L8  [control: the routine on the known case]  the same fit on FORM's image-gauge operator returns alpha = 203/320, b = 0,
      C = (c^2 + 10 c + 1)/12, g = 7/12 at the same K jets and c.
  Cases: K jets of a = 1/12 at phi0 = 1/2 and at two other phi0; a = 1/6; free jets with K1 = 0 and with K1 != 0, K2 varied;
  c in {1, 2, 1/2, -1}.  Six on-shell backgrounds per case (four fitted, two for the residual).  FORM: seconds.

WHAT A FAILURE WOULD MEAN.  L3 failing with L2 holding is either a gauge dependence of the on-shell divergence -- which the
paper's ref. 7 says cannot be -- or an error in point 2 or 3 above (a total derivative I have not listed, or a gamma'' I say is
absent).  A full stop on the author's side until which, and by how much.  L5 failing with L3 holding means F is not zero, and the
difference from the prediction is twice F.

WHAT A PASS WOULD AND WOULD NOT SHOW.
  Would: that both readings of "the gauge of the previous sections" give the same on-shell number for the section 7 class; and
  a gauge with a field-dependent coefficient, where the invariance holds only after total derivatives with non-constant
  coefficients are removed -- a stronger test than a constant c.
  Would not: any other gauge family; the off-shell result in either gauge; anything shared downstream (hv340, the ghost's -2,
  the reductions, the on-shell substitution, the count of four structures); the Jacobian.

RUNS.  One jet-probe run on the author's machine supplies L1's data (existing, reviewed code; about 25 minutes by the 13 Nov
timings; resumable; it ends status=INCOMPLETE by design -- it tests no prediction itself, the reader does):
      python3 class710_operator_derive.py --jobs 16 --gauge literal --c 1 --cases 1/12 --onshell 5 --offshell 1 --mode min --crosscheck 2 --log results/class710_literal_c1.log
Then FORM's literal output (already in results/ from 13 Nov) and the new reader, about a minute.
WHAT HAPPENED (sandbox, 13 Nov, FORM's outputs as pushed by the author 20:07 ET; this block written 13 Nov 20:53 ET, clock read by the
command that wrote it).  L2-L8 held on the first run, 71 s: seven sets of K jets x four values of c, six on-shell backgrounds
each.  alpha = 203/320 in all 28; b equal to the prediction in all 28 (e.g. 68/243 at K jets (3/2, -2/3, 1), c = 1; zero at
K1 = 0); C and g as predicted; residuals zero; the image-gauge control as it must be.  L1 NOT TESTED in the sandbox: it needs
the author's jet-probe run.  Until L1 the result stands on FORM's literal operator alone.

FRESH REVIEW (13 Nov, separate agent, no context; folded in before the patch).  The comparison is real and not circular; the
units of b, both total-derivative identities and the list of on-shell currents were re-derived by hand and hold; six mutations
(of the prediction, of three terms of FORM's scalar equation, of the conversion, of the ghost; the two FORM operators swapped)
were all caught.  L3 AND L5 CARRY THE TEST.  Corrections to the header above:
  * L2 is close to automatic: on shell, without derivatives of the curvature, there are no other dimension-four scalars, so
    any covariant operator gives zero residuals.  It guards covariance and the evaluation, nothing more.  No mutation moved it.
  * L7 is not a separate prediction: the predicted coefficients have no K2, so it follows from L3-L6 at the one pair of cases
    that differ in K2 (K2 = 1 and -4).  The content is that the literal operator does contain K2 and the fits do not.
  * "L5 failing with L3 holding means F is not zero" is wrong by the header's own point 4: that outcome would contradict the
    argument, not measure F.
  * Point 4 in full: alpha = 203/320 + F', b = 2F + H', C = H.  No gamma'' in alpha makes F a function of gamma alone; then
    no gamma'' in b makes H one; constant gamma (the author's runs) gives 2F = 0.  It also uses that F and H depend on chi only
    through gamma (the canonical theory has no other function of chi) -- not said above.
  * The total-derivative argument uses the field equations in a neighbourhood, which is right for the integral.
  Fixed in the code: a resumed jet-probe log reprints its totals and would have failed L1 spuriously (only the part after the
  last RESUMED line is read now); two crash paths (a degenerate fit; a log with no status line).
  Not checked by the reviewer: L1; the FORM program's literal branch beyond that it differs from the image one by c -> c K and
  the c K1 terms; point 2 (no gamma'') beyond the K2 pair.

Run (from graviton_loop; after the jet-probe run above; FORM's literal output is already in results/):
      python3 form_literal_gauge_check.py --log results/form_literal_gauge_check.log
About a minute.
"""
import sys, os, time, random, argparse, platform
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import form_class710_check as FC
Z = FC.Z

out = None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def predicted(Kj, c, on, g):
    """alpha, b, C, g in the fit's units, from the header's argument: gamma = c/K, b = gamma'(gamma + 5)/6 with
    gamma' = d gamma / d chi, i.e. -(c k1 / K)(gamma + 5)/6 per K u.Wt.u"""
    K0, K1, K2 = Kj
    gam = c / K0
    return {'v4': on, 'vWv': -(c * K1 / K0**2) * (gam + 5) / 6, 'ricid': (gam**2 + 10 * gam + 1) / 12, 'gb': g}


def fit_case(F, Kj, cvals, rng, nbg=6):
    """the on-shell four-structure fit of FORM's operator F at K jets Kj, for each c: {c: (coefficients, residuals)}"""
    bgs = []
    for i in range(nbg):
        R, u, W, cub = FC.sample(rng, Kj, True)
        T = FC.r3_prediction(R, u, W, Kj); T['uWu'] = (u.T * T['Wt'] * u)[0]
        assert T['onsh']
        bgs.append((R, u, F.pieces(R, u, W, Kj), T))
    res = {}
    for cv in cvals:
        rows = [FC.r8_row(FC.gilkey(P, R, u, Kj, cv)['tot'], T) for R, u, P, T in bgs]
        res[cv] = Z.fit_r8(rows)
    return res


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261115)
    ap.add_argument('--resdir', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results'))
    ap.add_argument('--probe-log', default='class710_literal_c1.log')
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    args.backgrounds = 3
    t0 = time.time()
    if args.log:
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; to start over:  del {args.log.replace('/', chr(92))}"); sys.exit(1)
        out = open(args.log, 'w')
    FC.out = out
    say(f"# form_literal_gauge_check: seed {args.seed}, reading {args.resdir}")
    say(f"# versions: python {platform.python_version()} ({platform.system()} {platform.machine()}), sympy {sp.__version__}; {sys.executable}")
    for tag in ('literal', 'main'):
        fn = os.path.join(args.resdir, f"form_class710_{tag}_Eh.out")
        if not os.path.exists(fn):
            say(f"no FORM output for TAG {tag}: run form/class710_operator.frm (its header has the lines)"); say("status=INCOMPLETE"); sys.exit(1)
        say(f"#   {os.path.basename(fn)}: {os.path.getsize(fn)} bytes, written {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(fn)))}")
    struct = Z.step3_structure(args.seed + 7)
    on, g712 = struct[2], struct[0]['Riem2']
    FL = FC.FormOperator(args.resdir, 'literal'); FI = FC.FormOperator(args.resdir, 'main')
    if not (FL.readable and FI.readable):
        say("FORM's output is not readable as an operator"); say("status=FAIL"); sys.exit(1)
    ok = True
    # ---- L1
    say("## L1  FORM's literal operator against the jet probes' (--gauge literal), from the run's checkpoint")
    coll = []
    ok1, nseen = FC.check_logs(FL, args, args.resdir, lognames=(args.probe_log,), collect=coll)
    say(f"L1 {'PASS' if ok1 and nseen else 'FAIL' if nseen else 'NOT TESTED'} ({nseen} backgrounds)   [{time.time() - t0:.0f} s]")
    ok &= ok1 or nseen == 0
    if nseen:
        ons = [x for x in coll if x[0] == 'on shell']
        if len(ons) >= 5:
            rows = []
            for kind, R, u, W, cub, Kj, cval, G in ons:
                T = FC.r3_prediction(R, u, W, Kj); T['uWu'] = (u.T * T['Wt'] * u)[0]
                rows.append(FC.r8_row(G['tot'], T))
            cf, resid = Z.fit_r8(rows)
            pr = predicted(ons[0][5], ons[0][6], on, g712)
            good = cf is not None and all(cf[k] == pr[k] for k in pr) and all(x == 0 for x in resid)
            if cf is None:
                say(f"   on the run's own {len(ons)} on-shell backgrounds: the four structures are degenerate; no fit: FAIL")
            else:
                say(f"   on the run's own {len(ons)} on-shell backgrounds (totals identical by both methods): alpha = {cf['v4']}, b = {cf['vWv']}, C = {cf['ricid']}, g = {cf['gb']}, residuals {resid};  predicted {pr}: {'PASS' if good else 'FAIL'}")
            ok &= good
    # ---- L2-L8
    rng = random.Random(args.seed)
    half = sp.Rational(1, 2)
    cases = [("a = 1/12 at phi0 = 1/2", Z.K_jets_of_a(sp.Rational(1, 12), half)),
             ("a = 1/12 at phi0 = 1", Z.K_jets_of_a(sp.Rational(1, 12), sp.Integer(1))),
             ("a = 1/12 at phi0 = 3/2", Z.K_jets_of_a(sp.Rational(1, 12), sp.Rational(3, 2))),
             ("a = 1/6 at phi0 = 1/2", Z.K_jets_of_a(sp.Rational(1, 6), half)),
             ("free jets, K1 = 0", [sp.Rational(7, 3), sp.Integer(0), sp.Rational(5, 2)]),
             ("free jets (3/2, -2/3, 1)", [sp.Rational(3, 2), sp.Rational(-2, 3), sp.Integer(1)]),
             ("the same K, K1 with K2 = -4", [sp.Rational(3, 2), sp.Rational(-2, 3), sp.Integer(-4)])]
    cvals = [sp.Integer(1), sp.Integer(2), half, sp.Integer(-1)]
    L = dict(L2=True, L3=True, L4=True, L5=True, L6=True, L8=True)
    table = {}
    for name, Kj in cases:
        lit = fit_case(FL, Kj, cvals, rng); img = fit_case(FI, Kj, cvals, rng)
        say(f"## {name}: K jets ({Kj[0]}, {Kj[1]}, {Kj[2]})   [{time.time() - t0:.0f} s]")
        for cv in cvals:
            cf, resid = lit[cv]; pr = predicted(Kj, cv, on, g712)
            if cf is None:
                say(f"   c = {cv}: degenerate fit"); L['L2'] = False; continue
            l2 = all(x == 0 for x in resid); l3 = cf['v4'] == pr['v4']; l4 = cf['ricid'] == pr['ricid']; l5 = cf['vWv'] == pr['vWv']; l6 = cf['gb'] == pr['gb']
            ci, ri = img[cv]
            l8 = ci is not None and ci['v4'] == on and ci['vWv'] == 0 and ci['ricid'] == (cv**2 + 10 * cv + 1) / 12 and ci['gb'] == g712 and all(x == 0 for x in ri)
            say(f"   c = {cv} (gamma = {cv / Kj[0]}): alpha = {cf['v4']} [{'ok' if l3 else 'NOT ' + str(pr['v4'])}];  b = {cf['vWv']} [{'ok' if l5 else 'predicted ' + str(pr['vWv'])}];  C = {cf['ricid']} [{'ok' if l4 else 'predicted ' + str(pr['ricid'])}];  g = {cf['gb']} [{'ok' if l6 else 'NOT ' + str(pr['gb'])}];  residuals {resid};  image-gauge control {'as it must be' if l8 else 'WRONG: ' + str(ci)}")
            for k, v_ in (('L2', l2), ('L3', l3), ('L4', l4), ('L5', l5), ('L6', l6), ('L8', l8)):
                L[k] &= v_
            table[(name, cv)] = cf
    # L7: the last two cases differ only in K2 (and in their random backgrounds)
    l7 = all(table.get((cases[5][0], cv)) == table.get((cases[6][0], cv)) and table.get((cases[5][0], cv)) is not None for cv in cvals)
    L['L7'] = l7
    say("")
    say("   ".join(f"{k} {'PASS' if L[k] else 'FAIL'}" for k in ('L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8')))
    if not L['L3'] and L['L2']:
        say("L3 FAILS WITH THE FIT CLOSED: either the on-shell divergence depends on this gauge, or the header's argument (points 2, 3) has a hole.  A full stop until which, and by how much.")
    ok &= all(L.values())
    complete = nseen >= 6
    say("status=" + ("FAIL" if not ok else "PASS" if complete else "INCOMPLETE (PASS needs L1: the jet-probe run class710_literal_c1.log with its checkpoint)"))
    say(f"END form_literal_gauge_check {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
