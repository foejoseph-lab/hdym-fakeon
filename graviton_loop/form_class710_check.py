#!/usr/bin/env python3
"""
Rung 3, step 3b, second method for route 2 (FORM), stage F1: the reader for the class (7.10).

form/class710_operator.frm derives, symbolically and covariantly, the Euler-Lagrange expressions of the gauge-fixed second-order
Lagrangian of sqrt(g)(R - 1/2 K(phi)(d phi)^2) with the gauge coefficient c a symbol, the ghost operator, and the covariant
derivative D_n B^m of the first-derivative coefficient.  Its header carries the predictions F1.1-F1.10 (written 13 Nov 17:52 ET,
before any F1 code and before any of the author's R8 or main-run logs were seen; the author's yes 19:36 ET, after those runs).
This script evaluates FORM's monomials on backgrounds in exact rationals and does the 11 x 11 matrix algebra:
    S = diag(P, K) (tested, F1.1);  omega^m = 1/2 S^-1 B^m;  D_n omega^m = 1/2 S^-1 D_n B^m - 1/2 S^-1 (D_n S) S^-1 B^m;
    E = S^-1 C - D_m omega^m - omega.omega;  Omega_mn = [D_m, D_n] on h + D_m omega_n - D_n omega_m + [omega_m, omega_n];
    then (3.40) by hv1974_gravscalar_M.hv340, the ghost, and class710_operator_derive.fit_r8.
Typed here, not derived by FORM: D_n S = K1 v_n on the scalar entry, and [D_m, D_n] h_ab = R(a,e,m,n) h_eb + R(b,e,m,n) h_ae
(the convention of the FORM program).

  F1.1-F1.6  on fresh backgrounds (random K jets and a = 1/12; three off shell and three on shell each): the second-derivative
        part; the first-derivative part at c = 1; E against the prediction R3 of class710_operator_derive.py and Tr Omega^2
        against step 3's; K2 dropping out; the ghost; the total against R6 and, on shell, R7.
  F1.7  the R8 fit from FORM's operator at c = -1, 0, 1/2, 1, 2, 3 on six on-shell backgrounds.
  F1.8  FORM against the jet probes on the author's three runs: the backgrounds are regenerated from each log's seed by
        class710_operator_derive's own sampling, the probes' coefficient data are read from the run's checkpoint (read only;
        the key carries the code version printed in that log), and S, B^m, E (11 x 11), Tr Omega^2 are compared entry by
        entry, the ghost with ghost_operator_derive's on the same background; at c != 1 also the pointwise total printed in
        the log.  Nothing is typed, and no base field is recomputed by probes.
  F1.9  controls: D_n B^m left out of the conversion; FORM run without K in the gauge term; FORM run with the commutator reversed.

WHAT EACH CHECK WOULD AND WOULD NOT CATCH.
  F1.8 is the content: two derivations of different kinds, number against number, including c = 0 and c = 2 where the
  jet-probe run had no other test of E, omega, Omega than R8 itself.  It also tests what is typed here (D_n S, the curvature on
  h, the sign of the cross term between that curvature and D omega in Tr Omega^2, which is invisible at c = 1).
  F1.3 compares with the same target the jet probes were compared with; F1.6 then follows from F1.3 and F1.5.
  F1.7 rests on FORM's operator (tested at c = 0, 1, 2 by F1.8, at other c by nothing but the polynomial dependence on c) and on
  alpha(c) being a polynomial of degree at most four: E is quadratic in c, the total quadratic in E and in the ghost's E.  Five
  equal values would already force a quartic to be constant; there are six.
  NOT caught: anything downstream that both methods share (hv340, the ghost's -2, the reductions, the on-shell substitution, the
  R8 argument that alpha is the invariant); the literal gauge (run and reported as a control at c = 1 only); the Jacobian.
SANDBOX RECORD (13 Nov evening, FORM 4.3.1, python 3.13.18, sympy 1.13.1; this block written 13 Nov 19:57 ET, clock read by the command
that wrote it; the logs and checkpoints are the author's, pushed 19:31 ET, code version e631f8b77425).
  FORM 0.00 s, 39 + 13 + 7 + 8 + 10 terms.  This reader, about 30 s.  On the first run: F1.1-F1.6 hold on twelve fresh
  backgrounds.  F1.7: alpha = 203/320, b = 0, g = 7/12, residuals zero at c = -1, 0, 1/2, 1, 2, 3; for the record
  C = -2/3, 1/12, 25/48, 1, 25/12, 10/3, which is (c^2 + 10 c + 1)/12 at all six.  F1.8: identical on all 22 of the author's
  backgrounds (5 at c = 0, 5 at c = 2, 12 at c = 1).  F1.9: each control breaks what it must.
  F1.10, by mutation: none of the four things F0 could not test becomes testable here.  Reversing the order of the two by-parts
  derivatives, or the index order in DT's first term, changes FORM's written text and not one number of the operator; dropping
  the second-order term of GG, or the Leibniz lines for H3, DW, DRic, leaves FORM's output identical.  They remain right by
  reading only.
  Nothing here has run on the author's machine.

FRESH REVIEW (13 Nov evening, separate agent, no context, scratch copies; folded in before the patch).
  F1.8 is real: the two sides are fed from disjoint sources (FORM's files; the checkpoint), no base field is recomputed (a spy
  on the probe task saw none), the backgrounds are regenerated exactly as class710_operator_derive.main draws them (all 66 keys
  found on each of the 22), and changing one rational in one checkpoint line (-512/2209 to -513/2209) made exactly that
  background differ.  gilkey() and the assembly of D_n B^m are right by reading and by mutation; the degree bound holds.
  Its 13 mutations of the new parts were all caught.  FIVE ARE INVISIBLE TO EVERYTHING AT c = 1, the 12-background c = 1 log
  included, and are caught only at c != 1: the K rule and the jm/jn order in FORM's last module, the sign of the curvature on h,
  the order of the two D omega terms in Omega, and the m/n order of D_n B^m in this reader.  That is what the c = 0 and c = 2
  runs and their checkpoints bought.  F1.7 alone would not have been independent evidence for them.
  Folded in: (1) the S compared in F1.8 was the typed diag(P, K), not FORM's -- FORM's second-derivative part is now
  evaluated on every one of the 22 backgrounds too; (2) a missing pointwise total in a log would have gone uncompared while the
  summary said "identical" -- now the count must match; (3) the ghost entered F1.8 only through the total at c = 2 -- it is now
  compared on all 22 backgrounds with the probes' own ghost (ghost_operator_derive.ghost_operator, 0.3 s each); (4) the F1.7
  line printed its conclusion even on FAIL; (5) FORM's last module does not differentiate curvature -- the reader now refuses
  an output whose first-derivative coefficients contain R, Ric or W.
  Too strong above, corrected: the degree count should say that Omega is quadratic in c through [omega, omega], so Tr Omega^2
  is quartic like Tr E^2 (the bound of four stands); the flip control breaks the ghost and the total as well as F1.3.
  Not checked by the reviewer: more values of c; the inherited F0 code; hv340, fit_r8, step3_structure, r3_prediction; a
  Windows run.

Run (from graviton_loop, after the three FORM runs in WSL listed in form/class710_operator.frm):
      python3 form_class710_check.py --log results/form_class710_check.log
About a minute.
"""
import sys, os, re, time, random, argparse, platform, hashlib, json
from fractions import Fraction
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gravscalar_operator_derive as D
import hv1974_gravscalar_M as GM
import form_operator_check as FO
import class710_operator_derive as Z
import ghost_operator_derive as GH

idx = range(4)
out = None
SYMS = ('c', 'p', 'pp', 'K', 'K1', 'K2', 'K3')
KNOWN = {'h', 'H1', 'S2', 'p1', 'pS2', 'eta', 'SE2', 'R', 'Ric', 'v', 'W', 'd_', 'hh'}
ESS = sp.zeros(11, 11); ESS[10, 10] = 1


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def rat(x):
    return sp.Rational(x.numerator, x.denominator) if x else sp.Integer(0)


def parse(text):
    """FORM expression -> list of (rational coefficient, {symbol: power}, [(tensor name, (index names))])"""
    s = re.sub(r'\s+', '', text).replace('_?', '')
    terms = []
    for t in re.findall(r'[+-]?[^+-]+', s):
        coef = Fraction(-1 if t[0] == '-' else 1); t = t.lstrip('+-')
        syms = {}; tens = []
        for f in re.findall(r'[A-Za-z_0-9]+\([^)]*\)|[^*]+', t):
            m = re.fullmatch(r'([A-Za-z][A-Za-z0-9]*)(\^(\d+))?', f)
            if '(' in f:
                name, args = f[:-1].split('(')
                tens.append((name, tuple(args.split(','))))
            elif re.fullmatch(r'\d+(/\d+)?', f):
                coef *= Fraction(f)
            elif m and m.group(1) in SYMS:
                syms[m.group(1)] = syms.get(m.group(1), 0) + int(m.group(3) or 1)
            else:
                raise ValueError(f"unparsed factor {f!r} in term {t!r}")
        terms.append((coef, syms, tens))
    return terms


def evaluate(terms, tensors, symvals, free):
    """-> {power of c: {free index tuple: Fraction}}.  A tensor given as None, or a symbol given as 0, switches its terms off;
    a tensor or symbol that is not in the tables at all is an error."""
    res = {}
    for coef, syms, tens in terms:
        if any(tensors[name] is None for name, _ in tens):
            continue
        fac = coef; cp = 0; dead = False
        for s, k in syms.items():
            if s == 'c':
                cp = k
            else:
                val = symvals[s]
                if val == 0:
                    dead = True; break
                fac *= Fraction(val) ** k
        if dead:
            continue
        val = FO.contract([(labs, tensors[name]) for name, labs in tens], free)
        acc = res.setdefault(cp, {})
        for k, x in val.items():
            acc[k] = acc.get(k, 0) + fac * x
    return res


class FormOperator:
    """the five expressions one run of form/class710_operator.frm wrote, and their evaluation on a background"""
    def __init__(self, resdir, tag):
        rd = lambda nm: parse(open(os.path.join(resdir, f"form_class710_{tag}_{nm}.out")).read())
        self.Eh, self.Ep, self.Gh, self.DBh, self.DBp = rd('Eh'), rd('Ep'), rd('Gh'), rd('DBh'), rd('DBp')
        allt = self.Eh + self.Ep + self.Gh + self.DBh + self.DBp
        self.names = sorted({n for t in allt for n, _ in t[2]})
        self.syms = sorted({s for t in allt for s in t[1]})
        fl = lambda t, objs, syms: sum(1 for n, _ in t[2] if n in objs) + sum(t[1].get(s, 0) for s in syms)
        self.linear = (all(fl(t, {'h', 'H1', 'S2', 'p1', 'pS2'}, ('p',)) == 1 for t in self.Eh + self.Ep)
                       and all(fl(t, {'eta', 'SE2'}, ()) == 1 for t in self.Gh) and all(fl(t, {'hh'}, ('pp',)) == 1 for t in self.DBh + self.DBp))
        self.unknown = sorted(set(self.names) - KNOWN)
        self.first_clean = not any({'H1', 'p1'} & {n for n, _ in t[2]} and {'R', 'Ric', 'W'} & {n for n, _ in t[2]} for t in self.Eh + self.Ep)
        self.readable = not self.unknown and self.linear and 'K3' not in self.syms and self.first_clean
        self.has_K2 = 'K2' in self.syms
        self.nterms = (len(self.Eh), len(self.Ep), len(self.Gh), len(self.DBh), len(self.DBp))

    def pieces(self, Rd, v, W, Kj, rng=None):
        """S-test, B^m, C, D_n B^m as {power of c: 11 x 11 sympy matrix}; the ghost's matrices; everything exact"""
        Ric, Rs, inv = GM.invariants(Rd, v, W)
        bg = {'R': {k: FO.fr(x) for k, x in Rd.items() if x != 0}, 'Ric': FO.sparse(Ric.tolist(), 2), 'v': FO.sparse(list(v), 1), 'W': FO.sparse(W.tolist(), 2),
              'd_': {(i, i): Fraction(1) for i in idx}}
        off = {'h': None, 'H1': None, 'S2': None, 'p1': None, 'pS2': None, 'eta': None, 'SE2': None, 'hh': None}
        K0, K1, K2 = [FO.fr(x) for x in Kj]
        sv0 = {'K': K0, 'K1': K1, 'K2': K2, 'p': 0, 'pp': 0}

        def O(jets, p=0):
            T = dict(bg); T.update(off); T.update(jets); sv = dict(sv0); sv['p'] = p
            eh = evaluate(self.Eh, T, sv, ('ia', 'ib')); ep = evaluate(self.Ep, T, sv, ())
            res = {}
            for cp in set(eh) | set(ep):
                a = eh.get(cp, {})
                res[cp] = sp.Matrix([rat((a.get((i, j), 0) + a.get((j, i), 0)) / 2) for i, j in D.PAIRS] + [rat(ep.get(cp, {}).get((), 0))])
            return res

        def cols_to_mats(colres):
            cps = set()
            for r in colres:
                cps |= set(r)
            z = sp.zeros(11, 1)
            return {cp: sp.Matrix.hstack(*[r.get(cp, z) for r in colres]) for cp in cps}
        basis = []
        for k in range(11):
            if D.COMPS[k] == 's':
                basis.append(None)
            else:
                a, b = D.COMPS[k]
                basis.append({(a, b): Fraction(1), (b, a): Fraction(1)})
        Cm = cols_to_mats([O({}, 1) if basis[k] is None else O({'h': basis[k]}) for k in range(11)])
        Bm = []
        for m in idx:
            colres = []
            for k in range(11):
                if basis[k] is None:
                    colres.append(O({'p1': {(m,): Fraction(1)}}))
                else:
                    colres.append(O({'H1': {key + (m,): x for key, x in basis[k].items()}}))
            Bm.append(cols_to_mats(colres))
        # D_n B^m
        DB = [[{} for m in idx] for n in idx]
        percol = []
        for k in range(11):
            T = dict(bg); T.update(off); sv = dict(sv0)
            if basis[k] is None:
                sv['pp'] = 1
            else:
                T['hh'] = basis[k]
            dh_ = evaluate(self.DBh, T, sv, ('ia', 'ib', 'jm', 'jn')); dp_ = evaluate(self.DBp, T, sv, ('jm', 'jn'))
            percol.append((dh_, dp_))
        cps = set()
        for dh_, dp_ in percol:
            cps |= set(dh_) | set(dp_)
        for n in idx:
            for m in idx:
                for cp in cps:
                    cols = []
                    for dh_, dp_ in percol:
                        a = dh_.get(cp, {})
                        cols.append(sp.Matrix([rat((a.get((i, j, m, n), 0) + a.get((j, i, m, n), 0)) / 2) for i, j in D.PAIRS] + [rat(dp_.get(cp, {}).get((m, n), 0))]))
                    DB[n][m][cp] = sp.Matrix.hstack(*cols)
        # second-derivative part on random symmetrised second derivatives: must be P Box h, K Box p
        s_ok = None
        if rng is not None:
            rq = lambda: Fraction(rng.choice([k for k in range(-9, 10) if k]), rng.randint(1, 5))
            S2 = {}; pS2 = {}
            for a in idx:
                for b in idx:
                    for m in idx:
                        for n in idx:
                            if a <= b and m <= n:
                                x = rq()
                                for key in ((a, b, m, n), (b, a, m, n), (a, b, n, m), (b, a, n, m)):
                                    S2[key] = x
                    if a <= b:
                        x = rq(); pS2[(a, b)] = x; pS2[(b, a)] = x
            boxh = sp.Matrix(4, 4, lambda a, b: sum(rat(S2[(a, b, m, m)]) for m in idx))
            Ph = boxh / 2 - sp.eye(4) * boxh.trace() / 4
            want = sp.Matrix([Ph[a, b] for a, b in D.PAIRS] + [sp.Rational(Kj[0]) * sum(rat(pS2[(m, m)]) for m in idx)])
            r = O({'S2': S2, 'pS2': pS2})
            s_ok = r.get(0, sp.zeros(11, 1)) == want and all(r[k] == sp.zeros(11, 1) for k in r if k)
        # ghost
        G0 = sp.zeros(4, 4); G1 = sp.zeros(4, 4); ghigh = False
        for k in idx:
            T = dict(bg); T.update(off); T['eta'] = {(k,): Fraction(1)}
            g = evaluate(self.Gh, T, sv0, ('ia',))
            for cp, dct in g.items():
                for (m,), x in dct.items():
                    if cp == 0:
                        G0[m, k] = rat(x)
                    elif cp == 1:
                        G1[m, k] = rat(x)
                    elif x:
                        ghigh = True
        gbox = None
        if rng is not None:
            SE2 = {}
            for a in idx:
                for m in idx:
                    for n in idx:
                        if m <= n:
                            x = rq(); SE2[(a, m, n)] = x; SE2[(a, n, m)] = x
            T = dict(bg); T.update(off); T['SE2'] = SE2
            g2 = evaluate(self.Gh, T, sv0, ('ia',))
            gbox = all(g2.get(0, {}).get((m,), 0) == sum(SE2[(m, n, n)] for n in idx) for m in idx) and all(not any(g2[k].values()) for k in g2 if k)
        return dict(C=Cm, B=Bm, DB=DB, s_ok=s_ok, G0=G0, G1=G1, ghigh=ghigh, gbox=gbox, Ric=Ric, Rs=Rs, inv=inv)


PBIG = sp.diag(D.op_matrix10(lambda h: h / 2 - sp.eye(4) * h.trace() / 4), sp.Matrix([[1]]))


def at_c(pc, cval):
    return sum((M * sp.Rational(cval) ** k for k, M in pc.items()), sp.zeros(11, 11))


def gilkey(P, Rd, v, Kj, cval, drop_DB=False):
    """omega, E, Tr Omega^2, the ghost's E and the total by (3.40), at one value of c, from FORM's pieces"""
    K0, K1, K2 = [sp.Rational(x) for x in Kj]
    S = PBIG.copy(); S[10, 10] = K0
    Si = S.inv()
    B = [at_c(P['B'][m], cval) for m in idx]
    C = at_c(P['C'], cval)
    w = [Si * B[m] / 2 for m in idx]
    z = sp.zeros(11, 11)
    Dw = [[(z if drop_DB else Si * at_c(P['DB'][n][m], cval) / 2) - Si * (K1 * v[n] * ESS) * Si * B[m] / 2 for m in idx] for n in idx]       # Dw[n][m] = D_n omega^m
    E = (Si * C - sum((Dw[m][m] for m in idx), z) - sum((w[m] * w[m] for m in idx), z)).applyfunc(sp.expand)
    Om2 = sp.Integer(0)
    for m in idx:
        for n in idx:
            Rmn = sp.Matrix(4, 4, lambda a, e: Rd[(a, e, m, n)])
            RH = sp.diag(D.op_matrix10(lambda h: Rmn * h + h * Rmn.T), sp.Matrix([[0]]))        # [D_m, D_n] on h_ab; nothing on the scalar
            Om = RH + Dw[m][n] - Dw[n][m] + w[m] * w[n] - w[n] * w[m]
            Om2 += (Om * Om).trace()
    Om2 = sp.expand(Om2)
    Egh = P['G0'] + sp.Rational(cval) * P['G1']
    tot = sp.expand(GM.hv340(11, Om2, E, P['Rs'], P['inv']['Ric2']) - 2 * GM.hv340(4, -P['inv']['Riem2'], Egh, P['Rs'], P['inv']['Ric2']))
    return dict(S=S, B=B, w=w, E=E, Om2=Om2, Egh=Egh, tot=tot)


def r3_prediction(Rd, u, W, Kj):
    """the prediction R3 of class710_operator_derive.py for E at c = 1 (same construction: build_E by homogeneity, plus delta)"""
    K0, K1, K2 = [sp.Rational(x) for x in Kj]
    k1 = K1 / K0; u2 = (u.T * u)[0]
    Wt = W + k1 * u * u.T / 2; beta = Wt.trace(); delta = sp.expand(k1 * beta / 2)
    Ric, Rs, inv = GM.invariants(Rd, u, Wt)
    E_u, Om2_s3 = GM.build_E(Rd, Ric, Rs, u, Wt)
    E_0, _ = GM.build_E(Rd, Ric, Rs, sp.zeros(4, 1), sp.zeros(4, 4))
    Ep = sp.zeros(11, 11)
    Ep[:10, :10] = E_0[:10, :10] + K0 * (E_u[:10, :10] - E_0[:10, :10])
    Ep[:10, 10] = K0 * E_u[:10, 10]; Ep[10, :10] = E_u[10, :10]; Ep[10, 10] = K0 * E_u[10, 10] + delta
    img = dict(inv); img['Rv2'] = K0 * inv['Rv2']; img['Ricvv'] = K0 * inv['Ricvv']; img['v4'] = K0**2 * inv['v4']; img['trW2'] = K0 * inv['trW2']; img['WW'] = K0 * inv['WW']
    onsh = (Ric - K0 * u * u.T / 2) == sp.zeros(4, 4) and beta == 0
    return dict(Ep=Ep.applyfunc(sp.expand), Om2=Om2_s3, delta=delta, img=img, inv=inv, Rs=Rs, u2=u2, K0=K0, Wt=Wt, onsh=onsh, Ric=Ric)


def r8_row(tot, T):
    return (tot, {'v4': T['img']['v4'], 'vWv': T['K0'] * (T['uWu']), 'ricid': T['img']['WW'] - T['img']['trW2'] + T['img']['Ricvv'],
                  'gb': T['inv']['Riem2'] - 4 * T['inv']['Ric2'] + T['inv']['R2']})


def sample(rng, Kj, onshell):
    R = D.riemann_sample(rng); u, W, cub = D.scalar_sample(rng)
    if onshell:
        k1 = Kj[1] / Kj[0]
        W = W - (W.trace() + k1 * (u.T * u)[0] / 2) * sp.eye(4) / 4
        R = D.with_ricci(R, Kj[0] * u * u.T / 2)
    return R, u, W, cub


class ReadOnlyCheckpoint:
    """a run's <log>.ckpt read as data; nothing is ever written or computed from here"""
    def __init__(self, path):
        self.done = {}
        with open(path) as fh:
            for line in fh:
                try:
                    rec = json.loads(line)
                except ValueError:
                    continue
                self.done[rec['k']] = ([(sp.Rational(a), [sp.Rational(x) for x in g]) for a, g in rec['r']], rec['dt'])

    def put(self, *a):
        raise RuntimeError("the reader must not compute base fields")


def check_main(F, args, struct, tag='main', quiet=False, drop_DB=False):
    """F1.1-F1.6 on fresh backgrounds: random K jets and a = 1/12, three off shell and three on shell each"""
    rng = random.Random(args.seed)
    res = dict(f11=True, f12=True, f13=True, f14=True, f15=True, f16=True)
    phi0 = sp.Rational(1, 2)
    c_, red, on, _ = struct
    for name in ('free K jets', 'a = 1/12'):
        for onshell in (False, True):
            for i in range(args.backgrounds):
                if name == 'a = 1/12':
                    Kj = Z.K_jets_of_a(sp.Rational(1, 12), phi0)
                else:
                    Kj = [sp.Rational(rng.randint(2, 7), rng.randint(1, 3)), sp.Rational(rng.choice((-3, -2, -1, 1, 2, 3)), rng.randint(1, 3)), sp.Rational(rng.choice((-3, -2, -1, 1, 2, 3)), rng.randint(1, 3))]
                    if Kj[0] == 1:
                        Kj[0] = sp.Rational(5, 2)
                R, u, W, cub = sample(rng, Kj, onshell)
                P = F.pieces(R, u, W, Kj, rng)
                T = r3_prediction(R, u, W, Kj)
                G = gilkey(P, R, u, Kj, 1, drop_DB)
                K0, K1, K2 = Kj
                # F1.2: at c = 1 the only first-derivative term is K1 v.D p in the scalar equation; B is linear in c
                f12 = all(G['B'][m] == K1 * u[m] * ESS for m in idx) and all(set(P['B'][m]) <= {0, 1} for m in idx)
                nmix = sum(1 for m in idx for x in P['B'][m].get(1, sp.zeros(11, 11)) if x != 0)
                f13 = G['E'] == T['Ep'] and G['Om2'] == T['Om2']
                Kj2 = [K0, K1, K2 + 1]
                f14 = F.has_K2 and gilkey(F.pieces(R, u, W, Kj2), R, u, Kj2, 1, drop_DB)['E'] == G['E']
                f15 = P['G0'] == P['Ric'] and P['G1'] == -K0 * u * u.T and not P['ghigh'] and bool(P['gbox'])
                pred = sp.expand(sum(c_[k] * T['img'][k] for k in GM.KEYS) + T['delta'] * (T['Rs'] / 6 - K0 * T['u2']) / 2 + T['delta']**2 / 4)
                f16 = G['tot'] == pred
                r7 = ''
                if onshell:
                    assert T['onsh']
                    gb = T['inv']['Riem2'] - 4 * T['inv']['Ric2'] + T['inv']['R2']
                    ricid = T['img']['WW'] - T['img']['trW2'] + T['img']['Ricvv']
                    rem = sp.expand(G['tot'] - c_['Riem2'] * gb - c_['WW'] * ricid)
                    ok7 = rem == on * K0**2 * T['u2']**2
                    f16 = f16 and ok7
                    r7 = f"; R7 remainder = ({on}) K^2 (u^2)^2: {ok7}"
                if not quiet:
                    say(f"   {name}, {'on' if onshell else 'off'} shell {i}: F1.1 {P['s_ok']};  F1.2 {f12} ({nmix} mixing entries at order c);  F1.3 E == R3's prediction and Tr Omega^2 == step 3's: {f13};  F1.4 K2 drops out of E: {f14};  F1.5 ghost: {f15};  F1.6 total == R6's formula: {G['tot'] == pred}{r7}")
                res['f11'] &= bool(P['s_ok']); res['f12'] &= f12; res['f13'] &= f13; res['f14'] &= f14; res['f15'] &= f15; res['f16'] &= f16
    return res


def check_r8(F, args, struct, cvals):
    """F1.7: the R8 fit from FORM's operator at several c, on six on-shell backgrounds (a = 1/12)"""
    rng = random.Random(args.seed + 1)
    Kj = Z.K_jets_of_a(sp.Rational(1, 12), sp.Rational(1, 2))
    bgs = []
    for i in range(6):
        R, u, W, cub = sample(rng, Kj, True)
        T = r3_prediction(R, u, W, Kj); T['uWu'] = (u.T * T['Wt'] * u)[0]
        bgs.append((R, u, F.pieces(R, u, W, Kj), T))
    ok = True; table = []
    for cv in cvals:
        rows = [r8_row(gilkey(P, R, u, Kj, cv)['tot'], T) for R, u, P, T in bgs]
        cf, resid = Z.fit_r8(rows)
        good = cf is not None and cf['v4'] == struct[2] and cf['vWv'] == 0 and all(x == 0 for x in resid)
        say(f"   c = {cv}: alpha = {cf['v4']}, b = {cf['vWv']}, C = {cf['ricid']}, g = {cf['gb']}, residuals {resid}: {'PASS' if good else 'FAIL'}")
        ok &= good; table.append((cv, cf))
    return ok, table


def check_logs(F, args, resdir, lognames=('class710_r8_c0.log', 'class710_r8_c2.log', 'class710_operator_derive.log'), collect=None):
    """F1.8: FORM's operator against the jet-probe runs, on those runs' own backgrounds, from their logs and checkpoints.
    `collect`, if a list, receives (kind, R, u, W, cub, Kj, c, FORM's gilkey data) for every background compared."""
    ok = True; nseen = 0
    saved_say, saved_ver = Z.say, Z.CODE_VERSION
    for logname in lognames:
        path = os.path.join(resdir, logname)
        if not (os.path.exists(path) and os.path.exists(path + '.ckpt')):
            say(f"   {logname}: log or checkpoint not found; not compared"); continue
        text = open(path).read()
        head = re.search(r"seed (\d+), cases (\S+), phi0 = (\S+), (\d+) on-shell \+ (\d+) off-shell backgrounds per case, c = (\S+), gauge (\w+), mode (\w+)", text)
        ver = re.search(r"# code version (\w+)", text).group(1)
        seed, cases, phi0, non, noff, cval, gauge, mode = head.groups()
        cval = sp.Rational(cval.rstrip(',')); phi0 = sp.Rational(phi0.rstrip(','))
        last = text.rsplit("# RESUMED", 1)[-1]          # a resumed run reprints every background below its last RESUMED line
        totals = re.findall(r"\(3\.40\) total, pointwise: ([-0-9/]+)", last)
        st = re.search(r"^status=(\S+)", last, re.M)
        status = st.group(1) if st else "none (the run did not finish)"
        ck = ReadOnlyCheckpoint(path + '.ckpt')
        rng = random.Random(int(seed)); nb = 0; allsame = True
        Z.say = lambda *a: None; Z.CODE_VERSION = ver
        try:
            for case in cases.rstrip(',').split(','):
                if case == 'free':
                    Kj = [sp.Rational(rng.randint(1, 5), rng.randint(1, 3)), sp.Rational(rng.choice((-3, -2, -1, 1, 2, 3)), rng.randint(1, 3)), sp.Rational(rng.randint(-3, 3), rng.randint(1, 3))]
                else:
                    a = sp.Rational(rng.choice((-2, -1, 1, 2, 3)), rng.randint(2, 7)) if case == 'rand' else sp.Rational(case)
                    Kj = Z.K_jets_of_a(a, phi0)
                for kind, n in (('on shell', int(non)), ('off shell', int(noff))):
                    for i in range(n):
                        R, u, W, cub = sample(rng, Kj, kind == 'on shell')
                        bk = Z.Checkpoint.bgkey(R, u, W, cub, Kj, cval, gauge == 'image')
                        J = Z.task_list('min', list(range(11)), 'expand', 0)
                        missing = [j for j in J if f"{bk}|{j[0]}|{j[1]}" not in ck.done]
                        if missing:
                            say(f"   {logname}: {case} {kind} {i}: {len(missing)} of {len(J)} base fields are not in the checkpoint under this key; not compared"); allsame = False; nb += 1; continue
                        dd = Z.derive(R, u, W, cub, Kj, cval, gauge == 'image', 'min', list(range(11)), 'expand', 0, 1, '', time.time(), ck)
                        gam1 = Z.Background(R, u, W, cub, Kj, 0, 1).gam1()
                        w0, Ep_, Om2p = Z.operator_data(dd, gam1, 'none')
                        P = F.pieces(R, u, W, Kj, random.Random(nb))
                        G = gilkey(P, R, u, Kj, cval)
                        c_eff = cval * Kj[0] if gauge == 'image' else cval
                        gh = GH.ghost_operator(D.Background(R, u, W, cub), c_eff, 1)          # the probes' ghost, derived here from the gauge condition (seconds)
                        same = dict(S=bool(P['s_ok']) and dd['S0'] == G['S'], B=all(dd['B0'][m] == G['B'][m] for m in idx), E=Ep_.applyfunc(sp.expand) == G['E'], Om2=Om2p == G['Om2'],
                                    ghost=gh[3] is not None and gh[3].applyfunc(sp.expand) == G['Egh'].applyfunc(sp.expand))
                        if totals:
                            same['total'] = nb < len(totals) and sp.Rational(totals[nb]) == G['tot']
                        if collect is not None:
                            collect.append((kind, R, u, W, cub, Kj, cval, G))
                        nb += 1; nseen += 1
                        if not all(same.values()):
                            allsame = False
                            say(f"   {logname}: {case} {kind} {i}: DIFFERS {same}")
        finally:
            Z.say, Z.CODE_VERSION = saved_say, saved_ver
        if totals and len(totals) != nb:
            allsame = False
            say(f"   {logname}: {len(totals)} pointwise totals found in the log for {nb} backgrounds")
        say(f"   {logname} (c = {cval}, its own status={status}, code version {ver}): {nb} backgrounds regenerated from its seed; S (FORM's second-derivative part), the first-derivative coefficient, E (11 x 11), Tr Omega^2, the ghost"
            + (" and the logged pointwise total" if totals else "") + f" by FORM == by the jet probes: {'all identical' if allsame and nb else 'NO'}")
        ok &= allsame and nb > 0
        if nb and not allsame and not nseen:
            say(f"   {logname}: found, but nothing in its checkpoint matches the backgrounds regenerated from its header (a different code version or options?)")
    return ok, nseen


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261114)
    ap.add_argument('--backgrounds', type=int, default=3)
    ap.add_argument('--resdir', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results'))
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; to start over:  del {args.log.replace('/', chr(92))}"); sys.exit(1)
        out = open(args.log, 'w')
    say(f"# form_class710_check: seed {args.seed}, {args.backgrounds} backgrounds per kind, reading {args.resdir}")
    say(f"# versions: python {platform.python_version()} ({platform.system()} {platform.machine()}), sympy {sp.__version__}; {sys.executable}")
    frm = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'form', 'class710_operator.frm')
    if os.path.exists(frm):
        say(f"#   form/class710_operator.frm as it is now: sha1 {hashlib.sha1(open(frm, 'rb').read().replace(bytes([13, 10]), bytes([10]))).hexdigest()[:12]}")
    have = lambda tag: all(os.path.exists(os.path.join(args.resdir, f"form_class710_{tag}_{nm}.out")) for nm in ('Eh', 'Ep', 'Gh', 'DBh', 'DBp'))
    if not have('main'):
        say("no FORM output for TAG main: run form/class710_operator.frm first (its header has the lines)"); say("status=INCOMPLETE"); sys.exit(1)
    struct = Z.step3_structure(args.seed + 7)
    say(f"# step 3's total on the eight structures: {struct[0]}; on shell {struct[2]} (v^2)^2   [{time.time() - t0:.0f} s]")
    F = FormOperator(args.resdir, 'main')
    say(f"## TAG main: {F.nterms} terms (h equation, scalar equation, ghost, D B for h, D B for the scalar); symbols {F.syms}; tensors {F.names}")
    if not F.readable:
        say(f"NOT READABLE AS AN OPERATOR: unknown objects {F.unknown}; linear in one fluctuation: {F.linear}; K3 present: {'K3' in F.syms}; first-derivative coefficients free of R, Ric, W (the last FORM module assumes it): {F.first_clean}"); say("status=FAIL"); sys.exit(1)
    R = check_main(F, args, struct)
    say("F1.1 {}   F1.2 {}   F1.3 {}   F1.4 {}   F1.5 {}   F1.6 {}   [{:.0f} s]".format(*['PASS' if R[k] else 'FAIL' for k in ('f11', 'f12', 'f13', 'f14', 'f15', 'f16')], time.time() - t0))
    ok = all(R.values())
    say("## F1.7  the gauge coefficient, by FORM: six on-shell backgrounds (a = 1/12), four fitted and two for the residual")
    ok7, table = check_r8(F, args, struct, [sp.Integer(-1), sp.Integer(0), sp.Rational(1, 2), sp.Integer(1), sp.Integer(2), sp.Integer(3)])
    say(f"F1.7 {'PASS: alpha is the same at six values of c; alpha(c) is at most a quartic (E and Omega are quadratic in c, the ghost linear), so it is constant in this family -- given FORM operator and that degree count' if ok7 else 'FAIL'}   [{time.time() - t0:.0f} s]")
    ok &= ok7
    say("## F1.8  FORM against the jet probes on the author's runs")
    ok8, nseen = check_logs(F, args, args.resdir)
    say(f"F1.8 {'PASS' if ok8 else 'FAIL' if nseen else 'NOT TESTED'} ({nseen} backgrounds)   [{time.time() - t0:.0f} s]")
    complete = nseen >= 22 and args.backgrounds >= 3
    ok &= ok8 or nseen == 0
    say("## F1.9  controls")
    Rc = check_main(F, args, struct, quiet=True, drop_DB=True)
    broke = not Rc['f13']
    say(f"F1.9 control: D_n B^m left out of the conversion: F1.3 {'breaks, as it must' if broke else 'STILL HOLDS'}")
    ok &= broke
    for tag, what, must in (('literal', 'the gauge without K', ('f12', 'f13')), ('flip', "the commutator's sign reversed", ('f13',))):
        if not have(tag):
            say(f"F1.9 control {tag} ({what}): no FORM output; not tested"); complete = False; continue
        Fc = FormOperator(args.resdir, tag)
        Rc = check_main(Fc, args, struct, tag, quiet=True) if Fc.readable else {k: False for k in R}
        broke = all(not Rc[k] for k in must)
        say(f"F1.9 control {tag} ({what}): breaks {[k for k in sorted(Rc) if not Rc[k]]}: {'as it must' if broke else 'NOT ENOUGH -- ' + str(must) + ' must break'}")
        ok &= broke
    say("status=" + ("FAIL" if not ok else "PASS" if complete else "INCOMPLETE (PASS needs the control outputs, three backgrounds per kind, and the three jet-probe logs with their checkpoints)"))
    say(f"END form_class710_check {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
