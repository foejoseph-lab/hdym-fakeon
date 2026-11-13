#!/usr/bin/env python3
"""
Rung 3, step 3b, second method for route 2 (FORM), stage F0: the reader of FORM's output.

form/gravscalar_operator.frm derives, symbolically and covariantly, the Euler-Lagrange expressions of step 3's gauge-fixed
second-order Lagrangian (K = 1, the gauge coefficient c a symbol) and the ghost operator, and writes them as sums of tensor
monomials to results/form_gravscalar_<TAG>_{Eh,Ep,Gh}.out.  Its header carries the predictions F0.1-F0.6 (written 13 Nov
17:11 ET, before any FORM code; the author's yes 17:25 ET) and how the derivation is done.  This script reads those files,
evaluates the monomials on random backgrounds in exact rationals (plain index sums; no numpy), and tests the predictions:

  F0.1  second-derivative terms: with only the symmetrised second derivatives S2, pS2 switched on, the equations are
        P Box h and Box phi (P = 1/2 - 1/4 g tr, the weights of gravscalar_operator_derive.py), for every c.
  F0.2  first-derivative terms: with only H1 = D h and p1 = D phi on, zero at c = 1 and proportional to (1 - c).
  F0.3  the algebraic part at c = 1, as an 11 x 11 matrix, equals hv1974_gravscalar_M.build_E exactly.
  F0.4  the ghost: Box eta + (Ric - c v v) eta, no first-derivative terms.
  F0.5  controls, from FORM runs with the commutator's sign reversed (TAG flip) and with the 1/8 h^2 of sqrt(g) dropped
        (TAG nosq): each must break F0.3.
  For the record, not predicted: how the algebraic part depends on c (it is not Gilkey's E at c != 1: the first-derivative
  terms then contribute through omega).

WHAT EACH CHECK WOULD AND WOULD NOT CATCH.  F0.3 is the content: every block of the operator, the curvature convention, the
expansions of sqrt(g), the inverse metric and the Ricci tensor, the integration by parts.  F0.1 the kinetic normalisation and the
cancellation of the non-minimal terms.  F0.2 the gauge condition's scalar term.  NOT caught: an error common to FORM's
convention (iii) and build_E's target (the flip control shows the sign is seen, not that it is right -- that rests on the jet-probe
derivation having matched the same target); anything about traces, (3.40), the ghost's -2, the reductions.  S2, H1 and h are
independent data at a point only because second derivatives were symmetrised with the commutator made explicit; a wrong
commutator therefore shows up in F0.3, not in F0.1.
SANDBOX RECORD (13 Nov, FORM 4.3.1 built from source, python 3.13.18, sympy 1.13.1; this block written 13 Nov 17:31 ET, clock read by
the command that wrote it).  FORM: 0.00 s per run, 35 + 8 + 7 terms.  This script, 2 s, three backgrounds per output set:
F0.1, F0.2, F0.3, F0.4 hold on TAG main; TAG flip breaks F0.3 (the h-h block) and the ghost; TAG nosq breaks F0.3 (the h-h
block) and nothing else.  FORM also reads the program correctly with Windows line endings (same output).  For the record: at
c != 1 the algebraic part changes in (E phi)_ab (linear in c) and in the scalar-scalar entry (c^2); the (E h)_phi row does not.
Nothing here has run on the author's machine.

FRESH REVIEW (13 Nov evening, separate agent, no context, scratch copies).  The pass is real: an independent brute-force
evaluator agrees with this one on all 50 terms; the scalar equation was re-derived by hand and matches FORM's 8 terms; the
FORM program is typed correctly by reading.  Folded in: an unknown object in FORM's output used to count as zero (a term
"7*Foo(ia,ib)*h" appended to the file passed everything) -- now only the twelve known names are accepted and every term must be
linear in exactly one fluctuation; the random data can no longer be zero (on one seed with one background a vanishing p1.v hid
a dropped term); PASS needs three backgrounds.
Its 21 mutations of the FORM program: 17 caught (every expansion coefficient, the gauge condition, the by-parts sign, the
commutator, the ghost).  Four NOT caught, because on this Lagrangian they change nothing -- so F0 does not test them, and
"F0.3 tests the integration by parts and the Ricci tensor" above is too strong by exactly these:
  * the order of the two derivatives moved off D D dh (the coefficient is symmetric in them here);
  * the order of the derivative indices in the first term of DT (the commutator terms it would add cancel);
  * the second-order term of the difference of connections GG (it only enters at third order);
  * the Leibniz lines for H2 -> H3, W -> D W, Ric -> D Ric never fire.
  They are right by reading and would be exercised for the first time in F1, if at all.
Also too strong: "a formula in c" is compared with an independent target only at c = 1; at other c only the (1 - c) of the
first-derivative terms and the ghost are tested, and a common factor on the mixing terms would pass F0.2.  And the non-minimal
second-derivative terms cancel only after this reader imposes the symmetries of S2 -- FORM is told none, so its own output
still shows them (and five curvature terms that vanish by symmetry).
This script needs sympy and python-flint (through gravscalar_operator_derive), no numpy.

Run (from graviton_loop, after the three FORM runs in WSL listed in form/gravscalar_operator.frm):
      python3 form_operator_check.py --log results/form_operator_check.log
"""
import sys, os, re, time, random, argparse, platform, hashlib
from fractions import Fraction
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gravscalar_operator_derive as D
import hv1974_gravscalar_M as GM

idx = range(4)
out = None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def fr(x):
    x = sp.Rational(x)
    return Fraction(int(x.p), int(x.q))


def parse(text):
    """FORM expression -> list of (rational coefficient, power of c, power of p, [(tensor name, (index names))])"""
    s = re.sub(r'\s+', '', text).replace('_?', '')
    terms = []
    for t in re.findall(r'[+-]?[^+-]+', s):
        coef = Fraction(-1 if t[0] == '-' else 1); t = t.lstrip('+-')
        cpow = ppow = 0; tens = []
        for f in re.findall(r'[A-Za-z_0-9]+\([^)]*\)|[^*]+', t):
            if '(' in f:
                name, args = f[:-1].split('(')
                tens.append((name, tuple(args.split(','))))
            elif re.fullmatch(r'\d+(/\d+)?', f):
                coef *= Fraction(f)
            elif re.fullmatch(r'c(\^\d+)?', f):
                cpow += int(f[2:]) if '^' in f else 1
            elif re.fullmatch(r'p(\^\d+)?', f):
                ppow += int(f[2:]) if '^' in f else 1
            else:
                raise ValueError(f"unparsed factor {f!r} in term {t!r}")
        terms.append((coef, cpow, ppow, tens))
    return terms


def sparse(nested, rank):
    """nested lists (or a sympy Matrix) -> {index tuple: Fraction}, zeros left out"""
    d = {}

    def walk(x, key):
        if len(key) == rank:
            v = fr(x)
            if v:
                d[key] = v
            return
        for i in idx:
            walk(x[i], key + (i,))
    walk(nested, ())
    return d


def contract(factors, free):
    """product of sparse tensors with labelled indices, summed over every label not in `free`; returns {free index tuple: value}"""
    labels = (); cur = {(): Fraction(1)}
    for i, (labs, data) in enumerate(factors):
        later = set(free)
        for l2, _ in factors[i + 1:]:
            later |= set(l2)
        new_labels = labels + tuple(l for j, l in enumerate(labs) if l not in labels and l not in labs[:j])
        pos = {l: j for j, l in enumerate(new_labels)}
        nxt = {}
        for key, val in cur.items():
            for k2, v2 in data.items():
                full = list(key) + [None] * (len(new_labels) - len(key)); ok = True
                for l, ix in zip(labs, k2):
                    j = pos[l]
                    if full[j] is None:
                        full[j] = ix
                    elif full[j] != ix:
                        ok = False; break
                if ok:
                    full = tuple(full)
                    nxt[full] = nxt.get(full, 0) + val * v2
        keep = tuple(l for l in new_labels if l in later)
        kp = [pos[l] for l in keep]
        cur = {}
        for key, val in nxt.items():
            k = tuple(key[j] for j in kp)
            cur[k] = cur.get(k, 0) + val
        labels = keep
    order = [labels.index(f) for f in free]
    res = {}
    for key, val in cur.items():
        k = tuple(key[j] for j in order)
        res[k] = res.get(k, 0) + val
    return res


def evaluate(terms, tensors, pval, free):
    """-> {power of c: {free index tuple: value}}; a term containing a tensor that is switched off (None) is zero"""
    res = {}
    for coef, cpow, ppow, tens in terms:
        if ppow and not pval:
            continue
        if any(tensors[name] is None for name, _ in tens):        # a name that is not in the table at all is an error, not zero
            continue
        val = contract([(labs, tensors[name]) for name, labs in tens], free)
        acc = res.setdefault(cpow, {})
        for k, x in val.items():
            acc[k] = acc.get(k, 0) + coef * Fraction(pval) ** ppow * x
    return res


def check_tag(tag, resdir, seed, nbg, quiet=False):
    """-> dict of results for one FORM output set"""
    read = lambda nm: parse(open(os.path.join(resdir, f"form_gravscalar_{tag}_{nm}.out")).read())
    Eh, Ep, Gh = read('Eh'), read('Ep'), read('Gh')
    names = sorted({n for t in Eh + Ep + Gh for n, _ in t[3]})
    leftover = sorted({'H2', 'H3', 'DW', 'DRic', 'p2', 'dh', 'dH1', 'dH2', 'dp1', 'E2'} & set(names))
    known = {'h', 'H1', 'S2', 'p1', 'pS2', 'eta', 'SE2', 'R', 'Ric', 'v', 'W', 'd_'}
    unknown = sorted(set(names) - known)
    fluct = {'h', 'H1', 'S2', 'p1', 'pS2'}
    lin = all(sum(1 for n, _ in t[3] if n in fluct) + t[2] == 1 for t in Eh + Ep) and all(sum(1 for n, _ in t[3] if n in ('eta', 'SE2')) == 1 and t[2] == 0 for t in Gh)
    if unknown or not lin:
        say(f"## TAG {tag}: NOT READABLE AS AN OPERATOR: unknown objects {unknown}; every term linear in one fluctuation: {lin}")
        return dict(f01=False, f02=False, f03=False, f04=False, left=False, cdep=set())
    if not quiet:
        say(f"## TAG {tag}: {len(Eh)} + {len(Ep)} + {len(Gh)} terms (h equation, scalar equation, ghost); tensors: {names}; leftover objects that must not be there: {leftover}")
    rng = random.Random(seed)
    rq = lambda: Fraction(rng.choice([k for k in range(-9, 10) if k]), rng.randint(1, 5))          # never zero (review: accidental zeros hid a term on one seed)
    Pbig = sp.diag(D.op_matrix10(lambda h: h / 2 - sp.eye(4) * h.trace() / 4), sp.Matrix([[1]]))
    Pinv = Pbig.inv()
    R = dict(f01=True, f02=True, f03=True, f04=True, left=not leftover, cdep=set())
    for trial in range(nbg):
        Rd = D.riemann_sample(rng); v, W, cub = D.scalar_sample(rng)
        Ric, Rs, inv = GM.invariants(Rd, v, W)
        bg = {'R': {k: fr(x) for k, x in Rd.items() if x != 0}, 'Ric': sparse(Ric.tolist(), 2), 'v': sparse(list(v), 1), 'W': sparse(W.tolist(), 2),
              'd_': {(i, i): Fraction(1) for i in idx}}
        off = {'h': None, 'H1': None, 'S2': None, 'p1': None, 'pS2': None, 'eta': None, 'SE2': None}

        def O(jets, pval):
            """the 11 Euler-Lagrange components, per power of c, as sympy column vectors"""
            T = dict(bg); T.update(off); T.update(jets)
            eh = evaluate(Eh, T, pval, ('ia', 'ib')); ep = evaluate(Ep, T, pval, ())
            res = {}
            for cp in set(eh) | set(ep):
                a = eh.get(cp, {})
                vec = [(a.get((i, j), 0) + a.get((j, i), 0)) / 2 for i, j in D.PAIRS] + [ep.get(cp, {}).get((), 0)]
                res[cp] = sp.Matrix([sp.Rational(x.numerator, x.denominator) if x else 0 for x in vec])
            return res
        tot = lambda res: sum(res.values(), sp.zeros(11, 1))
        # F0.1
        S2 = {}; pS2 = {}
        for a in idx:
            for b in idx:
                for m in idx:
                    for n in idx:
                        if a <= b and m <= n:
                            x = rq()
                            for k in ((a, b, m, n), (b, a, m, n), (a, b, n, m), (b, a, n, m)):
                                S2[k] = x
                if a <= b:
                    x = rq(); pS2[(a, b)] = x; pS2[(b, a)] = x
        boxh = sp.Matrix(4, 4, lambda a, b: sum(sp.Rational(S2[(a, b, m, m)].numerator, S2[(a, b, m, m)].denominator) for m in idx))
        Ph = boxh / 2 - sp.eye(4) * boxh.trace() / 4
        want = sp.Matrix([Ph[a, b] for a, b in D.PAIRS] + [sum(sp.Rational(pS2[(m, m)].numerator, pS2[(m, m)].denominator) for m in idx)])
        r = O({'S2': S2, 'pS2': pS2}, 0)
        f01 = r.get(0, sp.zeros(11, 1)) == want and all(r[k] == sp.zeros(11, 1) for k in r if k != 0)
        # F0.2
        H1 = {}
        for a in idx:
            for b in idx:
                for m in idx:
                    if a <= b:
                        x = rq(); H1[(a, b, m)] = x; H1[(b, a, m)] = x
        p1 = {(m,): rq() for m in idx}
        r = O({'H1': H1, 'p1': p1}, 0)
        z = sp.zeros(11, 1)
        f02 = tot(r) == z and r.get(0, z) == -r.get(1, z) and all(r[k] == z for k in r if k > 1)
        nfirst = sum(1 for x in r.get(1, z) if x != 0)
        # F0.3
        cols = {}
        for k in range(11):
            if D.COMPS[k] == 's':
                r = O({}, 1)
            else:
                a, b = D.COMPS[k]
                r = O({'h': {(a, b): Fraction(1), (b, a): Fraction(1)}}, 0)
            for cp, vec in r.items():
                cols.setdefault(cp, {})[k] = vec
        Ec = {cp: Pinv * sp.Matrix.hstack(*[cols[cp].get(k, z) for k in range(11)]) for cp in cols}
        E1 = sum(Ec.values(), sp.zeros(11, 11)).applyfunc(sp.expand)
        Eprint = GM.build_E(Rd, Ric, Rs, v, W)[0].applyfunc(sp.expand)
        blocks = dict(hh=E1[:10, :10] == Eprint[:10, :10], col=E1[:10, 10] == Eprint[:10, 10], row=E1[10, :10] == Eprint[10, :10], ss=E1[10, 10] == Eprint[10, 10])
        f03 = all(blocks.values())
        for cp, M in Ec.items():
            if cp and M != sp.zeros(11, 11):
                for nm, sub in (('h-h', M[:10, :10]), ('(E phi)_ab', M[:10, 10]), ('(E h)_phi', M[10, :10]), ('phi-phi', M[10:, 10:])):
                    if sub != sp.zeros(*sub.shape):
                        R['cdep'].add(f"{nm} ~ c^{cp}")
        # F0.4
        SE2 = {}
        for a in idx:
            for m in idx:
                for n in idx:
                    if m <= n:
                        x = rq(); SE2[(a, m, n)] = x; SE2[(a, n, m)] = x
        T = dict(bg); T.update(off); T['SE2'] = SE2
        g2 = evaluate(Gh, T, 0, ('ia',))
        box_ok = all(g2.get(0, {}).get((m,), 0) == sum(SE2[(m, n, n)] for n in idx) for m in idx) and all(not any(g2[k].values()) for k in g2 if k)
        M0 = sp.zeros(4, 4); M1 = sp.zeros(4, 4); higher = False
        for k in idx:
            T = dict(bg); T.update(off); T['eta'] = {(k,): Fraction(1)}
            g = evaluate(Gh, T, 0, ('ia',))
            for cp, dct in g.items():
                for (m,), x in dct.items():
                    if cp == 0:
                        M0[m, k] = sp.Rational(x.numerator, x.denominator)
                    elif cp == 1:
                        M1[m, k] = sp.Rational(x.numerator, x.denominator)
                    elif x:
                        higher = True
        f04 = box_ok and M0 == Ric and M1 == -(v * v.T) and not higher and 'E1' not in names
        if not quiet:
            say(f"   background {trial}: F0.1 {f01};  F0.2 zero at c = 1 and proportional to (1 - c): {f02} ({nfirst} nonzero entries);  F0.3 algebraic part at c = 1 == build_E: {f03} {blocks};  F0.4 ghost = Box eta + (Ric - c v v) eta: {f04}")
        R['f01'] &= f01; R['f02'] &= f02; R['f03'] &= f03; R['f04'] &= f04
    return R


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261113)
    ap.add_argument('--backgrounds', type=int, default=3)
    ap.add_argument('--resdir', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results'))
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; to start over:  del {args.log.replace('/', chr(92))}"); sys.exit(1)
        out = open(args.log, 'w')
    say(f"# form_operator_check: seed {args.seed}, {args.backgrounds} random off-shell backgrounds per FORM output set, reading {args.resdir}")
    say(f"# versions: python {platform.python_version()} ({platform.system()} {platform.machine()}), sympy {sp.__version__}; {sys.executable}")
    have = lambda tag: all(os.path.exists(os.path.join(args.resdir, f"form_gravscalar_{tag}_{nm}.out")) for nm in ('Eh', 'Ep', 'Gh'))
    frm = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'form', 'gravscalar_operator.frm')
    if os.path.exists(frm):
        say(f"#   form/gravscalar_operator.frm as it is now: sha1 {hashlib.sha1(open(frm, 'rb').read().replace(bytes([13, 10]), bytes([10]))).hexdigest()[:12]} (the outputs carry no record of the program that wrote them: run FORM, then this, without editing in between)")
    for tag in ('main', 'flip', 'nosq'):
        if have(tag):
            fn = os.path.join(args.resdir, f"form_gravscalar_{tag}_Eh.out")
            say(f"#   {os.path.basename(fn)}: {os.path.getsize(fn)} bytes, written {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(fn)))}")
    if not have('main'):
        say("no FORM output for TAG main: run form/gravscalar_operator.frm first (its header has the lines)"); say("status=INCOMPLETE"); sys.exit(1)
    Rm = check_tag('main', args.resdir, args.seed, args.backgrounds)
    ok = all(Rm[k] for k in ('f01', 'f02', 'f03', 'f04', 'left'))
    say(f"   for the record (not predicted): blocks of the algebraic part that depend on c: {sorted(Rm['cdep']) or 'none'}")
    say(f"F0.1 {'PASS' if Rm['f01'] else 'FAIL'}   F0.2 {'PASS' if Rm['f02'] else 'FAIL'}   F0.3 {'PASS' if Rm['f03'] else 'FAIL'}   F0.4 {'PASS' if Rm['f04'] else 'FAIL'}")
    complete = True
    for tag, what in (('flip', "the commutator's sign reversed"), ('nosq', 'the 1/8 h^2 of sqrt(g) dropped')):
        if not have(tag):
            say(f"F0.5 control {tag} ({what}): no FORM output; not tested"); complete = False; continue
        Rc = check_tag(tag, args.resdir, args.seed, args.backgrounds)
        broke = not Rc['f03']
        say(f"F0.5 control {tag} ({what}): F0.3 {'breaks, as it must' if broke else 'STILL HOLDS -- the check does not see it'}" + (f"; also breaks: {[k for k in ('f01', 'f02', 'f04') if not Rc[k]]}" if broke else ""))
        ok &= broke
    if args.backgrounds < 3:
        complete = False
    say("status=" + ("FAIL" if not ok else "PASS" if complete else "INCOMPLETE (PASS needs both control outputs and at least three backgrounds)"))
    say(f"END form_operator_check {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
