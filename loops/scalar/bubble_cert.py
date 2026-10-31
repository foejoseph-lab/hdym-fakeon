"""bubble_cert.py -- the trusted checker for rung 1.  Recomputes the published closed forms that bubble_checks.py
claims to hit in ball arithmetic (Arb via python-flint: every number is an interval that rigorously contains the
true value) and compares with the generator's numbers in results/*.json.  Nothing here imports bubble_lib.py.

What is certified
  lw.json      1801.00915 eq. 5.6-5.8 (D = 4): f0(t) = int_0^1 ln[(1-2x)^2 + t^2 x^2 (1-x)^2], g(t) = int_0^1 ln[1 + t^2 x^2 (1-x)^2],
               f_AV = f0 + (2 pi/t)(sqrt(4-t^2) - 2) theta(t) theta(2-t) - (4 pi/t) theta(t-2),
               the D = 2 forms f0 = -2t int x(1-x)/K, g = -2t int x(1-x)/H, f_AV = f0 + 4 pi/sqrt(4-t^2) theta(t)theta(2-t),
               the D = 3 Euclidean closed forms of eq. 5.12;  the generator's f_AV, g_AV, |f_d| must lie within tol of them.
  fakeon.json  the exact T +- i0 limits for ordinary masses (D = 4: ln T + sum_j [(1-x_j) ln(1-x_j) + x_j ln x_j - 1] -+ i pi (x2-x1);
               D = 3: (2/sqrt T)[asinh sqrt(x1/(x2-x1)) + asinh sqrt((1-x2)/(x2-x1))] + i pi/sqrt T;  D = 2: the PV log plus i pi
               residues), 2001.01942 eq. 6.6 and 4.1, 1801.00915 eq. 5.4 (= the D = 4 real part);  the generator's f_+, f_AV and
               route-(b) differences must agree.
  fw.json      2001.01942 eq. 2.6 (coefficient), 4.3, 3.1 and the sec. 7 imaginary parts: closed forms are algebraic, the check is
               that the generator's deviations are within the tolerances it declares (1/Lambda^2 for 2.6 and 3.1).
Rigorous integration: acb.integral on [0, 1]; the integrands are holomorphic in a neighbourhood of the segment for the t
used (K and H have no zero on [0, 1] for real t != 0).  The D = 4 real-mass check is algebraic.

Run (Spyder):  !"{sys.executable}" -B -u bubble_cert.py > results\\cert.log
Log ends '# END cert elapsed=<s>s status=PASS|FAIL'.
"""
import json
import os
import sys
import time

from flint import acb, arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
ctx.dps = 40
PI = arb.pi()
TOL = arb("1e-11")       # the generator's declared accuracy (1e-12 plus the O(sqrt delta) threshold terms)
t0 = time.time()
fails = []
worst = arb(0)
worst_what = ""


def say(*a):
    print(*a)
    sys.stdout.flush()


def C(v):
    """json value -> acb: [re, im] decimal strings (or floats), or a single string/float."""
    if isinstance(v, (list, tuple)):
        return acb(arb(str(v[0])), arb(str(v[1])))
    return acb(arb(str(v)))


def check(what, gen, cert, tol=TOL):
    """gen: json value (see C) or an arb/acb; cert: arb/acb ball.  Passes if the rigorous upper bound of
    |gen - cert| is <= tol.  The generator's value is a point (decimal string at its working precision); cert
    carries its own radius, reported."""
    global worst, worst_what
    g = gen if isinstance(gen, (arb, acb)) else C(gen)
    d = abs(acb(g) - acb(cert))        # arb ball
    dmax = d.upper()                   # rigorous upper bound
    ok = bool(dmax <= arb(tol))
    rad = acb(cert).rad()
    if bool(dmax > worst):
        worst, worst_what = dmax, what
    say(f"    [{'ok' if ok else 'FAIL'}] {what}: |gen - cert| <= {dmax.str(3)}  (cert radius {rad.str(2)})")
    if not ok:
        fails.append(what)


def load(name):
    p = os.path.join(RES, f"{name}.json")
    if not os.path.exists(p):
        say(f"  (no {name}.json; run bubble_checks.py --{name} first)")
        return None
    with open(p) as fh:
        return json.load(fh)


def integral01(f):
    return acb.integral(f, 0, 1)


# ------------------------------------------------------------------------------------------------------------------
say(f"# cert: start, ctx.dps = {ctx.dps}, tol = {TOL.str(2)}")

# ---- lw.json ------------------------------------------------------------------------------------------------------
lw = load("lw")
if lw:
    say("## 1801.00915 eqs. 5.6-5.8 (D = 4) and the D = 2 closed forms, rigorous integrals")
    for row in lw["D4"]:
        t = arb(row["T"])
        f0 = integral01(lambda x, analytic, t=t: ((1 - 2 * x) ** 2 + t * t * x * x * (1 - x) ** 2).log(analytic=analytic))
        g = integral01(lambda x, analytic, t=t: (1 + t * t * x * x * (1 - x) ** 2).log(analytic=analytic))
        f0, g = f0.real, g.real
        if t < 2:
            fav = f0 + (2 * PI / t) * ((4 - t * t).sqrt() - 2)
            fd = arb(0)
        else:
            fav = f0 - 4 * PI / t
            fd = (2 * PI / t) * (t * t - 4).sqrt()
        say(f"  D=4 T={row['T']}: f0 = {f0.str(20)}  g = {g.str(20)}  f_AV = {fav.str(20)}")
        check(f"D=4 T={row['T']} f_AV", row["f_AV"], acb(fav))
        check(f"D=4 T={row['T']} g_AV", row["g_AV"], acb(g))
        check(f"D=4 T={row['T']} |f_d| vs continued published term", abs(C(row["f_d"])), acb(fd))
        check(f"D=4 T={row['T']} f0 (published, vs generator f0_pub)", row["f0_pub"], f0)
    for row in lw["D2"]:
        t = arb(row["T"])
        f0 = -2 * t * integral01(lambda x, analytic, t=t: x * (1 - x) / ((1 - 2 * x) ** 2 + t * t * x * x * (1 - x) ** 2))
        g = -2 * t * integral01(lambda x, analytic, t=t: x * (1 - x) / (1 + t * t * x * x * (1 - x) ** 2))
        f0, g = f0.real, g.real
        if t < 2:
            fav = f0 + 4 * PI / (4 - t * t).sqrt()
            fd = arb(0)
        else:
            fav = f0
            fd = 4 * PI / (t * t - 4).sqrt()
        say(f"  D=2 T={row['T']}: f0 = {f0.str(20)}  g = {g.str(20)}  f_AV = {fav.str(20)}")
        check(f"D=2 T={row['T']} f_AV", row["f_AV"], acb(fav))
        check(f"D=2 T={row['T']} g_AV", row["g_AV"], acb(g))
        check(f"D=2 T={row['T']} |f_d| vs continued published term", abs(C(row["f_d"])), acb(fd))
    say("## 1801.00915 eq. 5.12, Euclidean closed forms (D = 3), against the generator's route (b)")
    e = lw["D3_euclid"]
    t = arb(e["t"])
    r = (-t).sqrt()
    I = acb(0, 1)
    f0 = (2 * I / r) * ((acb(2).sqrt() - I * r) / (acb(2).sqrt() + I * r)).log()
    si, smi = I.sqrt(), (-I).sqrt()
    g = (I / r) * (((2 * si + r) / (2 * si - r)).log() + ((2 * smi - r) / (2 * smi + r)).log())
    tot = f0 + g
    gen = C(e["J_mom"]) * 32 * PI / acb(0, 1)
    say(f"  t={e['t']}: f0+g = {tot.str(20)}   generator route (b) = {gen.str(20)}")
    check("D=3 Euclidean eq. 5.12 vs route (b)", gen, tot)

# ---- fakeon.json ----------------------------------------------------------------------------------------------------
fk = load("fakeon")
if fk:
    say("## ordinary masses: exact T +- i0 limits and 2001.01942 eqs. 6.6, 4.1")
    for key, (a1, a2) in {"D4_m1=m2=1": (arb(1), arb(1)), "D4_m1=1,m2=1.5": (arb(1), arb("2.25"))}.items():
        m1, m2 = a1.sqrt(), a2.sqrt()
        thr = (m1 + m2) ** 2
        for row in fk[key]:
            T = arb(row["T"])
            if T > thr:
                # roots of T x^2 + (a2 - a1 - T) x + a1
                B = a2 - a1 - T
                disc = (B * B - 4 * T * a1).sqrt()
                x1, x2 = (-B - disc) / (2 * T), (-B + disc) / (2 * T)
                Lc = lambda c: (1 - c) * (1 - c).log() + c * c.log() - 1
                fplus = acb(T.log() + Lc(x1) + Lc(x2), -PI * (x2 - x1))
                fav = fplus.real
                up, um = (m1 + m2) ** 2 - T, (m1 - m2) ** 2 - T
                d66 = (up * um).sqrt() / (16 * PI * T)
                im41 = (up * um).sqrt() / (8 * PI * T)
            else:
                # below threshold: rigorous integral of ln Delta (Delta > 0 on [0, 1])
                fplus = integral01(lambda x, analytic, a1=a1, a2=a2, T=T: (a1 * (1 - x) + a2 * x - T * x * (1 - x)).log(analytic=analytic))
                fav = fplus.real
                d66 = arb(0)
                im41 = arb(0)
            say(f"  {key} T={row['T']}: f_+ = {fplus.str(18)}")
            check(f"{key} T={row['T']} f_+", row["f_plus"], fplus)
            check(f"{key} T={row['T']} f_AV", row["f_AV"], acb(fav))
            dB = C(row["B_PV_mom"]) - C(row["B_F_mom"])
            check(f"{key} T={row['T']} eq. 6.6 (B_PV - B_F, route b)", dB.real, d66)
            check(f"{key} T={row['T']} eq. 6.6 imaginary part zero (identity by construction)", dB.imag, arb(0))
            bf = C(row["B_F_mom"])
            check(f"{key} T={row['T']} eq. 4.1 (2 Im(-i B_F))", 2 * (acb(0, -1) * bf).imag, im41)
            if "eq54" in row:
                check(f"{key} T={row['T']} eq. 5.4 generator value vs exact f_AV", row["eq54"], fav)
    say("## D = 3 and D = 2 ordinary masses: the exact limits against the generator's route (b)")
    for D in (3, 2):
        pref = {3: 1 / (8 * PI), 2: 1 / (4 * PI)}[D]
        for row in fk[f"D{D}"]:
            a1, a2 = (arb(1), arb(1)) if row["lab"] == "m1=m2=1" else (arb(1), arb("2.25"))
            T = arb(row["T"])
            thr = (a1.sqrt() + a2.sqrt()) ** 2
            if T > thr:
                B = a2 - a1 - T
                disc = (B * B - 4 * T * a1).sqrt()
                x1, x2 = (-B - disc) / (2 * T), (-B + disc) / (2 * T)
                if D == 3:
                    rT = T.sqrt()
                    re = (2 / rT) * ((x1 / (x2 - x1)).sqrt().asinh() + ((1 - x2) / (x2 - x1)).sqrt().asinh())
                    fplus = acb(re, PI / rT)
                else:
                    pv = (((1 - x1) / (1 - x2)).log() - (x1 / x2).log()) / (T * (x1 - x2))
                    fplus = acb(pv, 2 * PI / (T * (x2 - x1)))
            else:
                # below threshold Delta > 0 on [0, 1]: rigorous integral, no imaginary part
                if D == 3:
                    fplus = integral01(lambda x, analytic, a1=a1, a2=a2, T=T: (a1 * (1 - x) + a2 * x - T * x * (1 - x)).rsqrt(analytic=analytic))
                else:
                    fplus = integral01(lambda x, analytic, a1=a1, a2=a2, T=T: 1 / (a1 * (1 - x) + a2 * x - T * x * (1 - x)))
            BF = acb(0, 1) * pref * fplus          # B_F = i * pref * f_+
            BPV = acb(0, 1) * pref * fplus.real
            check(f"D={D} {row['lab']} T={row['T']} B_F route (b) vs exact", row["B_F_mom"], BF)
            check(f"D={D} {row['lab']} T={row['T']} B_PV route (b) vs exact", row["B_PV_mom"], BPV)

# ---- fw.json --------------------------------------------------------------------------------------------------------
fw = load("fw")
if fw:
    say("## Feynman-Wheeler: 2001.01942 eqs. 2.6, 4.3, 3.1 and sec. 7 (algebraic closed forms)")
    a1, a2 = arb(1), arb("0.25")
    m1, m2 = a1.sqrt(), a2.sqrt()
    for row in fw["m1_m2"]:
        T = arb(row["T"])
        c26 = acb(0, (a1 - a2) / ((4 * PI) ** 2 * T))
        # declared: converges as 1/Lambda^2; at Lambda = 500..5000 the generator's tolerance is 2e-5 relative
        check(f"T={row['T']} eq. 2.6 coefficient (tol 2e-5 rel)", row["slope"], c26, tol=arb("2e-5") * abs(c26.imag))
        up, um = (m1 + m2) ** 2 - T, (m1 - m2) ** 2 - T
        c43 = -(up * um).sqrt() / (8 * PI * T) if um > 0 else arb(0)
        check(f"T={row['T']} eq. 4.3", row["im"], c43)
    for row in fw["m2_0"]:
        T = arb(row["T"])
        Lam = arb(row["Lam"])
        rhs = -(1 / T) * (4 * Lam * Lam).log() - ((T - 1) / T) * acb(T - 1).log() + ((T - 1) / T) * T.log()
        # acb(T-1).log() has imaginary part pi for T < 1 (principal branch of a negative real; the +i eps side)
        check(f"T={row['T']} Lambda={row['Lam']} eq. 3.1 (tol 10/Lambda^2)", row["lhs"], rhs, tol=10 / (Lam * Lam))
    for row in fw["D3_D2"]:
        T = arb(row["T"])
        up, um = (m1 + m2) ** 2 - T, (m1 - m2) ** 2 - T
        if um > 0:
            pub = -(1 / (4 * T.sqrt())) if row["D"] == 3 else -(1 / (up * um).sqrt())
        else:
            pub = arb(0)
        check(f"D={row['D']} T={row['T']} sec. 7 imaginary part", row["im"], pub)

el = time.time() - t0
status = "PASS" if not fails else "FAIL(" + ",".join(fails) + ")"
say(f"# worst |gen - cert| {worst.str(3)} at {worst_what}")
say(f"# END cert elapsed={el:.1f}s status={status}")
