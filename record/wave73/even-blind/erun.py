"""EVEN-BLIND runner.

  python erun.py one  <n> <mult> <l> <tag> [prec]
  python erun.py plan <planname> <workers>

A case = one lattice vacuum (D2, n, N = mult*n, l) + a list of specs (set by tag):
  ("q2", conv)            second-order coefficient for momentum perturbation with G of `conv`
  ("ex", conv, g, geff)   exact affinity, N' = N + 2 pi geff a G(conv); Q normalised by g^2
  ("odd", g)              hinge-blind's primary Z_odd (mid set), exact affinity, plus Q2_mid
Output: results/raw/E_n<n>_N<N>_l<l>_p<prec>_<tag>.json, never overwritten.
"""
import json, math, os, sys, time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "results", "raw")


def est_digits(n, l):
    return n * (3.7 + 1.15 * math.log2(max(l, 2) / 2.0)) + 20


def default_prec(n, l):
    return int(math.ceil((est_digits(n, l) + 220) / 0.30103 / 64.0) * 64)


PM1 = (F(1), F(-1))
TAGS = {
    "main": [("q2", "ghost"), ("q2", "trunc"), ("q2", "site")] + [("ex", "ghost", g, g) for g in PM1],
    "odd": [("odd", g) for g in PM1],
    "ctrl": [("q2", c) for c in ("ghost", "local", "neumann", "nocent")]
            + [("ex", "ghost", g, g) for g in (F(1, 4), F(-1, 4), F(1, 16), F(-1, 16))]
            + [("ex", "trunc", g, g) for g in PM1] + [("ex", "site", g, g) for g in PM1]
            + [("ex", "local", g, g) for g in (F(1), F(-1), F(1, 4), F(-1, 4))]
            + [("ex", "neumann", g, g) for g in PM1] + [("ex", "nocent", g, g) for g in PM1]
            + [("ex", "ghost", g, 2 * g) for g in PM1],
    "out": [("q2", "ghost")] + [("ex", "ghost", g, g) for g in PM1],
    "prec": [("q2", "ghost")] + [("ex", "ghost", g, g) for g in PM1],
    "local": [("q2", "local")] + [("ex", "local", g, g) for g in PM1],
    "q2only": [("q2", "ghost"), ("q2", "trunc"), ("q2", "site")],
}


def run_case(n, mult, l, tag, prec=None, verbose=False):
    import hb, even
    from flint import arb
    N = mult * n
    prec = prec or default_prec(n, l)
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, "E_n%d_N%d_l%d_p%d_%s.json" % (n, N, l, prec, tag))
    if os.path.exists(path):
        return json.load(open(path))
    t0 = time.time()
    hb.setprec(prec)
    XA, PA, dg = hb.vacuum_blocks(n, N, l, "D2", prec)
    mo = hb.modular(XA, PA)
    D = -mo["diag"]["log10_min_delta"]
    margin = prec * 0.30103 - D
    res = dict(n=n, N=N, l=l, prec=prec, tag=tag, D=D, margin_digits=margin, vac_diag=dg, mod_diag=mo["diag"],
               entropy=mo["entropy"].str(25, radius=False), t_base=time.time() - t0, cases=[], q2={})
    if margin < 100:
        res["error"] = "precision margin %.0f digits < 100" % margin
        json.dump(res, open(path + ".lowprec", "w"), indent=1)
        return res
    nu, eps = mo["nu"], mo["eps"]
    A1 = hb.mid(hb.scale_cols(mo["W"], [1 / x.sqrt() for x in nu]))
    A1i = hb.mid(hb.diagmat([x.sqrt() for x in nu]) * mo["Winv"])
    ll = l * (l + 1)
    Gs = {}

    def getG(conv):
        if conv not in Gs:
            Gs[conv] = even.kb_matrix(n, l, conv)
        return Gs[conv]

    for spec in TAGS[tag]:
        tc = time.time()
        if spec[0] == "q2":
            a, G = getG(spec[1])
            res["q2"][spec[1]] = (even.q2_momentum(eps, A1i, a, G) * ll).str(30, radius=False)
            continue
        if spec[0] == "odd":
            g = spec[1]
            a = F(2, 2 * n + 1)
            if "odd_mid" not in res["q2"]:
                res["q2"]["odd_mid"] = (hb.second_order(eps, A1, a, [1] * n) * ll).str(30, radius=False)
            c = dict(set="odd_mid", g=str(g))
            Mp, Np = hb.perturbed_M(mo["M"], a, [1] * n, g), mo["N"]
            geff = g
        else:
            _, conv, g, geff = spec
            a, G = getG(conv)
            c = dict(set=conv if g == geff else conv + "_x2", g=str(g))
            Mp, Np = mo["M"], even.perturbed_N(mo["N"], a, G, geff)
        try:
            E2, A2, A2i, d2 = hb.williamson(Mp, Np)
            lA = hb.affinity_I(eps, A1, A1i, E2, A2, A2i)
            oma = hb.one_minus_A(lA)
            Q = (oma * ll / (arb(g.numerator) / g.denominator) ** 2).mid()
            c.update(one_minus_A=oma.str(30, radius=False), Q=Q.str(30, radius=False),
                     epsp_min=E2[0].str(20, radius=False), epsp_max=E2[-1].str(20, radius=False), wdiag=d2)
        except (ArithmeticError, ValueError) as ex:
            c.update(error=str(ex)[:200])
        c["t"] = time.time() - tc
        res["cases"].append(c)
        if verbose:
            print("   %s g=%s 1-A=%s Q=%s (%.0fs) %s" % (c["set"], g, c.get("one_minus_A"), c.get("Q"), c["t"], c.get("error", "")), flush=True)
    res["t_total"] = time.time() - t0
    json.dump(res, open(path, "w"), indent=1)
    return res


def _job(args):
    try:
        r = run_case(*args)
        qs = ["%s%+.3g:%s" % (c["set"], float(F(c["g"])), c.get("Q", "ERR")[:13]) for c in r["cases"]]
        q2 = ["q2_%s:%s" % (k, v[:13]) for k, v in r["q2"].items()]
        return "n=%d N=%d l=%d %s prec=%d D=%.0f margin=%.0f t=%.0fs  %s %s %s" % (
            r["n"], r["N"], r["l"], r["tag"], r["prec"], r["D"], r["margin_digits"], r.get("t_total", 0),
            " ".join(q2), " ".join(qs), r.get("error", ""))
    except Exception as ex:
        return "FAILED %s: %r" % (args, ex)


NS = (24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 144, 160, 176, 192)
LS = (2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48)
C3 = [(n, l) for n in (24, 48, 96) for l in (2, 8, 32)]


def grid(ns, ls=LS, tag="main"):
    js = [(n, 4, l, tag) for n in ns for l in ls if n >= 2 * l]
    return sorted(js, key=lambda j: -(j[0] ** 3) * default_prec(j[0], j[2]) ** 1.5)


PLANS = {
    "odd": [(24, 4, 8, "odd"), (32, 4, 4, "odd"), (48, 4, 16, "odd")],
    "pilot": grid((24, 32)),
    "mid": grid((40, 48, 56, 64, 80, 96)),
    "big": grid((112, 128, 144, 160, 176, 192)),
    "ctrl": [(n, 4, l, "ctrl") for n, l in C3],
    "out": [(n, m, l, "out") for n in (24, 48, 96) for m in (2, 8) for l in (2, 4, 8, 32)],
    "prec": [(n, 4, l, "prec", int(default_prec(n, l) * 1.4 // 64 * 64)) for n, l in ((24, 8), (48, 16), (96, 32))],
    "local192": [(192, 4, l, "local") for l in (32, 8, 2)],
    "big224": grid((224,), (48, 40, 32, 24, 16, 8, 4)),
    "hil": grid((128, 144, 160, 176, 192), (64,)),
}

if __name__ == "__main__":
    if sys.argv[1] == "one":
        n, mult, l, tag = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
        prec = int(sys.argv[6]) if len(sys.argv) > 6 else None
        r = run_case(n, mult, l, tag, prec, verbose=True)
        print(json.dumps({k: v for k, v in r.items() if k != "cases"}, indent=1))
    else:
        import multiprocessing as mp
        jobs = PLANS[sys.argv[2]]
        workers = int(sys.argv[3])
        print("plan %s: %d jobs, %d workers, pid %d" % (sys.argv[2], len(jobs), workers, os.getpid()), flush=True)
        with mp.Pool(workers, maxtasksperchild=1) as pool:
            for line in pool.imap_unordered(_job, jobs):
                print(time.strftime("%H:%M:%S"), line, flush=True)
        print("plan %s done" % sys.argv[2], flush=True)
