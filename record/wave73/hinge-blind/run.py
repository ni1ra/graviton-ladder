"""HINGE-BLIND production runner.

  python run.py one  <disc> <n> <mult> <l> [prec] [sets]     one case, prints and writes JSON
  python run.py plan <planname> <workers>                     run a named plan with a process pool

A case = one lattice vacuum (disc, n, N = mult*n, l) and a list of perturbation sets applied to
its modular Hamiltonian. Output: results/raw/<disc>_n<n>_N<N>_l<l>_p<prec>.json (skipped if present).

Perturbation sets (Z = (g/2) a sum_j w_j phi_j^2):
  mid    a = 1/(n+1/2), w = 1                 g = +1, -1, +1/4, -1/4      (primary, PREREG)
  site   a = 1/n,       w_n = 1/2             g = +1, -1                  (boundary-site variant)
  w0     a = 1/(n+1/2), w_n = 0               g = +1, -1                  (planted error)
  w2     a = 1/(n+1/2), w_n = 2               g = +1, -1                  (planted error)
  half   a = 1/(n+1/2), w = 1, Z without 1/2  g = +1, -1  (i.e. g -> 2g)  (planted error)
  small  a = 1/(n+1/2), w = 1                 g = +-1/16, +-1/2           (extra g-scaling points)
  pm1    a = 1/(n+1/2), w = 1                 g = +1, -1                  (primary, large runs)
Every case also stores Q2 = l(l+1) lim_{g->0}(1-A)/g^2 from second-order perturbation theory.
"""
import json, math, os, sys, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "results", "raw")


def est_digits(n, l):
    """Empirical size of -log10 min(nu^2 - 1/4), from pilots (see logs/pilot.log)."""
    return n * (3.7 + 1.15 * math.log2(max(l, 2) / 2.0)) + 20


def default_prec(n, l):
    return int(math.ceil((est_digits(n, l) + 220) / 0.30103 / 64.0) * 64)


def pert_sets(n, names):
    mid_a = Fraction(2, 2 * n + 1)
    out = []
    for nm in names:
        if nm == "mid":
            out += [("mid", g, mid_a, [1] * n, g) for g in (Fraction(1), Fraction(-1), Fraction(1, 4), Fraction(-1, 4))]
        elif nm == "pm1":
            out += [("mid", g, mid_a, [1] * n, g) for g in (Fraction(1), Fraction(-1))]
        elif nm == "site":
            out += [("site", g, Fraction(1, n), [1] * (n - 1) + [Fraction(1, 2)], g) for g in (Fraction(1), Fraction(-1))]
        elif nm == "w0":
            out += [("w0", g, mid_a, [1] * (n - 1) + [0], g) for g in (Fraction(1), Fraction(-1))]
        elif nm == "w2":
            out += [("w2", g, mid_a, [1] * (n - 1) + [2], g) for g in (Fraction(1), Fraction(-1))]
        elif nm == "half":
            out += [("half", g, mid_a, [1] * n, 2 * g) for g in (Fraction(1), Fraction(-1))]
        elif nm == "small":
            out += [("mid", g, mid_a, [1] * n, g) for g in (Fraction(1, 16), Fraction(-1, 16), Fraction(1, 2), Fraction(-1, 2))]
        else:
            raise ValueError(nm)
    return out   # (label, nominal g, a, w, effective g)


def run_case(disc, n, mult, l, prec=None, sets=("mid", "site"), verbose=False, force=False):
    import hb
    from flint import arb
    N = mult * n
    prec = prec or default_prec(n, l)
    os.makedirs(RAW, exist_ok=True)
    path = os.path.join(RAW, "%s_n%d_N%d_l%d_p%d_%s.json" % (disc, n, N, l, prec, "-".join(sets)))
    if os.path.exists(path) and not force:
        return json.load(open(path))
    log = (lambda s: print(s, flush=True)) if verbose else None
    t0 = time.time()
    hb.setprec(prec)
    XA, PA, dg = hb.vacuum_blocks(n, N, l, disc, prec, log)
    t1 = time.time()
    mo = hb.modular(XA, PA)
    t2 = time.time()
    D = -mo["diag"]["log10_min_delta"]
    margin = prec * 0.30103 - D
    res = dict(disc=disc, n=n, N=N, l=l, prec=prec, D=D, margin_digits=margin,
               vac_diag=dg, mod_diag=mo["diag"], entropy=mo["entropy"].str(25, radius=False),
               eps_min=mo["eps"][-1].str(20, radius=False), eps_max=mo["eps"][0].str(20, radius=False),
               t_vac=t1 - t0, t_mod=t2 - t1, cases=[])
    if margin < 100:
        res["error"] = "precision margin %.0f digits < 100" % margin
        json.dump(res, open(path + ".lowprec", "w"), indent=1)
        return res
    nu, eps = mo["nu"], mo["eps"]
    A1 = hb.mid(hb.scale_cols(mo["W"], [1 / x.sqrt() for x in nu]))
    A1i = hb.mid(hb.diagmat([x.sqrt() for x in nu]) * mo["Winv"])
    ll = l * (l + 1)
    res["Q2_mid"] = (hb.second_order(eps, A1, Fraction(2, 2 * n + 1), [1] * n) * ll).str(30, radius=False)
    res["Q2_site"] = (hb.second_order(eps, A1, Fraction(1, n), [1] * (n - 1) + [Fraction(1, 2)]) * ll).str(30, radius=False)
    for label, g, a, w, geff in pert_sets(n, sets):
        tc = time.time()
        c = dict(set=label, g=str(g), a=str(a))
        try:
            Mp = hb.perturbed_M(mo["M"], a, w, geff)
            E2, A2, A2i, d2 = hb.williamson(Mp, mo["N"])
            lA = hb.affinity_I(eps, A1, A1i, E2, A2, A2i)
            oma = hb.one_minus_A(lA)
            Q = (oma * ll / (arb(g.numerator) / g.denominator) ** 2).mid()
            c.update(one_minus_A=oma.str(30, radius=False), Q=Q.str(30, radius=False),
                     epsp_min=E2[0].str(20, radius=False), epsp_max=E2[-1].str(20, radius=False), wdiag=d2)
        except ArithmeticError as ex:
            c.update(error=str(ex))
        c["t"] = time.time() - tc
        res["cases"].append(c)
        if verbose:
            print("   %s g=%s  1-A=%s  Q=%s  (%.0fs)" % (label, g, c.get("one_minus_A"), c.get("Q"), c["t"]), flush=True)
    res["t_total"] = time.time() - t0
    json.dump(res, open(path, "w"), indent=1)
    return res


def run_pt2(disc, n, mult, l, prec=None):
    """Second-order coefficient Q2 = l(l+1) lim_{g->0} (1-A)/g^2 for the mid and site weight sets."""
    import hb
    N = mult * n
    prec = prec or default_prec(n, l)
    d = os.path.join(HERE, "results", "pt2")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "%s_n%d_N%d_l%d_p%d.json" % (disc, n, N, l, prec))
    if os.path.exists(path):
        return json.load(open(path))
    t0 = time.time()
    hb.setprec(prec)
    XA, PA, dg = hb.vacuum_blocks(n, N, l, disc, prec)
    mo = hb.modular(XA, PA)
    nu, eps = mo["nu"], mo["eps"]
    A1 = hb.mid(hb.scale_cols(mo["W"], [1 / x.sqrt() for x in nu]))
    ll = l * (l + 1)
    q_mid = hb.second_order(eps, A1, Fraction(2, 2 * n + 1), [1] * n) * ll
    q_site = hb.second_order(eps, A1, Fraction(1, n), [1] * (n - 1) + [Fraction(1, 2)]) * ll
    res = dict(disc=disc, n=n, N=N, l=l, prec=prec, D=-mo["diag"]["log10_min_delta"],
               Q2_mid=q_mid.str(30, radius=False), Q2_site=q_site.str(30, radius=False), t_total=time.time() - t0)
    json.dump(res, open(path, "w"), indent=1)
    return res


def _job_pt2(args):
    try:
        r = run_pt2(*args)
        return "pt2 %s n=%d N=%d l=%d prec=%d t=%.0fs  Q2_mid=%s Q2_site=%s" % (
            r["disc"], r["n"], r["N"], r["l"], r["prec"], r["t_total"], r["Q2_mid"][:16], r["Q2_site"][:16])
    except Exception as ex:
        return "FAILED pt2 %s: %r" % (args, ex)


def _job(args):
    try:
        r = run_case(*args)
        qs = ["%s%+.3g:%s" % (c["set"], float(Fraction(c["g"])), c.get("Q", c.get("error", "?"))[:14]) for c in r["cases"]]
        return "%s n=%d N=%d l=%d prec=%d D=%.0f margin=%.0f t=%.0fs  %s %s" % (
            r["disc"], r["n"], r["N"], r["l"], r["prec"], r["D"], r["margin_digits"], r.get("t_total", 0),
            " ".join(qs), r.get("error", ""))
    except Exception as ex:   # keep the pool alive, report the failure
        return "FAILED %s: %r" % (args, ex)


LS = (2, 3, 4, 6, 8, 12, 16, 24, 32)
PLANS = {
    # primary grid, small to large
    "grid_small": [("D2", n, 4, l, None, ("mid", "site")) for n in (24, 32, 48, 64) for l in LS],
    "grid_96": [("D2", 96, 4, l, None, ("mid", "site")) for l in LS],
    "grid_128": [("D2", 128, 4, l, None, ("mid", "site")) for l in LS],
    "grid_192": [("D2", 192, 4, l, None, ("mid", "site")) for l in (2, 4, 8, 16, 32, 3, 6, 12, 24)],
    # secondary discretisation
    "sred": [("Sred", n, 4, l, None, ("mid", "site")) for n in (24, 48, 96) for l in (2, 4, 8, 16, 32)],
    "sred_big": [("Sred", n, 4, l, None, ("mid",)) for n in (128, 192) for l in (2, 4, 8, 16, 32)],
    # outside size
    "outside": [("D2", n, m, l, None, ("mid",)) for n in (24, 48, 96) for m in (2, 8, 16) for l in (2, 4, 8, 32)],
    # planted errors and extra g points
    "planted": [("D2", n, 4, l, None, ("w0", "w2", "half", "small")) for n in (24, 48, 96) for l in (2, 8, 32)],
    # precision reruns (about 1.4x the bits)
    "precision": [("D2", n, 4, l, int(default_prec(n, l) * 1.4 // 64 * 64), ("mid",)) for n in (24, 48, 96) for l in (2, 8, 32)],
    # full l ladder at fixed n for direct partial sums
    "pt2_small": [("pt2", "D2", n, 4, l) for n in (24, 32, 48, 64, 96) for l in LS] + [("pt2", "Sred", n, 4, l) for n in (24, 48, 96) for l in (2, 4, 8, 16, 32)],
    "pt2_128": [("pt2", "D2", 128, 4, l) for l in LS],
    "pt2_192": [("pt2", "D2", 192, 4, l) for l in LS],
    "pt2_256": [("pt2", "D2", 256, 4, l) for l in (2, 4, 8, 16, 32)],
    "grid_256": [("D2", 256, 4, l, None, ("pm1",)) for l in (8, 32, 16, 2, 4, 24, 12, 48, 64)],
    "extra_n": [("D2", n, 4, l, None, ("pm1",)) for n in (40, 56, 80, 112) for l in LS + (48, 64)],
    "extra_160": [("D2", 160, 4, l, None, ("pm1",)) for l in LS + (48, 64)],
    "hi_l": [("D2", n, 4, l, None, ("pm1",)) for n in (64, 96, 128, 192) for l in (48, 64)],
    "ladder48": [("D2", 48, 4, l, None, ("mid",)) for l in list(range(0, 49)) if l not in LS],
}

if __name__ == "__main__":
    if sys.argv[1] == "one":
        disc, n, mult, l = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        prec = int(sys.argv[6]) if len(sys.argv) > 6 and sys.argv[6] != "-" else None
        sets = tuple(sys.argv[7].split(",")) if len(sys.argv) > 7 else ("mid", "site")
        r = run_case(disc, n, mult, l, prec, sets, verbose=True, force=True)
        print(json.dumps({k: v for k, v in r.items() if k != "cases"}, indent=1))
    else:
        import multiprocessing as mp
        jobs = PLANS[sys.argv[2]]
        workers = int(sys.argv[3])
        print("plan %s: %d jobs, %d workers, pid %d" % (sys.argv[2], len(jobs), workers, os.getpid()), flush=True)
        with mp.Pool(workers, maxtasksperchild=1) as pool:
            if jobs and jobs[0][0] == "pt2":
                it = pool.imap_unordered(_job_pt2, [j[1:] for j in jobs])
            else:
                it = pool.imap_unordered(_job, jobs)
            for line in it:
                print(time.strftime("%H:%M:%S"), line, flush=True)
        print("plan %s done" % sys.argv[2], flush=True)
