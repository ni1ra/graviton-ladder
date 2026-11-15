"""E5: on the ball lattice, even perturbation: affinity formula I against formula II (independent
determinant route), and N' positivity margin (smallest eigenvalue of N'^-1-normalised form)."""
from fractions import Fraction
import hb, even

for n, l, prec in ((12, 2, 1600), (16, 4, 2400), (24, 8, 4000)):
    hb.setprec(prec)
    XA, PA, dg = hb.vacuum_blocks(n, 4 * n, l, "D2", prec)
    mo = hb.modular(XA, PA)
    E1, A1, A1i, _ = hb.williamson(mo["M"], mo["N"])
    for conv in ("ghost", "site"):
        a, G = even.kb_matrix(n, l, conv)
        for g in (1, -1):
            Np = even.perturbed_N(mo["N"], a, G, Fraction(g))
            E2, A2, A2i, _ = hb.williamson(mo["M"], Np)
            lI = hb.affinity_I(E1, A1, A1i, E2, A2, A2i)
            lII, rel = hb.affinity_II(E1, A1, A1i, E2, A2, A2i)
            print("n=%d l=%d %s g=%+d  1-A(I)=%s  log10 rel|I-II|=%.1f (det radius %.1f)  eps' range [%s, %s]" % (
                n, l, conv, g, hb.one_minus_A(lI).str(20, radius=False), hb.l10((lI - lII) / lI), rel,
                E2[0].str(8, radius=False), E2[-1].str(8, radius=False)))
