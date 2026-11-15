"""Development probe: timing and accuracy of the building blocks (not a measurement)."""
import sys, time
import numpy as np
from flint import arb, arb_mat, acb_mat, ctx
import hb

n = int(sys.argv[1]); mult = int(sys.argv[2]); l = int(sys.argv[3]); prec = int(sys.argv[4])
alg = sys.argv[5] if len(sys.argv) > 5 else "approx"
N = mult * n
t = time.time()
XA, PA, dg = hb.vacuum_blocks(n, N, l, "D2", prec, log=print)
print("vacuum_blocks %.1fs" % (time.time() - t), dg)
# float64 comparison
d, e = hb.lattice_coeffs(N, l, "D2")
K = np.diag([float(x) for x in d]) + np.diag([float(x) for x in e], 1) + np.diag([float(x) for x in e], -1)
w, V = np.linalg.eigh(K)
Xd = (V * w ** -0.5) @ V.T / 2
Xf = np.array([[float(XA[i, j]) for j in range(n)] for i in range(n)])
print("max |X_arb - X_double| =", abs(Xf - Xd[:n, :n]).max())
hb.setprec(prec)
t = time.time()
L = hb.cholesky(XA)
print("cholesky %.1fs" % (time.time() - t))
Lm = hb.to_mat(L)
S = hb.mid(Lm.transpose() * PA * Lm)
t = time.time()
if alg == "approx":
    ev, U, dg2 = hb.eigh_arb(S)
    print("eig %.1fs" % (time.time() - t), dg2)
    nu2 = ev
else:
    E = acb_mat(S).eig(nonstop=True)
    print("eig(default) %.1fs" % (time.time() - t))
    nu2 = sorted([z.real for z in E], key=lambda z: float(z.mid()))
q = hb.arb(1) / 4
print("smallest nu^2-1/4:", [(x - q).str(8) for x in nu2[:4]])
print("largest  nu^2-1/4:", [(x - q).str(8) for x in nu2[-4:]])
