"""s4: a finite-dimensional sanity check of the KMS criterion (modflow Lemma 4) as THEOREM.md uses it.

A standard subspace K of C^(2n) is built from a modular pair (J, Delta): K = {x : J Delta^(1/2) x = x}.
The criterion says: if V(s) K = K and s -> <v, V(s) u> continues to the strip -pi <= Im z <= 0 with real
values on Im z = -pi, then V(s) = Delta^(i s / 2 pi).  Here we check the converse normalisation, which is
where a sign or a factor 2 pi could hide: for V(s) = Delta^(i s / 2 pi) the boundary value at s - i pi is
real for all u, v in K, and it is not real for the planted rates 1.1/(2 pi), -1/(2 pi) and 1/pi.
Inner product antilinear in the first slot, as in wave70/modflow.
"""
import numpy as np
from lib import Report

R = Report('s4_kms_toy')
rng = np.random.default_rng(20261110)
n = 4
lam = rng.uniform(0.2, 5.0, n)
# block model: H = C^n (+) C^n, J(a, b) = (conj b, conj a), Delta = diag(lam, 1/lam)
U = np.linalg.qr(rng.normal(size=(2 * n, 2 * n)) + 1j * rng.normal(size=(2 * n, 2 * n)))[0]
D = np.concatenate([lam, 1 / lam])


def in_K(a):
    """The vector (a, sqrt(lam) conj a), rotated by U."""
    return U @ np.concatenate([a, np.sqrt(lam) * np.conj(a)])


def J_Dhalf(xvec):
    w = np.conj(U).T @ xvec
    w = np.sqrt(D) * w
    a, b = w[:n], w[n:]
    return U @ np.concatenate([np.conj(b), np.conj(a)])


def V(zc, rate):
    return U @ np.diag(D**(1j * zc * rate)) @ np.conj(U).T


us = [in_K(rng.normal(size=n) + 1j * rng.normal(size=n)) for _ in range(6)]
R.check('S = J Delta^(1/2) fixes the sample vectors (they lie in K)', max(np.abs(J_Dhalf(u) - u).max() for u in us) < 1e-12)
R.check('K is real-linear but not complex: i u is not fixed', min(np.abs(J_Dhalf(1j * u) - 1j * u).max() for u in us) > 1e-3)


def worst_imag(rate):
    w = 0.0
    for s in (-1.3, 0.0, 0.7, 2.1):
        for u in us:
            for v in us:
                val = np.vdot(v, V(s - 1j * np.pi, rate) @ u)
                w = max(w, abs(val.imag) / (abs(val) + 1e-300))
    return w


good = worst_imag(1 / (2 * np.pi))
R.check('V(s) = Delta^(i s/2 pi): <v, V(s - i pi) u> is real on K', good < 1e-12, good)
R.value('worst relative imaginary part, correct rate', good)
for name, rate in (('rate 1.1/(2 pi)', 1.1 / (2 * np.pi)), ('opposite sign', -1 / (2 * np.pi)), ('rate 1/pi (pi for 2 pi)', 1 / np.pi)):
    bad = worst_imag(rate)
    R.plant('P ' + name + ': boundary value not real', bad > 1e-3, bad)
# invariance V(s) K = K
R.check('V(s) K = K', max(np.abs(J_Dhalf(V(0.9, 1 / (2 * np.pi)) @ u) - V(0.9, 1 / (2 * np.pi)) @ u).max() for u in us) < 1e-12)
R.write()
