"""Planted-error control: the same pipeline in float64, and the high-precision pipeline with the
entanglement energies capped, against the reference high-precision value.

  python control_double.py <n> <l>

float64 'naive': modes with nu^2 - 1/4 <= 0 make the modular Hamiltonian undefined (reported).
float64 'clip' : nu^2 - 1/4 clipped from below at 1e-14 (the usual double-precision patch).
arb 'cap E'    : exact covariances, entanglement energies eps_k capped at E.
"""
import sys
from fractions import Fraction
import numpy as np
import hb
from flint import arb
from run import default_prec

n, l = int(sys.argv[1]), int(sys.argv[2])
N = 4 * n
a = 1.0 / (n + 0.5)
ll = l * (l + 1)


def aff_np(eps1, A1, A1i, eps2, A2, A2i):
    def half(eps, A, Ai):
        nut = 0.5 / np.tanh(eps / 4)
        return (A * nut) @ A.T, (Ai.T * nut) @ Ai, np.sum(0.5 * np.log(2 * np.sinh(eps / 2)) - np.log(2 * np.sinh(eps / 4)))
    X1, P1, n1 = half(eps1, A1, A1i)
    X2, P2, n2 = half(eps2, A2, A2i)
    return n1 + n2 - 0.5 * (np.linalg.slogdet(X1 + X2)[1] + np.linalg.slogdet(P1 + P2)[1])


def double_pipeline(clip):
    d, e = hb.lattice_coeffs(N, l, "D2")
    K = np.diag([float(x) for x in d]) + np.diag([float(x) for x in e], 1) + np.diag([float(x) for x in e], -1)
    w, V = np.linalg.eigh(K)
    X = ((V * w ** -0.5) @ V.T / 2)[:n, :n]
    P = ((V * w ** 0.5) @ V.T / 2)[:n, :n]
    L = np.linalg.cholesky(X)
    S = L.T @ P @ L
    ev, U = np.linalg.eigh((S + S.T) / 2)
    delta = ev - 0.25
    nbad = int((delta <= 0).sum())
    if clip is None and nbad:
        return None, nbad
    if clip is not None:
        delta = np.maximum(delta, clip)
    nu = np.sqrt(delta + 0.25)
    eps = np.log((nu + 0.5) / (delta / (nu + 0.5)))
    W = L @ U
    Wi = U.T @ np.linalg.inv(L)
    Nm = (W * (eps / nu)) @ W.T
    Mm = (Wi.T * (eps * nu)) @ Wi
    A1 = W / np.sqrt(nu)
    A1i = (Wi.T * np.sqrt(nu)).T
    out = {}
    for g in (1.0, -1.0):
        Mp = Mm + 2 * np.pi * g * a * np.eye(n)
        LN = np.linalg.cholesky((Nm + Nm.T) / 2)
        T = LN.T @ Mp @ LN
        e2, U2 = np.linalg.eigh((T + T.T) / 2)
        if e2[0] <= 0:
            out[g] = float("nan")
            continue
        ep = np.sqrt(e2)
        A2 = (LN @ U2) / np.sqrt(ep)
        A2i = (U2 * np.sqrt(ep)).T @ np.linalg.inv(LN)
        out[g] = ll * (-np.expm1(aff_np(eps, A1, A1i, ep, A2, A2i)))
    return out, nbad


prec = default_prec(n, l)
hb.setprec(prec)
XA, PA, dg = hb.vacuum_blocks(n, N, l, "D2", prec)
mo = hb.modular(XA, PA)


def arb_Q(cap):
    eps = mo["eps"] if cap is None else [x if x < cap else arb(cap) for x in mo["eps"]]
    nu = mo["nu"]
    k = len(eps)
    W, Wi = mo["W"], mo["Winv"]
    # with capped eps the state is exp(-K_cap); its own nu follows from eps: nu = coth(eps/2)/2, same modes
    nuc = [((x / 2).cosh() / (x / 2).sinh() / 2).mid() for x in eps]
    # normal-mode map of the *original* covariances: phi = W nu^{-1/2} Q
    A1 = hb.mid(hb.scale_cols(W, [1 / x.sqrt() for x in nu]))
    A1i = hb.mid(hb.diagmat([x.sqrt() for x in nu]) * Wi)
    Mm = hb.mid(hb.scale_cols(A1i.transpose(), eps) * A1i)
    Nm = hb.mid(hb.scale_cols(A1, eps) * A1.transpose())
    out = {}
    for g in (1, -1):
        Mp = hb.perturbed_M(Mm, Fraction(2, 2 * n + 1), [1] * n, g)
        E2, A2, A2i, _ = hb.williamson(Mp, Nm)
        out[g] = float(hb.one_minus_A(hb.affinity_I(eps, A1, A1i, E2, A2, A2i)) * ll)
    return out


ref = arb_Q(None)
print("n=%d l=%d N=%d  reference (arb, %d bits):  Q(+1)=%.10f  Q(-1)=%.10f  eps_max=%.1f" % (
    n, l, N, prec, ref[1], ref[-1], float(mo["eps"][0])))
for cap in (16, 32, 64, 128, 256, 512):
    q = arb_Q(cap)
    print("  arb, eps capped at %4d:  Q(+1)=%.10f (%+.2e rel)  Q(-1)=%.10f (%+.2e rel)   modes capped: %d of %d" % (
        cap, q[1], q[1] / ref[1] - 1, q[-1], q[-1] / ref[-1] - 1, sum(1 for x in mo["eps"] if x > cap), n))
for clip in (None, 1e-14, 1e-12, 1e-10):
    out, nbad = double_pipeline(clip)
    if out is None:
        print("  float64 naive: %d of %d modes have nu^2 - 1/4 <= 0: modular Hamiltonian undefined, no number" % (nbad, n))
    else:
        print("  float64 clip=%s (%d modes had nu^2-1/4<=0):  Q(+1)=%.10f (%+.2e rel)  Q(-1)=%.10f (%+.2e rel)" % (
            clip, nbad, out[1.0], out[1.0] / ref[1] - 1, out[-1.0], out[-1.0] / ref[-1] - 1))
