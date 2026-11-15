"""Exploration only (closed-form side of the bound, limit model l = inf). No grid operator, no beta1hat.
Chooses (k, b, y0, eps) for u(y) = (1 - ((y-y0)/b)^2)^k by scanning the analytic lower bound
c = (8/15) A exp(-2 eps G^2/A) s* / (Q sqrt(2 S' C))   (see logs/01)."""
import math, numpy as np
TAU = math.pi
def consts(k, b, y0, n=40001):
    y = np.linspace(y0 - b, y0 + b, n); s = (y - y0) / b; dy = y[1] - y[0]
    w = 1 - s * s
    u = w ** k
    u1 = k * w ** (k - 1) * (-2 * s) / b
    u2 = (k * (k - 1) * w ** (k - 2) * 4 * s * s - 2 * k * w ** (k - 1)) / b ** 2
    V = np.exp(-2 * y)
    g = -u2 + V * u
    I = lambda f: float(np.sum(f) * dy)
    return I(u1 * u1 + V * u * u), math.sqrt(I(g * g)), I(np.exp(y) * (u1 * u1 + V * u * u))
P = np.linspace(1e-6, 40, 400001)
def SC(eps):
    lg = -2 * eps * P * P + np.log(np.sinh(np.minimum(2 * math.pi * P, 700)))
    return float(np.max(P * np.exp(lg))), float(np.max(P / np.tanh(math.pi * P) * np.exp(-2 * eps * P * P)))
def bound(A, G, Q, eps):
    J0 = A * math.exp(-eps * G * G / A); J1 = min(G * G, A / (math.e * eps))
    ss = math.sqrt(2 * J0 / J1)
    if ss <= TAU: M = 8 * J0 * J0 * ss / 15 / Q
    else: M = (J0 * J0 * TAU - J0 * J1 * TAU ** 3 / 3 + J1 * J1 * TAU ** 5 / 20) / Q
    S, C = SC(eps)
    return M / (2 * math.sqrt(S * A) * math.sqrt(C * A / 2)), ss
best = (0,)
for k in (2, 3, 4):
  for b in (1, 1.5, 2, 3, 4, 6, 8):
    for y0 in (-0.5, 0, 0.5, 1, 1.5, 2, 3, 4):
        A, G, Q = consts(k, b, y0)
        for eps in (0.25, 0.5, 1, 1.5, 2, 3, 4, 6):
            c, ss = bound(A, G, Q, eps)
            if c > best[0]: best = (c, k, b, y0, eps, A, G, Q, ss)
print("best", best)
for k, b, y0 in ((3, 2, 0), (3, 3, 1), (3, 4, 1), (2, 3, 1)):
    A, G, Q = consts(k, b, y0)
    print(k, b, y0, "A %.4g G %.4g Q %.4g" % (A, G, Q), [(e, "%.3g" % bound(A, G, Q, e)[0]) for e in (0.5, 1, 2, 3, 4)])
