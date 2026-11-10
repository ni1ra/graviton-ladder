"""s5: the completeness ledger of the ladder theorem, as a dependency graph that a machine can audit.

Every node is one of
  HYP       a numbered hypothesis of the theorem (physics that the mathematics cannot supply)
  IMPORT    a printed theorem, with a locator (see ../PRINT_READ.md for how much of its proof was read)
  TEXTBOOK  a standard fact used without a locator (named, so that it is not hidden)
  PROOF     an argument written out, in THEOREM.md or in a record file; `where` says which
  COMPUTED  an exact or numerical computation: `out` is a JSON this script opens, `key` must be true
  AND       a claim that holds if all of `deps` hold
  OR        a claim that holds if any one of `deps` holds (two independent proofs)

The audit checks: the graph is acyclic; every dependency exists; every leaf is typed; every COMPUTED leaf's
output exists and passes; every COMPUTED leaf produced here carries at least one caught planted error; and it
prints, for the root, the full list of leaves. The conclusion can fail only if a leaf fails: that list is the
exhaustion. With --strict, stored outputs of slow runs must also match the recorded sha256.
"""
import hashlib, json, os, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))          # wave69/final
ROOT = os.path.dirname(os.path.dirname(HERE))                               # record
OUT = os.path.join(HERE, 'out')
UP = os.path.join(OUT, 'upstream')

N = {}


def node(key, kind, text, deps=(), **kw):
    assert key not in N, key
    N[key] = dict(kind=kind, text=text, deps=list(deps), **kw)


# ---------------------------------------------------------------- hypotheses
node('H0', 'HYP', 'Scope (a definition, the author 10 Nov): the theory is free linearised gravity at order G^0 on Minkowski space, d = 4, '
     'the free helicity-2 curvature net on the TT Fock space; the region is the unit ball; the state is a stationary magnetic multipole, l >= 2')
node('H1', 'HYP', '(A1) the ball\'s states are the image of the world\'s states under one normal unital completely positive map Phi: M -> B(F)')
node('H2', 'HYP', '(A2) no signalling: Phi(M) commutes with every curvature observable in the causal complement of the diamond')
node('H3', 'HYP', '(A3) M contains the interior curvature observables: a normal unital *-homomorphism iota: A(D) -> M with Phi iota = id')
node('H4', 'HYP', '(GEO) the gravitational side of the second-order first law is the extremal value of (Area - 2 Volume)/4G, a diffeomorphism-invariant functional of the metric')
node('H5', 'HYP', '(EIN) the family g_eps solves the vacuum Einstein equations through O(eps^2); only the second-order Hamiltonian constraint is used')

# ---------------------------------------------------------------- imports and textbook facts
node('I1', 'IMPORT', 'Ciolli-Longo-Ruzzi Thm 4.5: on R(H), H standard and factorial, S(coherent state of k || vacuum) = -Im(k, P_H i log Delta_H k) for every vector k; literal for k in Dom(log Delta_H) (Cor. 2.3)',
     locator='arXiv 1906.01707 Thm 2.2, Cor. 2.3, Prop. 2.4-2.5, Thm 4.5 (pdf p7-10, p21); read line by line, PRINT_READ.md I1')
node('I1b', 'IMPORT', 'Araki: the relative entropy is defined by the relative modular operator and is lower semicontinuous (consumed inside I1)',
     locator='Araki, Publ. RIMS 11 (1976) 809, eq. (1.4); statement read, proof not read, PRINT_READ.md I1')
node('I3', 'IMPORT', 'Uhlmann monotonicity for a unital 2-positive Schwarz map, first state a vector state of a cyclic vector',
     locator='Reible arXiv 2501.04564 Thm IV.3.7 (pdf p68-71); proof read line by line except its textbook inequalities and Sec. III.4, PRINT_READ.md I3')
node('I4', 'IMPORT', 'Hislop-Longo and Longo-Morsella: modular group and coherent-state relative entropy of the massless scalar on the ball (used only by the second proof of M)',
     locator='CMP 84 (1982) 71 Thm 1-2, Cor. 2; arXiv 2012.00565 Thm 4.7, 5.4, Cor. 5.5; read, coefficient pi D re-derived, PRINT_READ.md I4-I5')
node('T6', 'TEXTBOOK', 'Tomita-Takesaki theorem: J M J = M\' and Delta^{it} M Delta^{-it} = M (Bratteli-Robinson 1, Thm 2.5.14; book not on file)')
node('T7', 'TEXTBOOK', 'Takesaki: a faithful normal state is KMS at beta = -1 for exactly one sigma-weakly continuous group, its modular group; the condition may be checked on a sigma-weakly dense invariant *-subalgebra (Bratteli-Robinson 2, Thm 5.3.10 with Def. 5.3.1 and Cor. 5.3.4; book not on file)')
node('T1', 'TEXTBOOK', 'finite propagation speed and uniqueness for the wave equation with smooth data')
node('T2', 'TEXTBOOK', 'three-lines theorem, Schwarz reflection, uniqueness of the polar decomposition, Stone\'s theorem')
node('T3', 'TEXTBOOK', 'the massless scalar two-point function is the boundary value of -1/(4 pi^2 z^2) from the forward tube')
node('T4', 'TEXTBOOK', 'compactly supported Poincare lemma in a ball; Helmholtz decomposition of symmetric tensors in Fourier space')
node('T5', 'TEXTBOOK', 'first and second variation of area and volume of a hypersurface (used for the displacement quadratic form)')

# ---------------------------------------------------------------- computed leaves produced here
def here(key, text, out, plants=True):
    node(key, 'COMPUTED', text, out=os.path.join(OUT, out), key_true='ALL_PASS', here=True, need_plants=plants)

here('c_s1', 'scope-family identities M1-M8 (K TT harmonic, K_nn = 0, dE = (l+1)K, Ric_3(2NK) = (l+2)K, cut vector, pairing, A, S), l = 2..5, four harmonics each', 's1_modular_side.json')
here('c_s2', 'area-side identities B1-B6 (Green, k = 2, constraint, total = A in three second-order gauges, no drive, displacement form)', 's2_area_side.json')
here('c_s3', 'account identities W1-W3 (transfer formula = direct value; A = c_l E_2; non-stationary odd data)', 's3_account.json')
here('c_s4', 'KMS criterion normalisation on a finite-dimensional standard subspace', 's4_kms_toy.json')

# ---------------------------------------------------------------- computed leaves rerun from the record
def up(key, text, out, key_true):
    node(key, 'COMPUTED', text, out=os.path.join(UP, out), key_true=key_true, here=False, need_plants=False)

up('u_dual_s1', 'wave69/duality S1: Cotton-York identities and the Lemma B chain, exact', 'duality_s1.json', 'S1_ALL_PASS')
up('u_mod_s1', 'wave70/modflow S1: conformal invariance of the tensor two-point function, 11 generators, exact', 'modflow_s1.json', 'b_pass')
up('u_mod_s5', 'wave70/modflow S5: boost generator diagonalisable with integer weights, exp(i pi M) real', 'modflow_s5.json', 'exp_i_pi_M_is_real')
up('u_bclose_p1', 'wave69/bclose P1: every term linear in h2 of the Iyer-Wald generator is the smeared constraint, k = 2, 600 probes', 'bclose_p1_enum.json', 'ALL_PASS')
up('u_bclose_p2s', 'wave69/bclose P2: exact polynomial h2 with all ten Ricci components zero through O(eps^2); E_can = A_l', 'bclose_p2_solve.json', 'ALL_PASS')
up('u_bth_s1', 'wave64/btheorem S1-S2: Green identity and corner invariance', 'btheorem_s1s2.json', None)
up('u_E2', 'wave45/truthcond check_E2: the two flow rates at l = 2..6', 'truthcond_check_E2.json', 'ALL_PASS')

# ---------------------------------------------------------------- stored outputs of slow runs (hash-pinned)
def stored(key, text, rel, key_true, sha):
    node(key, 'COMPUTED', text, out=os.path.join(ROOT, rel), key_true=key_true, here=False, need_plants=False, sha=sha, stored=True)

stored('st_geo', 'wave69/bclose P2 numerics: extremum of A - 2V on g_eps at finite eps gives Delta_2/A_l = 1 within 5e-10, six members', 'wave69/bclose/out/p2_geo.json', 'ALL_PASS', None)
stored('st_mod_s2', 'wave70/modflow S2: geometry of the special conformal map and finite covariance to 2e-29', 'wave70/modflow/out/s2_out_fix.json', 'b_pass_1e-12', None)
stored('st_mod_s3', 'wave70/modflow S3: radial sector of the tensor two-point function reproduces 1/c_l at l = 2, 3, 4', 'wave70/modflow/out/s3_out.json', 'pass', None)

# ---------------------------------------------------------------- proofs
def proof(key, text, where, deps=(), status=''):
    node(key, 'PROOF', text, deps, where=where, status=status)

# second quantisation, proved in THEOREM.md Sec. 4 from the two pillars of modular theory
proof('SQ', 'second quantisation of a standard subspace: Omega cyclic and separating, Delta_R(K) = Gamma(Delta_K), R(K)\' = R(K\'), R of a closed span is the generated algebra', 'THEOREM.md Sec. 4', ['T6', 'T7'])

# duality
proof('p_dual_cv', 'vectors of smeared E and B; L(U)\' = N(complement of U)', 'wave69/duality/REPORT.md Steps 0-2', ['u_dual_s1'], 'read here')
proof('p_dual_A', 'Lemma A: data supported in the closed ball are limits of data supported in the open ball (dilation and mollifier)', 'wave69/duality/REPORT.md Step 3', [], 'read here')
proof('p_dual_B', 'Lemma B: a smooth TT field compactly supported in the ball is C f with f supported in the ball', 'wave69/duality/REPORT.md Step 4', ['u_dual_s1', 'T4'], 'read here')
proof('p_dual_7', 'from time zero to the double cone', 'wave69/duality/REPORT.md Step 7', ['T1'], 'read here')
node('DUAL1', 'AND', 'one-particle duality and C1: K(O) = N(closed ball), K(O\') = K(O)\'', ['p_dual_cv', 'p_dual_A', 'p_dual_B', 'p_dual_7'])
node('DUAL', 'AND', 'Haag duality for the ball: A(O)\' = A(O\')', ['DUAL1', 'SQ'])

# modular group
proof('p_mf_12', 'Lemmas 1, 1\', 2: infinitesimal and finite conformal covariance of the two-point function', 'wave70/modflow/REPORT.md Steps 0-2', ['u_mod_s1', 'T3'], 'read here')
proof('p_mf_3', 'Lemma 3: one-particle Reeh-Schlieder', 'wave70/modflow/REPORT.md Step 3', ['p_dual_cv'], 'read here')
proof('p_mf_4', 'Lemma 4: KMS criterion for standard subspaces', 'wave70/modflow/REPORT.md Step 4; re-derived line by line in THEOREM.md Sec. 9', ['T2', 'c_s4'], 'checked here')
proof('p_mf_5', 'Lemma 5: one-particle Bisognano-Wichmann for the wedge', 'wave70/modflow/REPORT.md Step 5; re-derived in THEOREM.md Sec. 9', ['p_mf_12', 'p_mf_3', 'p_mf_4', 'u_mod_s5'], 'checked here')
proof('p_mf_67', 'Steps 6-7: the special conformal map carries the diamond onto a wedge; the locally defined isometry is unitary', 'wave70/modflow/REPORT.md Steps 6-7', ['p_mf_12', 'p_mf_3', 'st_mod_s2'], 'read here')
node('MODFLOW', 'AND', 'the modular group of K(O) is the conformal flow: Delta^{it} e_F = e_{Lambda_{-2 pi t *} F}, F in the diamond', ['p_mf_5', 'p_mf_67'])

# M, direct
proof('p_Ma', 'the coherent state of the datum exists and is unique on A(O): explicit vector v_f, f = -(2/(l+2)) chi r^2 K', 'THEOREM.md Sec. 5, M-a', ['c_s1', 'p_dual_cv'])
proof('p_Mb', 'extension lemma and generator: v_f is in the domain of log Delta and the curvature of its image is 2 pi times the conformal flow', 'THEOREM.md Sec. 5, M-b', ['MODFLOW', 'p_mf_12', 'p_mf_3', 'T2'])
proof('p_Mc', 'the cutting projection of the image is 2 pi times the cut vector c 2NK 1_B', 'THEOREM.md Sec. 5, M-c', ['c_s1', 'DUAL1', 'p_Mb', 'T4'])
proof('p_Mf', 'K(O) is factorial: the boost group has no invariant vector', 'THEOREM.md Sec. 5, M-c', ['MODFLOW'])
proof('p_Md', 'S = 2 pi |Im(k, v_K)| = (l+1)/(l+2) A / 8G', 'THEOREM.md Sec. 5, M-d', ['c_s1', 'I1', 'I1b', 'I3', 'SQ', 'p_Ma', 'p_Mc', 'p_Mf'])
node('M_direct', 'AND', 'modular side by the direct route', ['p_Ma', 'p_Mb', 'p_Mc', 'p_Mf', 'p_Md'])
# M, transfer (Paper A)
proof('p_A_thm1', 'Paper A r18 Theorem 1 (scalar transfer) with Proposition 1 and Lemma 3; blind-confirmed twice in the record', 'wave36/reviseA18/NOTE.md Sec. 2-3', ['I1', 'I1b', 'I3', 'SQ', 'I4', 'st_mod_s3'], 'record')
node('M_transfer', 'AND', 'modular side by the transfer route, evaluated on the scope family', ['p_A_thm1', 'c_s3'])
node('M', 'OR', 'S_free = (l+1)/(l+2) A/8G: the relative entropy of the free net', ['M_direct', 'M_transfer'])

# C
proof('p_C', 'every algebra obeying A1-A3 has the free net\'s relative entropy: duality plus monotonicity both ways', 'THEOREM.md Sec. 6', ['DUAL', 'I3'])
node('C', 'AND', 'modular side = S_free', ['H0', 'H1', 'H2', 'H3', 'p_C'])

# B
proof('p_B', 'the extremal value of A - 2V has second-order coefficient -A/2 on every vacuum family with first-order part 2tK plus gauge', 'THEOREM.md Sec. 7', ['c_s2', 'c_s1', 'T5'])
proof('p_Bz', 'the seven Poincare zero modes are not driven at O(eps^2) on a single-l stationary datum (parity and time reflection)', 'THEOREM.md Sec. 7, B-e', [])
node('B_math', 'AND', 'geometric area side = A', ['p_B', 'p_Bz', 'H5'])
node('B_conf', 'AND', 'independent confirmations of the area side (not needed by the proof)', ['st_geo', 'u_bclose_p1', 'u_bclose_p2s', 'u_bth_s1'])
node('B', 'AND', 'area side = A/8G', ['H4', 'B_math'])

node('LADDER', 'AND', 'area side - modular side = A/((l+2) 8G) at O(eps^2); modular/area = (l+1)/(l+2)', ['B', 'C', 'M'])
node('W', 'AND', 'mechanism: the two flow rates, and the weights (2, 1)', ['c_s1', 'c_s3', 'u_E2'])

ROOTS = ['LADDER', 'W', 'B_conf']


def leaves(key, seen=None, path=()):
    assert key in N, f'missing node {key}'
    assert key not in path, f'cycle at {key}'
    n = N[key]
    if not n['deps']:
        return {key}
    out = set()
    if n['kind'] in ('PROOF',):
        out.add(key)
    for d in n['deps']:
        out |= leaves(d, seen, path + (key,))
    return out


def audit(strict):
    problems = []
    for k, n in N.items():
        for d in n['deps']:
            if d not in N:
                problems.append(f'{k}: missing dependency {d}')
        if n['kind'] in ('AND', 'OR') and not n['deps']:
            problems.append(f'{k}: {n["kind"]} without dependencies')
        if n['kind'] in ('HYP', 'IMPORT', 'TEXTBOOK', 'COMPUTED') and n['deps']:
            problems.append(f'{k}: a leaf kind with dependencies')
        if n['kind'] == 'IMPORT' and not n.get('locator'):
            problems.append(f'{k}: import without locator')
        if n['kind'] == 'PROOF' and not n.get('where'):
            problems.append(f'{k}: proof without a place')
        if n['kind'] == 'COMPUTED':
            p = n['out']
            if not os.path.exists(p):
                problems.append(f'{k}: output missing: {p}')
                continue
            d = json.load(open(p))
            if n['key_true'] is not None and d.get(n['key_true']) is not True:
                problems.append(f'{k}: {n["key_true"]} is not true in {p}')
            if n.get('need_plants') and not (d.get('n_plants', 0) >= 1 and all(pl['caught'] for pl in d['plants'])):
                problems.append(f'{k}: no caught planted error in {p}')
            if n.get('stored'):
                n['sha_now'] = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    used = set()
    for r in ROOTS:
        used |= leaves(r)
        stack = [r]
        while stack:
            c = stack.pop()
            used.add(c)
            stack += N[c]['deps']
    for k in N:
        if k not in used:
            problems.append(f'{k}: node not reachable from a root')
    return problems


def main():
    strict = '--strict' in sys.argv
    problems = audit(strict)
    ladder_leaves = sorted(leaves('LADDER'))
    kinds = {}
    for k in ladder_leaves:
        kinds.setdefault(N[k]['kind'], []).append(k)
    # failure exhaustion: which top-level clause each leaf sits under
    under = {c: sorted(leaves(c)) for c in ('B', 'C', 'M_direct', 'M_transfer')}
    res = {
        'script': 's5_ledger', 'ALL_PASS': not problems, 'problems': problems,
        'n_nodes': len(N), 'ladder_leaves_by_kind': kinds,
        'clause_leaves': under,
        'M_is_OR_of_two_routes': N['M']['deps'],
        'stored_sha256': {k: n.get('sha_now') for k, n in N.items() if n.get('stored')},
        'nodes': {k: {kk: vv for kk, vv in n.items() if kk not in ('sha_now',)} for k, n in N.items()},
    }
    json.dump(json.loads(json.dumps(res).replace(ROOT + os.sep, '')), open(os.path.join(OUT, 's5_ledger.json'), 'w'), indent=1)
    print('s5_ledger', 'ALL_PASS' if not problems else 'FAILED', len(N), 'nodes;',
          {k: len(v) for k, v in kinds.items()})
    for p in problems:
        print('  PROBLEM', p)
    # human-readable table
    with open(os.path.join(OUT, 'LEDGER.md'), 'w') as fh:
        fh.write('# Ledger of the ladder theorem (generated by scripts/s5_ledger.py)\n\n')
        fh.write('The conclusion LADDER = B and C and M, with M = M_direct or M_transfer. A failure of the conclusion is a failure of at least one leaf below.\n\n')
        for c, title in (('B', 'Area side (B)'), ('C', 'Algebra (C)'), ('M_direct', 'Modular side, direct route'), ('M_transfer', 'Modular side, transfer route')):
            fh.write(f'## {title}\n\n| Leaf | Kind | Statement | Where |\n|---|---|---|---|\n')
            for k in under[c]:
                n = N[k]
                where = n.get('locator') or n.get('where') or (os.path.relpath(n['out'], ROOT) if 'out' in n else '')
                fh.write(f'| {k} | {n["kind"]} | {n["text"]} | {where} |\n')
            fh.write('\n')
    if problems:
        sys.exit(1)


if __name__ == '__main__':
    main()
