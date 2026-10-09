#!/usr/bin/env python3
"""cmp_r.py A B [--label L] [--json]: compare two SAIGE step-2 text outputs row by row (by MarkerID).
Reports max relative differences of p.value / BETA / SE / Tstat / var, rows beyond the printed-digit
floor, rows crossing 5e-8 / 1e-5 on one side only, Is.SPA mismatches and Firth-route mismatches
(Firth applied <=> BETA * var != Tstat beyond print rounding)."""
import sys, math, json
def load(p):
    with open(p) as f:
        h = f.readline().rstrip('\n').split('\t')
        rows = {}
        for l in f:
            if not l.strip(): continue
            x = l.rstrip('\n').split('\t'); rows[x[h.index('MarkerID')]] = x
    return h, rows
def fl(x):
    try: return float(x)
    except: return float('nan')
def rel(a, b):
    if a != a and b != b: return 0.0
    if a != a or b != b: return float('inf')
    d = abs(a - b); m = max(abs(a), abs(b))
    return 0.0 if d == 0 else d / m
def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    label = next((sys.argv[i+1] for i, a in enumerate(sys.argv) if a == '--label'), f"{args[0]} vs {args[1]}")
    asjson = '--json' in sys.argv
    ha, A = load(args[0]); hb, B = load(args[1])
    ia = {k: j for j, k in enumerate(ha)}; ib = {k: j for j, k in enumerate(hb)}
    common = [m for m in A if m in B]
    onlyA = len(A) - len(common); onlyB = len(B) - len(common)
    cols = ['p.value', 'BETA', 'SE', 'Tstat', 'var']
    mx = {c: (0.0, None) for c in cols}; beyond = {c: 0 for c in cols}
    floor = {'p.value': 2e-6, 'BETA': 2e-5, 'SE': 2e-5, 'Tstat': 2e-5, 'var': 2e-5}   # printed digits: %.6E / %.6g
    cross = {5e-8: 0, 1e-5: 0}; spa = 0; firth = 0; firthA = firthB = 0; spaA = 0
    worst_cross = []
    for m in common:
        x, y = A[m], B[m]
        for c in cols:
            r = rel(fl(x[ia[c]]), fl(y[ib[c]]))
            if r > mx[c][0]: mx[c] = (r, (m, x[ia[c]], y[ib[c]]))
            if r > floor[c]: beyond[c] += 1
        pa, pb = fl(x[ia['p.value']]), fl(y[ib['p.value']])
        for th in cross:
            if (pa < th) != (pb < th):
                cross[th] += 1
                if len(worst_cross) < 5: worst_cross.append((m, th, x[ia['p.value']], y[ib['p.value']]))
        if 'Is.SPA' in ia and 'Is.SPA' in ib:
            if x[ia['Is.SPA']] != y[ib['Is.SPA']]: spa += 1
            if x[ia['Is.SPA']] == 'true': spaA += 1
        def fr(z, i):
            b, v, t = fl(z[i['BETA']]), fl(z[i['var']]), fl(z[i['Tstat']])
            if not (v > 0) or t != t or b != b: return False
            return abs(b * v - t) > 1e-4 * max(abs(t), 1e-300)
        fa, fb = fr(x, ia), fr(y, ib)
        firthA += fa; firthB += fb
        if fa != fb: firth += 1
    out = {'label': label, 'rows': len(common), 'onlyA': onlyA, 'onlyB': onlyB,
           'maxrel': {c: mx[c][0] for c in cols}, 'worst': {c: mx[c][1] for c in cols},
           'beyondPrint': beyond, 'cross5e8': cross[5e-8], 'cross1e5': cross[1e-5],
           'isSPAmismatch': spa, 'nSPA_A': spaA, 'firthMismatch': firth, 'firthA': firthA, 'firthB': firthB,
           'worstCross': worst_cross}
    if asjson: print(json.dumps(out)); return
    print(f"{label}: rows {len(common)} (onlyA {onlyA}, onlyB {onlyB})")
    for c in cols:
        w = mx[c][1]
        print(f"  {c:8s} maxrel {mx[c][0]:.3g}  beyond-print {beyond[c]:6d}" + (f"  worst {w[0]} {w[1]} vs {w[2]}" if w else ""))
    print(f"  cross 5e-8: {cross[5e-8]}  cross 1e-5: {cross[1e-5]}  Is.SPA mismatch: {spa} (SPA rows {spaA})"
          f"  Firth-route mismatch: {firth} (Firth rows {firthA} / {firthB})")
    for w in worst_cross: print(f"    cross {w[1]}: {w[0]} {w[2]} vs {w[3]}")
if __name__ == '__main__': main()
