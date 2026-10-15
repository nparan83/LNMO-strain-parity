#!/usr/bin/env python3
# gk2.py : E_F'den BAGIMSIZ bant araligi — sabit bant sayisi (n_up=88, n_dn=78)
import re, sys
F = re.compile(r"-?\d+\.\d+"); NOCC = {"up": 88, "dn": 78}
def parse(fn):
    A = {}; eig = {"up": [], "dn": []}; spin = None; mode = None; ef = None
    for L in open(fn, errors="ignore"):
        m = re.match(r"\s*a\((\d)\) = \((.*)\)", L)
        if m: A[int(m.group(1))] = [float(x) for x in F.findall(m.group(2))]; continue
        if "SPIN UP" in L:   spin, mode = "up", None; continue
        if "SPIN DOWN" in L: spin, mode = "dn", None; continue
        if "bands (ev)" in L and spin:
            k = [float(x) for x in F.findall(L.split("=", 1)[1].split("(")[0])]
            cur = [k, []]; eig[spin].append(cur); mode = "e"; continue
        if "Fermi energ" in L: ef = L.strip(); mode = None; continue
        if mode == "e":
            if re.search(r"[A-Za-z]", L): mode = None; continue
            cur[1].extend(float(x) for x in F.findall(L))
    return A, eig, ef
def cr(A, k): return "(%+.3f %+.3f %+.3f)" % tuple(sum(k[j]*A[i][j] for j in range(3)) for i in (1, 2, 3))
for fn in sys.argv[1:]:
    A, eig, ef = parse(fn)
    print("== %s   nk=%d   [%s]" % (fn, len(eig["up"]), ef))
    for s in ("up", "dn"):
        n = NOCC[s]
        vb = max(((b[n-1], k) for k, b in eig[s]), key=lambda t: t[0])
        cb = min(((b[n], k) for k, b in eig[s]), key=lambda t: t[0])
        print("   %s  VB(band %d)=%9.4f @ %-24s CB(band %d)=%9.4f @ %-24s Eg=%.4f" % (s, n, vb[0], cr(A, vb[1]), n+1, cb[0], cr(A, cb[1]), cb[0]-vb[0]))
