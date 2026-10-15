#!/usr/bin/env python3
# kedge.py : bant kenarlarinin k-noktasi (kristal koord.) -> mesh'e bagimlilik testi
import re, sys, os
R = sys.argv[1]; EPS = ["-0.05", "-0.01", "0.00", "0.01", "0.05"]
F = re.compile(r"-?\d+\.\d+")
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
        if "the Fermi energy is" in L: ef = float(F.findall(L)[0]); mode = None; continue
        if mode == "e":
            if re.search(r"[A-Za-z]", L): mode = None; continue
            cur[1].extend(float(x) for x in F.findall(L))
    return A, eig, ef
def cryst(A, k): return tuple(sum(k[j]*A[i][j] for j in range(3)) for i in (1, 2, 3))
print("%6s %3s %5s %9s %-26s %9s %-26s" % ("eps", "spn", "nbnd", "VB", "k_VB (kristal)", "CB", "k_CB (kristal)"))
for e in EPS:
    fn = os.path.join(R, "strain_" + e, "FM", "nscf_FM.out")
    if not os.path.exists(fn): print(e, "yok:", fn); continue
    A, eig, ef = parse(fn)
    for s in ("up", "dn"):
        vb = max(((x, k) for k, b in eig[s] for x in b if x <= ef), key=lambda t: t[0])
        cb = min(((x, k) for k, b in eig[s] for x in b if x > ef),  key=lambda t: t[0])
        f = lambda k: "(%+.4f %+.4f %+.4f)" % cryst(A, k)
        print("%6s %3s %5d %9.4f %-26s %9.4f %-26s" % (e, s, len(eig[s][0][1]), vb[0], f(vb[1]), cb[0], f(cb[1])))
