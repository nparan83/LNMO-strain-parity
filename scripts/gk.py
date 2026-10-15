#!/usr/bin/env python3
# gk.py : pw.x ciktilari -> spin-cozumlu bant araliklari + kenarlarin k-noktasi (kristal)
import re, sys
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
        if "highest occupied, lowest unoccupied" in L:
            v = [float(x) for x in F.findall(L)]; ef = 0.5*(v[0]+v[1]); mode = None; continue
        if mode == "e":
            if re.search(r"[A-Za-z]", L): mode = None; continue
            cur[1].extend(float(x) for x in F.findall(L))
    return A, eig, ef
def cr(A, k): return "(%+.2f %+.2f %+.2f)" % tuple(sum(k[j]*A[i][j] for j in range(3)) for i in (1, 2, 3))
print("%-58s %4s %8s %8s %8s  %-20s %-20s" % ("dosya", "nk", "Eg_up", "Eg_dn", "Eg_fund", "VB_up k", "CB_up k"))
for fn in sys.argv[1:]:
    try: A, eig, ef = parse(fn)
    except Exception as ex: print(fn, "HATA", ex); continue
    if not eig["up"] or ef is None or len(A) < 3: print("%-58s  ozdeger/E_F/eksen yok" % fn[-58:]); continue
    r = {}
    for s in ("up", "dn"):
        vb = max(((x, k) for k, b in eig[s] for x in b if x <= ef), key=lambda t: t[0])
        cb = min(((x, k) for k, b in eig[s] for x in b if x > ef),  key=lambda t: t[0])
        r[s] = (vb, cb)
    g = min(r["up"][1][0], r["dn"][1][0]) - max(r["up"][0][0], r["dn"][0][0])
    print("%-58s %4d %8.4f %8.4f %8.4f  %-20s %-20s" % (fn[-58:], len(eig["up"]),
          r["up"][1][0]-r["up"][0][0], r["dn"][1][0]-r["dn"][0][0], g, cr(A, r["up"][0][1]), cr(A, r["up"][1][1])))
