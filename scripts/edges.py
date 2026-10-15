#!/usr/bin/env python3
# edges.py : nscf_FM.out -> spin-cozumlu bant kenarlari + temel bant araligi
import re, glob, os, sys
R = sys.argv[1]
def epsof(p):
    m = re.search(r"strain_([-+]?\d*\.?\d+)", p); v = float(m.group(1))
    return 100*v if abs(v) < 0.2 else v
fl = sorted(glob.glob(os.path.join(R, "strain_*", "FM", "nscf_FM.out")), key=epsof)
if not fl: sys.exit("HATA: %s/strain_*/FM/nscf_FM.out bulunamadi" % R)
REF = {-5: (1.0763, 3.1580), 0: (1.2950, 3.2321), 5: (1.0836, 2.8969)}
FLT = re.compile(r"-?\d+\.\d+")
def parse(fn):
    eig = {"up": [], "dn": []}; spin = None; mode = None; ef = {}
    for L in open(fn, errors="ignore"):
        if "SPIN UP" in L:   spin, mode = "up", None; continue
        if "SPIN DOWN" in L: spin, mode = "dn", None; continue
        if "bands (ev)" in L and spin:
            cur = []; eig[spin].append(cur); mode = "e"; continue
        if "Fermi energies are" in L:
            v = FLT.findall(L); ef = {"up": float(v[0]), "dn": float(v[1])}; mode = None; continue
        if "the Fermi energy is" in L:
            v = float(FLT.findall(L)[0]); ef = {"up": v, "dn": v}; mode = None; continue
        if "highest occupied, lowest unoccupied" in L:
            v = [float(x) for x in FLT.findall(L)]; m = 0.5*(v[0]+v[1]); ef = {"up": m, "dn": m}; mode = None; continue
        if mode == "e":
            if re.search(r"[A-Za-z]", L): mode = None; continue
            cur.extend(float(x) for x in FLT.findall(L))
    return eig, ef
print("%6s %4s %4s %9s %9s %9s %9s %8s %8s %8s  %s" % ("eps", "n_up", "n_dn", "VB_up", "CB_up",
      "VB_dn", "CB_dn", "Eg_up", "Eg_dn", "Eg_fund", "kenar(VB/CB)"))
bad = 0
for fn in fl:
    e = epsof(fn); eig, ef = parse(fn)
    if not eig["up"] or not eig["dn"] or not ef:
        print("%6.1f  HATA: ozdeger/E_F okunamadi -> %s" % (e, fn)); bad += 1; continue
    r = {}
    for s in ("up", "dn"):
        a = [x for k in eig[s] for x in k]
        r[s] = (max(x for x in a if x <= ef[s]), min(x for x in a if x > ef[s]),
                sum(1 for x in eig[s][0] if x <= ef[s]), len(eig[s]))
    vb = max(r["up"][0], r["dn"][0]); cb = min(r["up"][1], r["dn"][1])
    tag = ("up" if r["up"][0] >= r["dn"][0] else "dn") + "/" + ("up" if r["up"][1] <= r["dn"][1] else "dn")
    gu, gd = r["up"][1]-r["up"][0], r["dn"][1]-r["dn"][0]
    print("%+6.1f %4d %4d %9.4f %9.4f %9.4f %9.4f %8.4f %8.4f %8.4f  %s  (nk=%d)" % (e, r["up"][2], r["dn"][2],
          r["up"][0], r["up"][1], r["dn"][0], r["dn"][1], gu, gd, cb-vb, tag, r["up"][3]))
    k = int(round(e))
    if k in REF and abs(e-k) < 1e-6:
        ok = abs(gu-REF[k][0]) < 5e-4 and abs(gd-REF[k][1]) < 5e-4
        bad += (not ok); print("        referans gap.dat (%+d%%): %s" % (k, "OK" if ok else "<-- SAPMA"))
print("dogrulama: %d sorun" % bad)
