#!/usr/bin/env python3
# afm_id.py -- (1) AFM_000 spin desenini bag-bag A/C tipi olarak teshis et
#              (2) strain_* icin FM referansli dogru dE(eps) tablosunu kur
import os, re, glob, sys
import numpy as np
sys.path.insert(0, ".")
from q3_geo import read_qe, neighbours          # test edilmis ayristirici

RY = 13605.693122994

def species_signs(fn):
    t = open(fn, errors="ignore").read()
    L = t.splitlines()
    sp = []
    for i, l in enumerate(L):
        if l.strip().lower().startswith("atomic_species"):
            j = i + 1
            while j < len(L):
                p = L[j].split()
                if len(p) < 3 or not re.match(r"^[A-Za-z]", p[0]): break
                sp.append(p[0]); j += 1
            break
    sm = {}
    for m in re.finditer(r"starting_magnetization\((\d+)\)\s*=\s*(-?[\d.eEdD+]+)", t):
        k = int(m.group(1)) - 1
        if k < len(sp): sm[sp[k]] = float(m.group(2).replace('d', 'e'))
    return sp, sm

print("=== (1) AFM_000 spin deseni ===")
cand = sorted(glob.glob("LNMO_Pilots/AFM_000/*.in"))
for fn in cand:
    sp, sm = species_signs(fn)
    cell, els, X = read_qe(fn)
    labs = [re.sub(r"[^A-Za-z0-9]", "", l) for l in open(fn, errors="ignore").read().split("ATOMIC_POSITIONS")[1].split("\n")[1:] if l.split() and re.match(r"^[A-Za-z]", l.split()[0])]
    labs = labs[:len(els)]
    iNi = [i for i, s in enumerate(els) if s == "Ni"]
    iMn = [i for i, s in enumerate(els) if s == "Mn"]
    iO  = [i for i, s in enumerate(els) if s == "O"]
    sig = []
    for i in range(len(els)):
        lab = labs[i] if i < len(labs) else els[i]
        sig.append(np.sign(sm.get(lab, sm.get(els[i], 0.0))) or 0.0)
    par_p = par_a = perp_p = perp_a = 0
    for i in iNi:
        nb = neighbours(cell, X, i, iO, 6)
        nb.sort(key=lambda t: abs(t[1][2]), reverse=True)
        for k, (d, v, jo) in enumerate(nb):
            Ox = X[i] + v
            mn = neighbours(cell, np.vstack([X, Ox]), len(X), iMn, 1)[0]
            same = sig[i] * sig[mn[2]] > 0
            if k < 2:  perp_p += same; perp_a += not same
            else:      par_p  += same; par_a  += not same
    print("  %s   turler=%s" % (fn, sp))
    print("    duzlem-ici : %d paralel / %d antiparalel" % (par_p, par_a))
    print("    duzlem-disi: %d paralel / %d antiparalel" % (perp_p, perp_a))
    if par_a and not perp_a:   print("    -> C-TIPI (duzlem-ici AF)")
    elif perp_a and not par_a: print("    -> A-TIPI (duzlem-disi AF)")
    elif par_a and perp_a:     print("    -> G-TIPI (hepsi AF)")
    else:                      print("    -> FM / teshis edilemedi")

print()
print("=== (2) strain_* dogru dE(eps) tablosu ===")
def lastE(pat):
    fs = sorted(glob.glob(pat))
    best = None
    for f in fs:
        t = open(f, errors="ignore").read()
        e = re.findall(r"^!\s+total energy\s*=\s*(-?[\d.]+)", t, re.M)
        v = re.findall(r"unit-cell volume\s*=\s*([\d.]+)", t)
        if e: best = (float(e[-1]), float(v[-1])*0.14818471 if v else None, f)
    return best

dirs = sorted(glob.glob("LNMO_Strain_Production/strain_*"))
print("%7s %15s %8s %15s %8s %11s %15s %8s %11s" % (
      "eps","E_FM(Ry)","V_FM","E_FiM(Ry)","V_FiM","dE_FiM","E_AFM(Ry)","V_AFM","dE_AFM"))
for d in dirs:
    eps = float(re.search(r"strain_([+-]?[\d.]+)", d).group(1)) * 100
    fm = lastE(d + "/FM/*.out"); fi = lastE(d + "/FiM/*.out"); af = lastE(d + "/AFM/*.out")
    if not fm: print("%7.1f  *** FM ciktisi yok: %s" % (eps, d)); continue
    def fmt(x):
        if not x: return ("%15s %8s %11s" % ("-", "-", "-"))
        return "%15.8f %8.2f %11.2f" % (x[0], x[1] or 0, (x[0]-fm[0])*RY/2.0)
    print("%7.1f %15.8f %8.2f %s %s" % (eps, fm[0], fm[1] or 0, fmt(fi), fmt(af)))
print()
print("paper-1 hedefleri: FiM 194.9 (-5%) | 192.2 (0%) | 156.8 (+5%) ; AFM 123-134 araligi")
