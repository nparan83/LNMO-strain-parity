#!/usr/bin/env python3
# q3lib.py -- ortak QE ayristirici (import-guvenli, govde kodu YOK)
import re
import numpy as np

BOHR = 0.529177210903
SH = np.array([[i,j,k] for i in (-1,0,1) for j in (-1,0,1) for k in (-1,0,1)], float)

def read_qe(fn):
    txt = open(fn, errors="ignore").read()
    low = txt.lower()
    alat = None
    m = re.search(r"celldm\(1\)\s*=\s*([\deEdD.+-]+)", low)
    if m: alat = float(m.group(1).replace('d','e')) * BOHR
    m = re.search(r"^\s*a\s*=\s*([\deEdD.+-]+)", low, re.M)
    if m and alat is None: alat = float(m.group(1).replace('d','e'))
    lines = txt.splitlines()
    cell = None; cunit = "alat"; pos = []; punit = "alat"; i = 0
    while i < len(lines):
        k = lines[i].strip().lower()
        if k.startswith("cell_parameters"):
            cunit = (re.sub(r"[{}()]", " ", k).split() + ["alat"])[1]
            cell = np.array([[float(x) for x in lines[i+1+j].split()[:3]] for j in range(3)])
            i += 4; continue
        if k.startswith("atomic_positions"):
            punit = (re.sub(r"[{}()]", " ", k).split() + ["alat"])[1]
            j = i + 1
            while j < len(lines):
                p = lines[j].split()
                if len(p) < 4 or not re.match(r"^[A-Za-z]", p[0]): break
                try: xyz = [float(p[1]), float(p[2]), float(p[3])]
                except ValueError: break
                pos.append((p[0], xyz)); j += 1
            i = j; continue
        i += 1
    if cell is None: raise RuntimeError("CELL_PARAMETERS yok: " + fn)
    if   cunit.startswith("bohr"): cell *= BOHR
    elif cunit.startswith("alat"):
        if alat is None: raise RuntimeError("CELL_PARAMETERS alat ama celldm(1)/A yok: " + fn)
        cell *= alat
    X = np.array([p[1] for p in pos], float)
    if   punit.startswith("crystal"): X = X @ cell
    elif punit.startswith("bohr"):    X = X * BOHR
    elif punit.startswith("alat"):    X = X * alat
    els = [re.sub(r"[^A-Za-z]", "", p[0]) for p in pos]
    return cell, els, X

def raw_labels(fn):
    L = open(fn, errors="ignore").read().splitlines()
    out = []
    for i, l in enumerate(L):
        if l.strip().lower().startswith("atomic_positions"):
            j = i + 1
            while j < len(L):
                p = L[j].split()
                if len(p) < 4 or not re.match(r"^[A-Za-z]", p[0]): break
                try: float(p[1]); float(p[2]); float(p[3])
                except ValueError: break
                out.append(p[0]); j += 1
            break
    return out

def species(fn):
    L = open(fn, errors="ignore").read().splitlines()
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
    for m in re.finditer(r"starting_magnetization\((\d+)\)\s*=\s*(-?[\d.eEdD+]+)",
                         open(fn, errors="ignore").read()):
        k = int(m.group(1)) - 1
        if k < len(sp): sm[sp[k]] = float(m.group(2).replace('d','e').replace('D','e'))
    return sp, sm

def neighbours(cell, X, ic, idx_pool, n):
    T = SH @ cell
    out = []
    for ja in idx_pool:
        V = X[ja] + T - X[ic]
        d = np.linalg.norm(V, axis=1)
        m = int(np.argmin(d))
        out.append((d[m], V[m], ja))
    out.sort(key=lambda t: t[0])
    return out[:n]
