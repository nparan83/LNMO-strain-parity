#!/usr/bin/env python3
# extract_p1.py - Paper-1 (La2NiMnO6 under biaxial strain): rebuilds the numbers of
# Sec. 3.1-3.4, Tables 1 and 3 and Figs. 2-4 directly from the repository files
# (inputs/*.in and gzip-compressed outputs/*.out.gz) and checks them against the manuscript.
# Usage: python scripts/extract_p1.py            (run from the repository root)
# Writes: data/energies.dat, data/cperp.dat, data/lowdin.dat, data/gap_check.dat
import glob, gzip, os, re, sys
import numpy as np

RY = 13605.693122994          # meV per Ry
BOHR = 0.529177210903         # Angstrom
Z = 2                         # formula units per 20-atom cell
NOCC = {"up": 88, "dn": 78}   # occupied bands per spin channel (fixed at every strain)
NUM = re.compile(r"-?\d+\.\d+")

def rd(p):
    return (gzip.open(p, "rt", errors="replace") if p.endswith(".gz") else open(p, errors="replace")).read()

def eps_of(d):
    return float(os.path.basename(d).split("_", 1)[1])

def first(paths):
    for p in paths:
        if os.path.exists(p): return p
    return None

def energy(p):
    t = rd(p); e = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", t, re.M)
    m = re.findall(r"total magnetization\s+=\s+(-?\d+\.\d+)", t)
    return (float(e[-1]) if e else np.nan), (float(m[-1]) if m else np.nan)

def cell_and_species(p):
    t = rd(p)
    hdr = re.search(r"^\s*CELL_PARAMETERS(.*)$", t, re.I | re.M)
    rows = t[hdr.end():].strip().splitlines()[:3]
    A = np.array([[float(x) for x in r.split()[:3]] for r in rows])
    unit = re.sub(r"[^a-z]", "", hdr.group(1).lower())
    if unit.startswith("bohr"): A *= BOHR
    elif unit.startswith("alat") or unit == "":
        c1 = re.search(r"celldm\(1\)\s*=\s*([\d.eEdD+-]+)", t)
        aa = re.search(r"\bA\s*=\s*([\d.eEdD+-]+)", t)
        A *= float(c1.group(1).replace("d", "e").replace("D", "e")) * BOHR if c1 else float(aa.group(1))
    nat = int(re.search(r"nat\s*=\s*(\d+)", t).group(1))
    pos = t[re.search(r"ATOMIC_POSITIONS", t).end():].splitlines()[1:nat + 1]
    sp = [re.sub(r"[^A-Za-z]", "", l.split()[0]) for l in pos]
    return A, sp

def lowdin(t):
    """Per-atom Lowdin charge and polarization from a projwfc.x output (last block)."""
    blk = t[t.rfind("Lowdin Charges"):]
    out = {}
    lines = blk.splitlines()
    for i, l in enumerate(lines):
        m = re.match(r"\s*Atom #\s*(\d+):\s*total charge\s*=\s*(-?\d+\.\d+)", l)
        if not m: continue
        n, q, pol, up, dn = int(m.group(1)), float(m.group(2)), None, None, None
        for l2 in lines[i + 1:i + 6]:
            if "Atom #" in l2: break
            if "polarization" in l2: pol = float(NUM.findall(l2)[0])
            if "spin up" in l2: up = float(NUM.findall(l2)[0])
            if "spin down" in l2: dn = float(NUM.findall(l2)[0])
        if pol is None and up is not None and dn is not None: pol = up - dn
        if n not in out: out[n] = (q, pol)
    return out

def gaps(t):
    """Spin-resolved gap with fixed band indices; eigenvalues from the last 'bands (ev)' listing."""
    res = {}
    for spin, tag in (("up", "SPIN UP"), ("dn", "SPIN DOWN")):
        i0 = t.rfind("------ " + tag)
        if i0 < 0: return None
        seg = t[i0:]
        if spin == "up": seg = seg[:seg.find("------ SPIN DOWN")]
        end = seg.find("the Fermi energy"); seg = seg[:end] if end > 0 else seg
        vbm, cbm, n = -1e9, 1e9, NOCC[spin]
        for b in re.split(r"\n\s*k =", seg)[1:]:
            body = b.split("bands (ev):", 1)
            if len(body) < 2: continue
            body = body[1].split("occupation numbers")[0]
            ev = [float(x) for x in NUM.findall(body)]
            if len(ev) > n: vbm, cbm = max(vbm, ev[n - 1]), min(cbm, ev[n])
        res[spin] = (cbm - vbm, vbm, cbm)
    return res

root = os.getcwd()
dirs = sorted(glob.glob(os.path.join(root, "inputs", "strain_*")), key=eps_of)
if len(dirs) != 11: sys.exit("expected 11 strain directories under inputs/, found %d" % len(dirs))
os.makedirs("data", exist_ok=True)
E, C, L, G = [], [], [], []
for d in dirs:
    e = eps_of(d); s = os.path.basename(d); od = os.path.join(root, "outputs", s)
    fFM = first([os.path.join(od, "FM", "scf_FM.out.gz")])
    fAF = first([os.path.join(od, "AFM", x) for x in ("scf_final_AFM.out.gz", "scf_AFM.out.gz")])
    fFi = first([os.path.join(od, "FiM", x) for x in ("scf_final_FiM.out.gz", "scf_FiM.out.gz")])
    eFM, mFM = energy(fFM) if fFM else (np.nan, np.nan)
    eAF, _ = energy(fAF) if fAF else (np.nan, np.nan)
    eFi, mFi = energy(fFi) if fFi else (np.nan, np.nan)
    E.append([e, eFM, eAF, eFi, mFM, mFi])
    print("eps %+.2f  FM %s | AFM %s | FiM %s" % (e, fFM and os.path.basename(fFM), fAF and os.path.basename(fAF), fFi and os.path.basename(fFi)))
    A, sp = cell_and_species(os.path.join(d, "FM", "scf_FM.in"))
    V = abs(np.linalg.det(A)); cp = V / np.linalg.norm(np.cross(A[0], A[1]))
    C.append([e] + list(np.linalg.norm(A, axis=1)) + [V, cp])
    lw = None
    for p in sorted(glob.glob(os.path.join(od, "FM", "*.out.gz"))):
        t = rd(p)
        if "Lowdin Charges" in t: lw = lowdin(t); break
    if lw:
        mn = [lw[i + 1] for i, x in enumerate(sp) if x == "Mn" and i + 1 in lw]
        ni = [lw[i + 1] for i, x in enumerate(sp) if x == "Ni" and i + 1 in lw]
        L.append([e, mFM, np.mean([q for q, _ in mn]), np.mean([m for _, m in mn]),
                  np.mean([q for q, _ in ni]), np.mean([m for _, m in ni])])
    fN = first([os.path.join(od, "FM", "nscf_FM.out.gz")])
    g = gaps(rd(fN)) if fN else None
    if g: G.append([e, g["up"][0], g["dn"][0], g["up"][1], g["up"][2], g["dn"][1], g["dn"][2]])

E = np.array(E); e = E[:, 0]; i0 = int(np.argmin(abs(e)))
dAF = (E[:, 2] - E[:, 1]) / Z * RY; dFi = (E[:, 3] - E[:, 1]) / Z * RY
rFM = (E[:, 1] - E[i0, 1]) / Z * RY; rAF = (E[:, 2] - E[i0, 1]) / Z * RY
np.savetxt("data/energies.dat", np.c_[e, E[:, 1:4], dAF, dFi, rFM, rAF, E[:, 4]], fmt="%14.8f",
           header="eps E_FM_Ry E_AFM_Ry E_FiM_Ry dE_AFM-FM_meV/fu dE_FiM-FM_meV/fu E_FM-E_FM(0)_meV/fu E_AFM-E_FM(0)_meV/fu Mtot_FM_muB/cell")
np.savetxt("data/cperp.dat", np.array(C), fmt="%12.6f", header="eps |a| |b| |c| V(A^3) c_perp=V/|axb|(A)")
if L: np.savetxt("data/lowdin.dat", np.array(L), fmt="%10.4f", header="eps Mtot Mn_pop Mn_mom Ni_pop Ni_mom (Lowdin, projwfc.x)")
if G: np.savetxt("data/gap_check.dat", np.array(G), fmt="%10.5f", header="eps Eg_up Eg_dn VBM_up CBM_up VBM_dn CBM_dn (eV; fixed n_up=88, n_dn=78)")

# ---------------- checks against the manuscript ----------------
Cp = np.array(C)[:, 5]
T3u = [1.0763, 1.1430, 1.1942, 1.2406, 1.2764, 1.2950, 1.2953, 1.2526, 1.1912, 1.1325, 1.0836]
T3d = [3.1580, 3.2022, 3.2304, 3.2437, 3.2439, 3.2321, 3.1812, 3.1164, 3.0474, 2.9753, 2.8969]
T1 = {"Mn_pop": [13.299, 13.298, 13.297, 13.296, 13.295, 13.295, 13.295, 13.296, 13.296, 13.296, 13.298],
      "Mn_mom": [3.211, 3.227, 3.245, 3.263, 3.282, 3.302, 3.324, 3.347, 3.370, 3.393, 3.414],
      "Ni_pop": [9.588, 9.588, 9.587, 9.586, 9.585, 9.583, 9.580, 9.577, 9.574, 9.570, 9.566],
      "Ni_mom": [1.649, 1.654, 1.658, 1.663, 1.668, 1.674, 1.680, 1.686, 1.692, 1.698, 1.703]}
chk = [("dE(AFM-FM) at eps=0, meV/f.u.", dAF[i0], 133.83, 0.01),
       ("min dE(AFM-FM), meV/f.u. ('at least 123')", np.nanmin(dAF), 123.0, 0.05),
       ("dE(FiM-FM) at -5 %", dFi[0], 194.9, 0.05), ("dE(FiM-FM) at 0", dFi[i0], 192.2, 0.05),
       ("dE(FiM-FM) at +5 %", dFi[-1], 156.8, 0.05),
       ("c_perp(0), A", Cp[i0], 7.84554, 0.0005), ("c_perp(-5 %), A", Cp[0], 8.20899, 0.0005), ("c_perp(+5 %), A", Cp[-1], 7.58622, 0.0005),
       ("c_perp(-5)/c_perp(0)-1, %", 100 * (Cp[0] / Cp[i0] - 1), 4.6, 0.05),
       ("c_perp(+5)/c_perp(0)-1, %", 100 * (Cp[-1] / Cp[i0] - 1), -3.3, 0.05),
       ("M_tot (FM), muB/cell, max |dev| from 10", np.nanmax(abs(E[:, 4] - 10)), 0.0, 0.005)]
if L:
    La = np.array(L)
    for k, col in (("Mn_pop", 2), ("Mn_mom", 3), ("Ni_pop", 4), ("Ni_mom", 5)):
        chk.append(("Table 1 %s, max |dev|" % k, np.max(abs(La[:, col] - T1[k])), 0.0, 0.0006))
if G:
    Ga = np.array(G)
    chk.append(("Table 3 Eg_up, max |dev| (eV)", np.max(abs(Ga[:, 1] - T3u)), 0.0, 0.00006))
    chk.append(("Table 3 Eg_dn, max |dev| (eV)", np.max(abs(Ga[:, 2] - T3d)), 0.0, 0.00006))
print("\n=== check against the manuscript ===")
bad = 0
for nm, got, want, tol in chk:
    ok = abs(got - want) <= tol + 1e-12; bad += (not ok)
    print("  %-46s computed=%10.4f  manuscript=%9.4f  %s" % (nm, got, want, "OK" if ok else "<-- MISMATCH"))
print("Lowdin blocks found: %d/11   NSCF gap blocks found: %d/11" % (len(L), len(G)))
print("check: %d/%d mismatches" % (bad, len(chk)))
