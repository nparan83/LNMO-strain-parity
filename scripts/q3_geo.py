#!/usr/bin/env python3
# q3_geo.py -- NiO6 / MnO6 Jahn-Teller Q3 (ve Q2) genlikleri, paper-1 konvansiyonunda
#   Q3 = (2/sqrt(6)) * ( d_apikal - d_ekvatoryal )      [Angstrom]
#   Q2 = (1/sqrt(2)) * ( d_eq,a - d_eq,b )
# Kaynak: LNMO_Strain_Production/strain_*/FM/scf_FM.in
# Cikti : q3_geo.dat  +  ekrana ~15 satir ozet
import os, re, sys, glob
import numpy as np

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/arf/scratch/nparan/lnmo/LNMO_Strain_Production"
PAT  = os.path.join(ROOT, "strain_*", "FM", "scf_FM.in")
BOHR = 0.529177210903

def eps_of(path):
    m = re.search(r"strain_([+-]?[pm]?\d+(?:[._]\d+)?)", path)
    if not m: return None
    s = m.group(1).replace("_", ".").replace("p", "+").replace("m", "-")
    try: return float(s)
    except ValueError: return None

def read_qe(fn):
    txt = open(fn, errors="ignore").read()
    low = txt.lower()
    alat = None
    m = re.search(r"celldm\(1\)\s*=\s*([\deEdD.+-]+)", low)
    if m: alat = float(m.group(1).replace('d','e')) * BOHR
    m = re.search(r"^\s*a\s*=\s*([\deEdD.+-]+)", low, re.M)
    if m and alat is None: alat = float(m.group(1).replace('d','e'))
    lines = txt.splitlines()
    cell = None; cunit = "alat"
    pos = []; punit = "alat"
    i = 0
    while i < len(lines):
        L = lines[i].strip()
        k = L.lower()
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
    labs = [p[0] for p in pos]; X = np.array([p[1] for p in pos], float)
    if   punit.startswith("crystal"): X = X @ cell
    elif punit.startswith("bohr"):    X = X * BOHR
    elif punit.startswith("alat"):    X = X * alat
    els = [re.sub(r"[^A-Za-z]", "", l) for l in labs]
    return cell, els, X

SH = np.array([[i,j,k] for i in (-1,0,1) for j in (-1,0,1) for k in (-1,0,1)], float)

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

K3 = 2.0/np.sqrt(6.0); K2 = 1.0/np.sqrt(2.0)

def octa(cell, X, ic, iO):
    nb = neighbours(cell, X, ic, iO, 6)
    nb.sort(key=lambda t: abs(t[1][2]), reverse=True)      # |dz| buyuk 2 tanesi apikal
    ap, eq = nb[:2], nb[2:]
    d_ap = np.mean([t[0] for t in ap]); d_eq = np.mean([t[0] for t in eq])
    u = [t[1]/t[0] for t in eq]; used = set(); pairs = []
    for a in range(4):
        if a in used: continue
        b = min([c for c in range(4) if c != a and c not in used], key=lambda c: u[a] @ u[c])
        used |= {a, b}; pairs.append((a, b))
    da = np.mean([eq[p][0] for p in pairs[0]]); db = np.mean([eq[p][0] for p in pairs[1]])
    ra = abs(eq[pairs[0][0]][1][0]); rb = abs(eq[pairs[1][0]][1][0])
    if rb > ra: da, db = db, da
    return d_ap, d_eq, K3*(d_ap-d_eq), K2*(da-db), ap, eq

rows = []
files = sorted(glob.glob(PAT))
if not files: sys.exit("HATA: dosya bulunamadi -> " + PAT)
bad_name = [f for f in files if eps_of(f) is None]
if bad_name:
    print("HATA: dizin adindan eps okunamadi. Ornekler:")
    for f in bad_name[:5]: print("   ", f)
    sys.exit("  -> eps_of() regexini dizin adina gore duzelt.")
for fn in files:
    e = eps_of(fn)
    cell, els, X = read_qe(fn)
    iNi = [i for i,s in enumerate(els) if s == "Ni"]
    iMn = [i for i,s in enumerate(els) if s == "Mn"]
    iO  = [i for i,s in enumerate(els) if s == "O"]
    r = {"eps": e, "nNi": len(iNi), "nMn": len(iMn), "nO": len(iO), "f": fn}
    for tag, idx in (("Ni", iNi), ("Mn", iMn)):
        A = np.array([octa(cell, X, i, iO)[:4] for i in idx])
        r[tag] = A.mean(axis=0); r[tag+"_sd"] = A.std(axis=0)
    tp, tq = [], []
    for i in iNi:
        d_ap, d_eq, q3, q2, ap, eq = octa(cell, X, i, iO)
        for grp, box in ((eq, tp), (ap, tq)):
            for dist, vec, jo in grp:
                Ox = X[i] + vec
                mn = neighbours(cell, np.vstack([X, Ox]), len(X), iMn, 1)[0]
                v1 = -vec; v2 = mn[1]
                c = v1 @ v2 / (np.linalg.norm(v1)*np.linalg.norm(v2))
                box.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    r["th_par"] = float(np.mean(tp)); r["th_perp"] = float(np.mean(tq))
    r["nbond"] = (len(tp), len(tq))
    rows.append(r)
rows.sort(key=lambda r: (r["eps"] if r["eps"] is not None else 0))
if max(abs(r["eps"]) for r in rows) <= 1.0:   # kesir olarak okunmus -> yuzde
    for r in rows: r["eps"] *= 100.0
    print("NOT: dizin adlarindaki eps kesir formatinda; yuzdeye cevrildi (x100).")

with open("q3_geo.dat", "w") as f:
    f.write("# eps  d_eq_Ni  d_ap_Ni  Q3_Ni  Q2_Ni  d_eq_Mn  d_ap_Mn  Q3_Mn  Q2_Mn  th_par  th_perp  nNi nMn nO n_par n_perp  file\n")
    f.write("# Q3 = (2/sqrt6)(d_ap - d_eq)  [paper-1 konvansiyonu, Angstrom]\n")
    for r in rows:
        f.write("%+5.2f %8.4f %8.4f %+9.5f %+9.5f %8.4f %8.4f %+9.5f %+9.5f %9.4f %9.4f %3d %3d %3d %3d %3d  %s\n" % (
            r["eps"], r["Ni"][1], r["Ni"][0], r["Ni"][2], r["Ni"][3],
            r["Mn"][1], r["Mn"][0], r["Mn"][2], r["Mn"][3],
            r["th_par"], r["th_perp"], r["nNi"], r["nMn"], r["nO"],
            r["nbond"][0], r["nbond"][1], r["f"]))

print("q3_geo.dat yazildi:", len(rows), "nokta")
bad = [r for r in rows if (r["nNi"],r["nMn"],r["nO"]) != (rows[0]["nNi"],rows[0]["nMn"],rows[0]["nO"])
       or r["nbond"] != rows[0]["nbond"]]
print("atom/bag sayimi:", rows[0]["nNi"], "Ni |", rows[0]["nMn"], "Mn |", rows[0]["nO"], "O |",
      rows[0]["nbond"][0], "duzlem-ici +", rows[0]["nbond"][1], "duzlem-disi bag",
      "| TUM NOKTALARDA AYNI" if not bad else "| *** TUTARSIZ: %d nokta ***" % len(bad))
sd = max(max(r["Ni_sd"][2], r["Mn_sd"][2]) for r in rows)
print("oktahedralar arasi en buyuk Q3 sacilimi: %.5f A  (%s)" % (sd, "esdeger" if sd < 5e-3 else "*** ESDEGER DEGIL, BAK ***"))
E  = np.array([r["eps"] for r in rows]); Q3 = np.array([r["Ni"][2] for r in rows])
Q3m= np.array([r["Mn"][2] for r in rows]); Q2 = np.array([r["Ni"][3] for r in rows])
c  = np.polyfit(E, Q3, 1); rr = np.corrcoef(E, Q3)[0,1]
cm = np.polyfit(E, Q3m,1)
print("Q3_Ni(eps) = %+.6f*eps %+.6f   r=%.6f   sifir gecis eps=%+.4f%%" % (c[0], c[1], rr, -c[1]/c[0]))
print("Q3_Mn(eps) = %+.6f*eps %+.6f   |egim orani Ni/Mn| = %.2f" % (cm[0], cm[1], abs(c[0]/cm[0])))
print("Q3_Ni araligi: %+.5f (eps=%+.0f) ... %+.5f (eps=%+.0f)   |Q2_Ni|max = %.5f" % (
      Q3[0], E[0], Q3[-1], E[-1], np.abs(Q2).max()))
print("th_par : %.3f ... %.3f   th_perp: %.3f ... %.3f" % (
      min(r["th_par"] for r in rows), max(r["th_par"] for r in rows),
      min(r["th_perp"] for r in rows), max(r["th_perp"] for r in rows)))

REF = {  # durum_La2NiMnO6_J.md  sec.7e  (geo2.py -> geo_theta.dat)
 -5:(1.9616,2.1697,153.205,163.799), -4:(1.9821,2.1465,154.137,162.099),
 -3:(2.0035,2.1253,154.887,160.304), -2:(2.0265,2.1061,155.408,158.532),
 -1:(2.0510,2.0886,155.691,156.815),  0:(2.0768,2.0727,155.780,155.156),
  1:(2.1042,2.0581,155.661,153.587),  2:(2.1325,2.0448,155.463,152.044),
  3:(2.1623,2.0321,155.125,150.608),  4:(2.1932,2.0200,154.718,149.239),
  5:(2.2256,2.0085,154.229,147.962)}
print("-- sec.7e referansi ile karsilastirma (tol: 0.002 A / 0.05 deg) --")
nbad = 0
for r in rows:
    k = int(round(r["eps"]))
    if k not in REF: continue
    rp, rq, tp, tq = REF[k]
    d = [r["Ni"][1]-rp, r["Ni"][0]-rq, r["th_par"]-tp, r["th_perp"]-tq]
    if max(abs(d[0]),abs(d[1])) > 2e-3 or max(abs(d[2]),abs(d[3])) > 5e-2:
        nbad += 1
        print("  SAPMA eps=%+d  dd_par=%+.4f dd_perp=%+.4f  dth_par=%+.3f dth_perp=%+.3f" % (k,*d))
print("  sapan nokta: %d / %d %s" % (nbad, len(rows), "-> sec.7e DOGRULANDI" if nbad==0 else "-> BAK"))
