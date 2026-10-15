#!/usr/bin/env python3
# gap3.py -- E_g(eps) her iki spin kanali icin; gap.dat + parite ayristirmasi.
# DUZELTME: 'total magnetization' yoksa (nscf ciktilarinda yok) 0 VARSAYILMAZ;
# dolu bant sayisi E_F'yi iceren araligi arayarak bulunur.
import os, sys, glob, re, math

ROOT = sys.argv[1] if len(sys.argv) > 1 else "LNMO_Strain_Production"
Q3F  = sys.argv[2] if len(sys.argv) > 2 else "q3_geo.dat"
KRE = re.compile(r"^\s+k =.*bands \(ev\):\s*$")
NUM = re.compile(r"-?\d+\.\d+")

def read(fn):
    L = open(fn, errors="ignore").read().splitlines()
    nelec = mag = ef = None
    for l in L:
        if "number of electrons" in l:
            m = NUM.search(l); nelec = float(m.group()) if m else nelec
        if "total magnetization" in l and "Bohr" in l:
            m = NUM.search(l.split("=")[-1]); mag = float(m.group()) if m else mag
        if "the Fermi energy is" in l:
            m = NUM.search(l); ef = float(m.group()) if m else ef
    chan, cur = {"up": [], "dn": []}, "up"
    i, n = 0, len(L)
    spinpol = any("SPIN UP" in x for x in L)
    while i < n:
        l = L[i]
        if "SPIN UP" in l: cur = "up"
        elif "SPIN DOWN" in l: cur = "dn"
        if KRE.match(l):
            i += 1; ev = []
            while i < n:
                s = L[i]
                if ("occupation numbers" in s or KRE.match(s) or "the Fermi" in s
                        or "highest occupied" in s or "SPIN" in s or "Writing" in s
                        or "End of" in s): break
                f = NUM.findall(s)
                if f: ev += [float(x) for x in f]
                elif ev: break
                i += 1
            if ev: chan[cur].append(sorted(ev))
            continue
        i += 1
    if not spinpol: chan["dn"] = list(chan["up"])
    return dict(f=fn, nelec=nelec, mag=mag, ef=ef, chan=chan)

def pick(B, ef):
    if not B or ef is None: return None
    nb = min(len(b) for b in B)
    best = None
    for nn in range(1, nb):
        v = max(b[nn - 1] for b in B); c = min(b[nn] for b in B)
        if v <= ef <= c and (best is None or c - v > best[3]):
            best = (nn, v, c, c - v)
    return best

def eps_of(d):
    m = re.search(r"strain_(-?\d+\.?\d*)", d)
    e = float(m.group(1)); return e * 100 if abs(e) <= 1 else e

Q3 = {}
if os.path.exists(Q3F):
    for l in open(Q3F):
        if l.startswith("#") or not l.strip(): continue
        t = l.split()
        try: e = float(t[0]); q = float(t[3])
        except Exception: continue
        if abs(e) <= 1.0: e *= 100.0
        Q3[round(e)] = q

rows = []
for d in sorted(glob.glob(os.path.join(ROOT, "strain_*", "FM"))):
    e = eps_of(d)
    got = None
    for f in sorted(glob.glob(os.path.join(d, "*.out"))):
        if not os.path.basename(f).startswith(("nscf_FM", "scf_FM")): continue
        r = read(f)
        if not r["chan"]["up"]: continue
        if r["ef"] is None:
            print("  UYARI: %s -> E_F yok, atlaniyor" % f); continue
        u = pick(r["chan"]["up"], r["ef"]); dn = pick(r["chan"]["dn"], r["ef"])
        if not u or not dn:
            print("  UYARI: %s -> E_F hicbir araligin icinde degil (METALIK?)" % f); continue
        got = (f, r, u, dn)
        if os.path.basename(f).startswith("nscf"): break
    if got is None:
        print("  %+.0f%% : cikarim BASARISIZ" % e); continue
    f, r, u, dn = got
    chk = ""
    if r["nelec"] is not None and r["mag"] is not None:
        nu = (r["nelec"] + r["mag"]) / 2; nd = (r["nelec"] - r["mag"]) / 2
        chk = "OK" if (abs(nu - u[0]) < 0.5 and abs(nd - dn[0]) < 0.5) else "<-- SAYIM UYUSMUYOR"
    rows.append((e, u, dn, r, f, chk))
rows.sort(key=lambda x: x[0])

print("\n%-7s %-5s %-9s %-5s %-9s %-8s %-8s %s" %
      ("eps", "n_up", "gap_up", "n_dn", "gap_dn", "E_F", "kontrol", "dosya"))
for e, u, dn, r, f, chk in rows:
    print("%+7.1f %-5d %-9.4f %-5d %-9.4f %-8.3f %-8s %s"
          % (e, u[0], u[3], dn[0], dn[3], r["ef"], chk, os.path.basename(f)))

with open("gap.dat", "w") as fh:
    fh.write("# eps Q3_Ni |Q3| Eg_up Eg_dn n_up n_dn E_F dosya\n")
    for e, u, dn, r, f, chk in rows:
        q = Q3.get(round(e), float('nan'))
        fh.write("%6.2f %10.5f %10.5f %10.4f %10.4f %5d %5d %10.4f %s\n"
                 % (e, q, abs(q), u[3], dn[3], u[0], dn[0], r["ef"], f))
print("\nyazildi: gap.dat (%d satir)" % len(rows))

E = [x[0] for x in rows]
if len(rows) == 11 and abs(E[5]) < 1e-6:
    for nm, idx in (("E_g cogunluk (up)", 1), ("E_g azinlik (dn)", 2)):
        v = [x[idx][3] for x in rows]
        print("\n=== %s : parite ayristirmasi f(e)=C(e)+T(e) ===" % nm)
        print(" |eps|      cift C        tek T")
        C, T = [], []
        for k in range(1, 6):
            a, b = v[5 + k], v[5 - k]
            C.append((a + b) / 2); T.append((a - b) / 2)
            print("   %d      %9.4f    %+9.4f" % (k, C[-1], T[-1]))
        ac = max(abs(c - v[5]) for c in C); at = max(map(abs, T))
        print("  f(0)=%.4f | cift genlik %.4f eV | TEK genlik %.4f eV | tek/cift %.3f"
              % (v[5], ac, at, at / ac if ac else float('nan')))
        print("  tek bilesen / f(0) = %%%.2f" % (100 * at / v[5]))

    print("\n=== paper-1 referans kontrolu (cogunluk) ===")
    for e0, want in ((0, 1.295), (-5, 1.076), (+5, 1.084)):
        g = [x[1][3] for x in rows if abs(x[0] - e0) < 0.5]
        ok = g and abs(g[0] - want) < 0.005
        print("  eps=%+3d hesap=%s paper-1=%.3f %s"
              % (e0, ("%.4f" % g[0]) if g else "-", want, "OK" if ok else "<-- SAPMA"))
else:
    print("\n(11 nokta + eps=0 yoksa parite ayristirmasi yapilmaz: %d nokta)" % len(rows))
