#!/usr/bin/env python3
# gap.py -- paper-1'in E_g(|Q3|) verisini kaynak .out dosyalarindan yeniden cikarir,
# lineer vs kuadratik model karsilastirmasi yapar ve artiklari basar.
import os, sys, glob, re, math

ROOT = sys.argv[1] if len(sys.argv) > 1 else "LNMO_Strain_Production"
Q3F  = sys.argv[2] if len(sys.argv) > 2 else "q3_geo.dat"

KRE = re.compile(r"^\s+k =.*bands \(ev\):\s*$")
NUM = re.compile(r"-?\d+\.\d+")

def parse(fn):
    try:
        L = open(fn, errors="ignore").read().splitlines()
    except Exception:
        return None
    nelec = mag = None
    for l in L:
        if "number of electrons" in l:
            m = NUM.search(l)
            if m: nelec = float(m.group())
        if "total magnetization" in l and "Bohr" in l:
            m = NUM.search(l.split("=")[-1])
            if m: mag = float(m.group())
    if nelec is None: return None
    if mag is None: mag = 0.0
    direct = None
    for l in L:
        if "highest occupied" in l and "lowest unoccupied" in l:
            v = [float(x) for x in NUM.findall(l)]
            if len(v) >= 2: direct = (v[-2], v[-1])
    chan, cur = {"up": [], "dn": []}, "up"
    i, n = 0, len(L)
    spinpol = any("SPIN UP" in l for l in L)
    while i < n:
        l = L[i]
        if "SPIN UP" in l: cur = "up"
        elif "SPIN DOWN" in l: cur = "dn"
        if KRE.match(l):
            i += 1
            ev = []
            while i < n:
                s = L[i]
                if ("occupation numbers" in s or KRE.match(s) or "the Fermi" in s
                        or "highest occupied" in s or "SPIN" in s or "Writing" in s):
                    break
                f = NUM.findall(s)
                if f: ev += [float(x) for x in f]
                elif ev: break
                i += 1
            if ev: chan[cur].append(sorted(ev))
            continue
        i += 1
    if not spinpol:
        chan["dn"] = chan["up"]
    if not chan["up"] or not chan["dn"]:
        return {"f": fn, "nelec": nelec, "mag": mag, "direct": direct, "gap": None}
    nup = (nelec + mag) / 2.0
    ndn = (nelec - mag) / 2.0
    out = {"f": fn, "nelec": nelec, "mag": mag, "direct": direct,
           "nup": nup, "ndn": ndn, "nk": len(chan["up"])}
    for key, nocc in (("up", nup), ("dn", ndn)):
        k = int(round(nocc))
        if abs(nocc - k) > 1e-6: out[key] = None; continue
        B = chan[key]
        if any(len(b) < k + 1 for b in B): out[key] = None; continue
        vbm = max(b[k - 1] for b in B)
        cbm = min(b[k] for b in B)
        out[key] = (vbm, cbm, cbm - vbm)
    g = [out[k][2] for k in ("up", "dn") if out.get(k)]
    out["gap"] = min(g) if g else None
    return out

def eps_of(d):
    m = re.search(r"strain_(-?\d+\.?\d*)", d)
    if not m: return None
    e = float(m.group(1))
    return e * 100.0 if abs(e) <= 1.0 else e

dirs = sorted(glob.glob(os.path.join(ROOT, "strain_*", "FM")))
if not dirs:
    dirs = sorted(glob.glob(os.path.join(ROOT, "strain_*")))
print("bulunan dizin sayisi:", len(dirs))

rows = []
for d in dirs:
    e = eps_of(d)
    outs = sorted(glob.glob(os.path.join(d, "*.out")))
    cands = []
    for f in outs:
        r = parse(f)
        if r and (r.get("gap") is not None or r.get("direct")):
            cands.append(r)
    if not cands:
        print("  %-40s -> bosluk cikarilamadi (%d .out)" % (d, len(outs)))
        continue
    def rank(r):
        b = os.path.basename(r["f"])
        return (0 if b.startswith("nscf") else 1 if b.startswith("scf") else 2, -r.get("nk", 0), b)
    cands.sort(key=rank)
    best = cands[0]
    gaps = sorted(set(round(c["gap"], 4) for c in cands if c.get("gap") is not None))
    flag = " **" if len(gaps) > 1 and (max(gaps) - min(gaps)) > 0.01 else ""
    rows.append((e, best, cands, flag))
rows.sort(key=lambda r: r[0])

print("\n%-6s %-10s %-10s %-10s %-8s %-6s %s" % ("eps", "gap_up", "gap_dn", "gap_min", "Mtot", "nk", "dosya"))
for e, b, c, fl in rows:
    gu = "%.4f" % b["up"][2] if b.get("up") else "nan"
    gd = "%.4f" % b["dn"][2] if b.get("dn") else "nan"
    gm = "%.4f" % b["gap"] if b.get("gap") is not None else "nan"
    print("%+6.1f %-10s %-10s %-10s %-8.2f %-6d %s%s"
          % (e, gu, gd, gm, b["mag"], b.get("nk", 0), os.path.basename(b["f"]), fl))
    if len(c) > 1:
        alt = ", ".join("%s=%s" % (os.path.basename(x["f"]),
                                   ("%.4f" % x["gap"]) if x.get("gap") is not None else "nan")
                        for x in c[1:])
        print("        adaylar: %s" % alt)

Q3 = {}
if os.path.exists(Q3F):
    for l in open(Q3F):
        if l.startswith("#") or not l.strip(): continue
        t = l.split()
        try:
            e = float(t[0]); q = float(t[3])
        except Exception:
            continue
        if abs(e) <= 1.0: e *= 100.0
        Q3[round(e)] = q
else:
    print("\nUYARI: %s yok, Q3 analizi atlaniyor." % Q3F)

def polyfit(x, y, deg):
    n = deg + 1
    A = [[sum(xi ** (i + j) for xi in x) for j in range(n)] for i in range(n)]
    b = [sum(yi * xi ** i for xi, yi in zip(x, y)) for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(A[r][i]))
        A[i], A[p] = A[p], A[i]; b[i], b[p] = b[p], b[i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            for cc in range(i, n): A[r][cc] -= f * A[i][cc]
            b[r] -= f * b[i]
    c = [0.0] * n
    for i in range(n - 1, -1, -1):
        c[i] = (b[i] - sum(A[i][j] * c[j] for j in range(i + 1, n))) / A[i][i]
    return c

def ev(c, x): return sum(ci * x ** i for i, ci in enumerate(c))
def rmse(x, y, c): return math.sqrt(sum((yi - ev(c, xi)) ** 2 for xi, yi in zip(x, y)) / len(x))

data = [(e, b["gap"]) for e, b, c, fl in rows if b.get("gap") is not None and round(e) in Q3]
if len(data) >= 5:
    q = [abs(Q3[round(e)]) for e, g in data]
    y = [g for e, g in data]
    xs = [e for e, g in data]
    print("\n=== E_g vs |Q3| ===")
    for d in (1, 2, 3):
        c = polyfit(q, y, d)
        print("derece %d : RMSE %.5f eV | katsayilar %s" % (d, rmse(q, y, c),
              " ".join("%+.4f" % ci for ci in c)))
    c1 = polyfit(q, y, 1); c2 = polyfit(q, y, 2)
    print("\nlineer model: E_g = %+.4f %+.4f*|Q3|   (paper-1: 1.327 - 1.371|Q3|)" % (c1[0], c1[1]))
    print("\n%-7s %-9s %-9s %-10s %-10s" % ("eps", "|Q3|", "E_g", "artik(d1)", "artik(d2)"))
    r1 = []
    for xi, qi, yi in zip(xs, q, y):
        a1 = yi - ev(c1, qi); a2 = yi - ev(c2, qi); r1.append(a1)
        print("%+7.1f %-9.5f %-9.4f %+10.4f %+10.4f" % (xi, qi, yi, a1, a2))
    npos = sum(1 for a in r1 if a > 0)
    print("\nlineer artiklarin isaret dagilimi: %d pozitif / %d negatif" % (npos, len(r1) - npos))
    print("lineer max|artik| = %.4f eV ; kuadratik max|artik| = %.4f eV"
          % (max(abs(a) for a in r1), max(abs(yi - ev(c2, qi)) for qi, yi in zip(q, y))))
    print("\n=== DAL ASIMETRISI (basma vs cekme, benzer |Q3|) ===")
    comp = sorted([(qi, yi, xi) for xi, qi, yi in zip(xs, q, y) if xi < -0.5])
    tens = sorted([(qi, yi, xi) for xi, qi, yi in zip(xs, q, y) if xi > 0.5])
    dif = []
    for qc, yc, xc in comp:
        if not tens: break
        qt, yt, xt = min(tens, key=lambda t: abs(t[0] - qc))
        d = yc - yt
        dif.append(d)
        print("  eps %+.0f (|Q3|=%.5f) vs eps %+.0f (|Q3|=%.5f) : dE_g = %+.4f eV"
              % (xc, qc, xt, qt, d))
    if dif:
        mn = sum(dif) / len(dif)
        print("  ortalama dal farki = %+.4f eV ; hepsi ayni isaretli mi: %s"
              % (mn, "EVET" if all(d > 0 for d in dif) or all(d < 0 for d in dif) else "hayir"))
        print("  -> |%.4f| eV sistematikse E_g |Q3|'un TEK DEGERLI fonksiyonu degildir;" % mn)
        print("     paper-1'in cift parite iddiasi 'bu mertebeye kadar' diye sinirlandirilmali.")
    imp = rmse(q, y, c1) / rmse(q, y, c2) if rmse(q, y, c2) > 0 else float('inf')
    print("RMSE iyilesmesi (d1/d2) = %.2f kat" % imp)
    print("\nKARAR: %s" % ("KUADRATIK GEREKLI - paper-1'in denklemi guncellenmeli"
          if imp > 1.8 else "LINEER YETERLI - paper-1'in denklemi ayakta"))
    print("\n=== paper-1 referans kontrolu ===")
    def near(e):
        for xi, yi in zip(xs, y):
            if abs(xi - e) < 0.5: return yi
        return None
    for e, want in ((0, 1.295), (-5, 1.076), (+5, 1.084)):
        got = near(e)
        ok = got is not None and abs(got - want) < 0.02
        print("  eps=%+3d  hesaplanan=%s  paper-1=%.3f  %s"
              % (e, ("%.4f" % got) if got is not None else "nan", want,
                 "OK" if ok else "<-- SAPMA"))
    with open("gap.dat", "w") as f:
        f.write("# eps |Q3| E_g_min E_g_up E_g_dn dosya\n")
        for (e, b, c, fl) in rows:
            if b.get("gap") is None or round(e) not in Q3: continue
            gu = b["up"][2] if b.get("up") else float('nan')
            gd = b["dn"][2] if b.get("dn") else float('nan')
            f.write("%6.2f %10.5f %10.4f %10.4f %10.4f %s\n"
                    % (e, abs(Q3[round(e)]), b["gap"], gu, gd, b["f"]))
    print("\nyazildi: gap.dat (%d satir)" % len(data))
else:
    print("\nYETERSIZ VERI: %d nokta." % len(data))
    print("Oneri: ozdegerler basilmamis olabilir; .dos dosyalarindan cikarim gerekir.")
