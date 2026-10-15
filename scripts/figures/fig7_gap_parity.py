# ============================================================
#  Fig. 7 - parity decomposition of the spin-resolved band gaps (180 mm wide)
#  Input : data/gap.dat   (columns used: 0 = strain, 3 = majority gap, 4 = minority gap; eV)
#  Output: fig_gap_parity.pdf / .svg / .png
#  Also reproduces Table 4 (even/odd components) and the peak positions of Sec. 3.5,
#  and checks them against the values quoted in the manuscript.
#  Usage : python scripts/figures/fig7_gap_parity.py   (run from the repository root)
# ============================================================
import os, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

PATH = os.environ.get("GAP_DAT", os.path.join("data", "gap.dat"))
OUT = os.environ.get("FIG_OUT", ".")
os.makedirs(OUT, exist_ok=True)

MM = 1/25.4
OK = {"blue": "#0072B2", "orange": "#D55E00", "grey": "#666666"}
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "path",
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 8, "axes.labelsize": 8.5, "legend.fontsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.8,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "legend.frameon": False, "legend.handlelength": 1.6,
    "figure.dpi": 200,
})

# ---------------- data ----------------
d = np.loadtxt(PATH, usecols=(0, 3, 4))
o = np.argsort(d[:, 0]); e, up, dn = d[o, 0], d[o, 1], d[o, 2]
if np.abs(e).max() <= 1.0: e = e*100.0
i0 = int(np.argmin(np.abs(e)))
assert len(e) == 11 and abs(e[i0]) < 1e-9, "11 strain points including eps=0 required, got %d" % len(e)

def parity(v):
    k = np.arange(1, 6)
    C = np.array([(v[i0+j] + v[i0-j]) / 2 for j in k])
    T = np.array([(v[i0+j] - v[i0-j]) / 2 for j in k])
    return k, C, T
k, Cu, Tu = parity(up); _, Cd, Td = parity(dn)

def peak(x, y, deg=4):
    c = np.polyfit(x, y, deg); xs = np.linspace(x.min(), x.max(), 2001)
    return xs[np.argmax(np.polyval(c, xs))], c

print("Table 4 (meV):  |eps|  dC_up   T_up  |T/dC|_up   dC_dn   T_dn  |T/dC|_dn")
for j in range(5):
    dcu, dcd = (Cu[j]-up[i0])*1e3, (Cd[j]-dn[i0])*1e3
    print("               %d %8.1f %6.1f %8.3f %8.1f %6.1f %8.2f" % (k[j], dcu, Tu[j]*1e3, abs(Tu[j]*1e3/dcu),
                                                               dcd, Td[j]*1e3, abs(Td[j]*1e3/dcd)))

fig, ax = plt.subplots(1, 2, figsize=(180*MM, 68*MM))

# ---------------- (a) E_g(eps) ----------------
a = ax[0]
pu, cu = peak(e, up); pd, cd = peak(e, dn)
xs = np.linspace(-5, 5, 400)
a.plot(xs, np.polyval(cu, xs), "-", color=OK["blue"], lw=1.0, zorder=1)
a.plot(e, up, "o", color=OK["blue"], ms=4, zorder=3, label="majority (left)")
a.axvline(pu, color=OK["blue"], lw=0.7, ls=":", zorder=0)
a.set_xlabel(r"$\varepsilon$  (%)")
a.set_ylabel(r"$E_g^{\uparrow}$  (eV)", color=OK["blue"])
a.tick_params(axis="y", colors=OK["blue"]); a.set_xlim(-5.6, 5.6)
a.set_ylim(0.985, 1.315)
a.xaxis.set_major_locator(MultipleLocator(2))
a2 = a.twinx()
a2.plot(xs, np.polyval(cd, xs), "--", color=OK["orange"], lw=1.0, zorder=1)
a2.plot(e, dn, "s", color=OK["orange"], ms=4, mfc="white", zorder=3, label="minority (right)")
a2.axvline(pd, color=OK["orange"], lw=0.7, ls=":", zorder=0)
a2.set_ylabel(r"$E_g^{\downarrow}$  (eV)", color=OK["orange"])
a2.tick_params(axis="y", colors=OK["orange"], direction="in"); a2.minorticks_on()
a2.set_ylim(2.79, 3.265)
fx = lambda x: (x - a.get_xlim()[0]) / (a.get_xlim()[1] - a.get_xlim()[0])
a.text(fx(pu)+0.012, 0.30, r"peak $%+.2f$ %%" % pu, transform=a.transAxes,
       color=OK["blue"], fontsize=7, ha="left")
a.text(fx(pd)-0.012, 0.30, r"peak $%+.2f$ %%" % pd, transform=a.transAxes,
       color=OK["orange"], fontsize=7, ha="right")
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
a.legend(h1+h2, l1+l2, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2,
         columnspacing=1.2, handletextpad=0.4)
a.text(0.02, 0.90, "(a)", transform=a.transAxes, fontweight="bold")

# ---------------- (b) odd component ----------------
b = ax[1]
b.axhline(0, color=OK["grey"], lw=0.7, zorder=0)
b.plot(np.r_[0, k], np.r_[0, Tu*1000], "o-",  color=OK["blue"],   ms=4, lw=1.2, label="majority")
b.plot(np.r_[0, k], np.r_[0, Td*1000], "s--", color=OK["orange"], ms=4, lw=1.2,
       mfc="white", label="minority")
b.set_xlabel(r"$|\varepsilon|$  (%)"); b.set_ylabel(r"odd component $T$  (meV)")
b.set_xlim(-0.25, 5.25); b.xaxis.set_major_locator(MultipleLocator(1))
b.legend(loc="lower left")
b.text(0.02, 0.90, "(b)", transform=b.transAxes, fontweight="bold")
ju, jd = int(np.argmax(np.abs(Tu))), int(np.argmax(np.abs(Td)))
au, ad = abs(Tu[ju])*1000, abs(Td[jd])*1000
b.annotate("max |T| = %.1f meV\nat $|\\varepsilon|$ = %d %%\n%.2f %% of $E_g(0)$" % (au+1e-6, k[ju], 100*au/1000/up[i0]),
           xy=(k[ju], Tu[ju]*1000), xytext=(1.55, -30), color=OK["blue"], fontsize=7, va="center",
           arrowprops=dict(arrowstyle="-", color=OK["blue"], lw=0.7, shrinkB=3))
b.annotate("max |T| = %.1f meV\nat $|\\varepsilon|$ = %d %%\n%.2f %% of $E_g(0)$" % (ad+1e-6, k[jd], 100*ad/1000/dn[i0]),
           xy=(k[jd], Td[jd]*1000), xytext=(0.30, -93), color=OK["orange"], fontsize=7, va="center",
           arrowprops=dict(arrowstyle="-", color=OK["orange"], lw=0.7, shrinkB=3))

fig.tight_layout(pad=0.6)
for ext, kw in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
    p = os.path.join(OUT, "fig_gap_parity." + ext)
    fig.savefig(p, **kw); print("written:", p)

# ---------------- check against the manuscript ----------------
REF = [("peak, majority (%)",   pu,                 +0.04,  0.05),
       ("peak, minority (%)",   pd,                 -1.58,  0.05),
       ("dC up at |eps|=5 (eV)", Cu[-1]-up[i0],     -0.2150, 0.001),
       ("dC dn at |eps|=5 (eV)", Cd[-1]-dn[i0],     -0.2047, 0.001),
       ("max|T| up (eV)",        au/1000,            0.0095, 0.001),
       ("max|T| up at |eps|",    k[ju],              1,      0),
       ("max|T| dn (eV)",        ad/1000,            0.1306, 0.001),
       ("odd/E_g(0) up (%)",     100*au/1000/up[i0], 0.73,   0.02),
       ("odd/E_g(0) dn (%)",     100*ad/1000/dn[i0], 4.04,   0.02)]
bad = 0
for nm, got, want, tol in REF:
    ok = abs(got - want) <= tol + 1e-12; bad += (not ok)
    print("  %-24s computed=%9.4f  manuscript=%8.4f  %s" % (nm, got, want, "OK" if ok else "<-- MISMATCH"))
print("check: %d/%d mismatches" % (bad, len(REF)))
