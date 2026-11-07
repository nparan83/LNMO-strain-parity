# tilt_gen.py - Glazer tilt deseninden baslangic yapisi; uretim hucresi ve B konumlari korunur.
# Kullanim: python3 tilt_gen.py scf_FM.in aaa|aac alfa gamma(derece) > pos.txt   (aac 9 8 -> uretim yapisina ~0.08 A)
import sys, numpy as np
src, pat, al, ga = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
L = open(src).read().splitlines()
i = [k for k, l in enumerate(L) if l.strip().upper().startswith('CELL_PARAMETERS')][0]
A = np.array([[float(x) for x in L[i + j].split()[:3]] for j in (1, 2, 3)])
j = [k for k, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')][0]
S = [l.split()[0] for l in L[j + 1:j + 21]]; X = np.array([[float(v) for v in l.split()[1:4]] for l in L[j + 1:j + 21]])
Ef = np.array([[.5, .5, 0], [-.5, .5, 0], [.5, 0, .5]]); E = Ef @ A; Ax = E / np.linalg.norm(E, axis=1)[:, None]  # psodokubik
Bi = [k for k, s in enumerate(S) if s[:2] in ('Ni', 'Mn')]; o = X[Bi[0]]
def R(w):
    v = Ax.T @ w; t = np.linalg.norm(v); k = v / t; K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K
def om(f):
    n = np.rint((f - o) @ np.linalg.inv(Ef)).astype(int); s = 1 - 2 * (n.sum() % 2)
    return np.radians([al * s, al * s, ga * (s if pat == 'aaa' else 1 - 2 * ((n[0] + n[1]) % 2))])
wrap = lambda f: np.mod(f + 1e-7, 1) - 1e-7
out, mis = [], []
for b in Bi:
    for k in range(3):                       # her B'nin +e_k komsusuyla paylastigi O
        f1 = X[b]; f2 = X[b] + Ef[k]; d = E[k]
        p1 = f1 @ A + R(om(f1)) @ (d / 2); p2 = f2 @ A - R(om(f2)) @ (d / 2)
        out.append(('O', wrap(0.5 * (p1 + p2) @ np.linalg.inv(A)))); mis.append(np.linalg.norm(p1 - p2))
la = []
for b in Bi:
    for sg in np.ndindex(2, 2, 2):
        f = wrap(X[b] + 0.5 * (2 * np.array(sg) - 1) @ Ef)
        if not any(np.linalg.norm((f - g) - np.rint(f - g)) < 1e-4 for g in la): la.append(f)
assert len(out) == 12 and len(la) == 4, (len(out), len(la))
print('# kose uyumsuzlugu max %.3f A (rijit oktahedron yaklasimi; relax giderir)' % max(mis), file=sys.stderr)
for f in la: print('La  %14.10f %14.10f %14.10f' % tuple(f))
for k in Bi: print('%-3s %14.10f %14.10f %14.10f' % (S[k], *X[k]))
for s, f in out: print('O   %14.10f %14.10f %14.10f' % tuple(f))
