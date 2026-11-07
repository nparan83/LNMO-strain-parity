# el_post.py - C_ij (GPa), Born kararliligi, epitaksiyel tahmin d ln c_perp / d eps. Kullanim: python3 el_post.py WORKDIR
import sys, re, glob, numpy as np
W = sys.argv[1]; d = 0.005; P = '/arf/scratch/nparan/lnmo/LNMO_Strain_Production'
def stress(f):
    t = open(f).read(); k = t.rfind('total   stress'); ok = 'JOB DONE' in t and k > 0
    S = np.array([[float(x) for x in l.split()[3:6]] for l in t[k:].splitlines()[1:4]]) / 10 if ok else None   # kbar -> GPa
    return S, ok
C = np.zeros((6, 6)); bad = []
for D in glob.glob(f'{W}/*/'):
    j, s = map(int, open(D + 'tag').read().split()); S, ok = stress(D + 'pw.out')
    if not ok: bad.append(D); continue
    v = np.array([S[0, 0], S[1, 1], S[2, 2], S[1, 2], S[0, 2], S[0, 1]])
    C[:, j] += -s * v / (2 * d)                                               # QE gerilmesi: pozitif = basinc -> isaret degisir
if bad: sys.exit('bitmemis: %s' % bad)
print('asimetri max %.2f GPa' % np.abs(C - C.T).max()); C = 0.5 * (C + C.T)
np.savetxt(f'{W}/Cij.dat', C, fmt='%9.2f', header='C_ij (GPa), Voigt xx yy zz yz xz xy; kartezyen: x||a, y||b'); print(np.round(C, 1))
ev = np.linalg.eigvalsh(C); print('ozdegerler (GPa):', np.round(ev, 1), '-> Born', 'KARARLI' if ev.min() > 0 else 'KARARSIZ')
f, fx = [2, 3, 4], [0, 1]                                                     # epitaksiyel: e_xx = e_yy = eps, e_xy = 0; sigma_zz,yz,xz = 0
e = np.linalg.solve(C[np.ix_(f, f)], -C[np.ix_(f, fx)].sum(axis=1)); print('tahmin d e_zz/d eps = %.3f' % e[0])
def cp(fn):
    L = open(fn).read().splitlines(); i = [k for k, l in enumerate(L) if l.strip().upper().startswith('CELL_PARAMETERS')][0]
    A = np.array([[float(x) for x in L[i + j].split()[:3]] for j in (1, 2, 3)]); return np.linalg.det(A) / np.linalg.norm(np.cross(A[0], A[1]))
try: print('DFT seri d ln c_perp / d eps (+-1 %%) = %.3f' % ((np.log(cp(f'{P}/strain_0.01/FM/scf_FM.in')) - np.log(cp(f'{P}/strain_-0.01/FM/scf_FM.in'))) / 0.02))
except Exception as x: print('c_perp kontrolu yapilamadi:', x)
