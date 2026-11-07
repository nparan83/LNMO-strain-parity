# tilt_post.py - a-a-a- vs uretim (a-a-c+): enerji farki, moment, c_perp, z-donme deseni, uzay grubu
import re, glob, numpy as np, spglib
P = '/arf/scratch/nparan/lnmo/LNMO_Strain_Production'; R = '/arf/scratch/nparan/R1/tilt'; RY = 13605.693
num = {'La': 57, 'Mn': 25, 'Ni': 28, 'O': 8}
def struct(f):
    L = open(f).read().splitlines(); i = [k for k, l in enumerate(L) if l.strip().upper().startswith('CELL_PARAMETERS')][0]
    A = np.array([[float(x) for x in L[i + j].split()[:3]] for j in (1, 2, 3)]); j = [k for k, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')][0]
    return A, [l.split()[0] for l in L[j + 1:j + 21]], np.array([[float(v) for v in l.split()[1:4]] for l in L[j + 1:j + 21]])
def rotz(A, S, X):                          # B oktahedronlarinin psodokubik z etrafinda donmesi (derece)
    E = np.array([[.5, .5, 0], [-.5, .5, 0], [.5, 0, .5]]) @ A; Ax = E / np.linalg.norm(E, axis=1)[:, None]; C = X @ A; w = {}
    for s in ('Mn1', 'Ni1'):
        c = C[S.index(s)]; v = np.array([C[k] + (np.array(n) - 1) @ A - c for k, t in enumerate(S) if t == 'O' for n in np.ndindex(3, 3, 3)])
        v = v[np.argsort(np.linalg.norm(v, axis=1))][:6]; u = v / np.linalg.norm(v, axis=1)[:, None]
        M = np.array([u[np.argmax(u @ ax)] for ax in Ax]) @ Ax.T; w[s] = np.degrees(np.arctan2(M[0, 1], M[0, 0]))
    return w
E = lambda f: float(re.findall(r'!\s+total energy\s+=\s+(\S+)', open(f).read())[-1])
mag = lambda f: re.findall(r'total magnetization\s+=\s+(\S+)', open(f).read())[-1]
print('eps   dE(aaa-aac) meV/f.u.  M   c_perp(aaa/aac) A   wz Mn1/Ni1 aaa | aac   grup')
for d in sorted(glob.glob(f'{R}/eps_*')):
    e = d.split('_')[-1]; p = f'{P}/strain_{e}/FM'
    try: dE = (E(f'{d}/scf.out') - E(f'{p}/scf_FM.out')) * RY / 2
    except Exception: print(e, 'bitmemis'); continue
    A, S, X = struct(f'{d}/scf.in'); A0, S0, X0 = struct(f'{p}/scf_FM.in')
    cp = lambda A: np.linalg.det(A) / np.linalg.norm(np.cross(A[0], A[1]))
    w, w0 = rotz(A, S, X), rotz(A0, S0, X0); g = spglib.get_spacegroup((A, X, [num[s[:2].rstrip('0123456789')] for s in S]), symprec=1e-3)
    print('%s %8.1f  %s  %.4f/%.4f  %+.1f/%+.1f | %+.1f/%+.1f  %s' % (e, dE, mag(f'{d}/scf.out'), cp(A), cp(A0), w['Mn1'], w['Ni1'], w0['Mn1'], w0['Ni1'], g))
