# ph_irr.py - Gamma modlarinin indirgenemez temsilleri (Ag/Bg Raman, Au/Bu IR). Kullanim: python3 ph_irr.py WORKDIR
import sys, re, numpy as np, phonopy, io, contextlib
W = sys.argv[1]; RYB = 25.71104309541616
def F(f):
    t = open(f).read(); k = t.rfind('Forces acting on atoms')
    return np.array([list(map(float, l.split('=')[1].split()[:3])) for l in t[k:].splitlines() if 'force =' in l and 'atom' in l][:20]) * RYB
n = int(open(f'{W}/ndisp').read()); R = [F(f'{W}/d{i}/pw.out') for i in range(n + 1)]
ph = phonopy.load(f'{W}/phonopy_disp.yaml', produce_fc=False, log_level=0)
ph.forces = np.array([r - R[0] for r in R[1:]]); ph.produce_force_constants(); ph.symmetrize_force_constants()
ph.set_irreps([0, 0, 0], degeneracy_tolerance=1e-3)
buf = io.StringIO()
with contextlib.redirect_stdout(buf): ph.irreps.show()
out = buf.getvalue(); open(f'{W}/gamma_irreps.txt', 'w').write(out)
for l in out.splitlines():
    m = re.match(r'\s*(\d+)\s+\(\s*(-?[\d.]+)\):\s+(\S+)', l)
    if m: print('%3s %9.2f cm-1  %s' % (m.group(1), float(m.group(2)) * 33.35641, m.group(3) if m.group(3)[0] in 'ABE' else '?'))
