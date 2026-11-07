# mkin.py - uretim sablonundan girdi. relax: mkin.py relax scf_FM.in pos.txt | scf: mkin.py scf scf_FM.in vcrelax.out
# restart: mkin.py restart vcrelax.in vcrelax.out (SON adimin hucre+konumu, nstep=200)
import sys, re
mode, tpl, src = sys.argv[1:4]; T = open(tpl).read()
def setpos(T, cell, pos):
    out, skip = [], 0
    for l in T.splitlines():
        u = l.strip().upper()
        if skip: skip -= 1; continue
        if u.startswith('CELL_PARAMETERS') and cell: out += ['CELL_PARAMETERS (angstrom)'] + cell; skip = 3; continue
        if u.startswith('ATOMIC_POSITIONS'): out += ['ATOMIC_POSITIONS (crystal)'] + pos; skip = 20; continue
        out.append(l)
    return '\n'.join(out) + '\n'
if mode == 'relax':
    T = setpos(T, None, [l for l in open(src).read().splitlines() if l.strip()][:20])
    T = T.replace("'scf'", "'vc-relax'", 1).replace('1.0d-10', '1.0d-8', 1).replace("'vc-relax'", "'vc-relax'\n    nstep         = 200", 1)
    T = T.replace('ATOMIC_SPECIES', "&IONS\n    ion_dynamics    = 'bfgs'\n/\n&CELL\n    cell_dynamics   = 'bfgs'\n    press           = 0.0\n"
                  "    press_conv_thr  = 0.1\n    cell_dofree     = 'epitaxial_ab'\n/\nATOMIC_SPECIES", 1)
else:
    t = open(src).read(); a, b = t.rfind('Begin final coordinates'), t.rfind('End final coordinates')
    if mode == 'restart': a, b = t.rfind('CELL_PARAMETERS'), len(t)
    elif a < 0 or b < a: sys.exit('vc-relax bitmemis')
    L = t[a:b].splitlines(); i = [k for k, l in enumerate(L) if l.strip().startswith('CELL_PARAMETERS')][0]
    m = re.search(r'alat=\s*([\d.]+)', L[i]); f = float(m.group(1)) * 0.529177210903 if m else (0.529177210903 if 'bohr' in L[i] else 1.0)
    cell = ['%16.10f%16.10f%16.10f' % tuple(float(x) * f for x in L[i + k].split()) for k in (1, 2, 3)]
    j = [k for k, l in enumerate(L) if l.strip().startswith('ATOMIC_POSITIONS')][0]; assert 'crystal' in L[j]
    T = setpos(T, cell, L[j + 1:j + 21])
if mode == 'restart' and 'nstep' not in T: T = T.replace("'vc-relax'", "'vc-relax'\n    nstep         = 200", 1)
print(T, end='')
