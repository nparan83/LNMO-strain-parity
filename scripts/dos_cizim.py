import matplotlib.pyplot as plt
import numpy as np
import os

def read_dos(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    ef = float(lines[0].split()[-2])
    data = np.loadtxt(lines[1:])
    e = data[:, 0] - ef  # Fermi seviyesini sıfıra kaydır
    dos_up = data[:, 1]
    dos_down = data[:, 2]
    
    if np.mean(dos_down) > 0:
        dos_down = -dos_down
        
    return e, dos_up, dos_down

e0, u0_up, u0_down = read_dos("U0/lmno_U0.dos")
e3, u3_up, u3_down = read_dos("U3/lmno_U3.dos")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# ===== U=0 PANEL =====
ax1.plot(e0, u0_up, color='royalblue', linewidth=1.5)
ax1.plot(e0, u0_down, color='royalblue', linewidth=1.5)
ax1.fill_between(e0, 0, u0_up, color='royalblue', alpha=0.3)
ax1.fill_between(e0, 0, u0_down, color='royalblue', alpha=0.3)
ax1.axvline(0, color='black', linestyle='--', linewidth=1.5)
ax1.set_title("Standard DFT (U = 0 eV)", fontsize=14, fontweight='bold')
ax1.set_ylabel("Density of States", fontsize=12)
ax1.grid(True, linestyle=':', alpha=0.6)

# ===== U=3 PANEL =====
ax2.plot(e3, u3_up, color='crimson', linewidth=1.5)
ax2.plot(e3, u3_down, color='crimson', linewidth=1.5)
ax2.fill_between(e3, 0, u3_up, color='crimson', alpha=0.3)
ax2.fill_between(e3, 0, u3_down, color='crimson', alpha=0.3)
ax2.axvline(0, color='black', linestyle='--', linewidth=1.5)
ax2.set_title("DFT+U (U = 3 eV)", fontsize=14, fontweight='bold')
ax2.set_xlabel("Energy relative to Fermi level (eV)", fontsize=12)
ax2.set_ylabel("Density of States", fontsize=12)
ax2.grid(True, linestyle=':', alpha=0.6)

max_dos = max(np.max(u0_up), np.max(u3_up))
limit = min(max_dos + 5, 40)
ax1.set_ylim(-limit, limit)
ax2.set_ylim(-limit, limit)

plt.xlim(-7.2, 5.0)
plt.tight_layout()
plt.savefig("DOS_U0_vs_U3.png", dpi=300, bbox_inches='tight')
print("Muazzam! 'DOS_U0_vs_U3.png' dosyası basariyla olusturuldu.")
