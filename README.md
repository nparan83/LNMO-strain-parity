# LNMO-strain-parity

Data and code for

> N. Paran, *Channel-selective parity of the strain response in La₂NiMnO₆: a sign-blind majority gap and a displaced minority gap*, submitted to Computational Materials Science (manuscript COMMAT-D-26-04486).

All first-principles calculations were performed with Quantum ESPRESSO 7.2 (PBE+U, Dudarev, ortho-atomic projectors, U(Mn) = 4.0 eV, U(Ni) = 5.1 eV; PSlibrary PAW datasets La.pbe-spfn-kjpaw_psl.1.0.0, Mn.pbe-spn-kjpaw_psl.0.3.1, Ni.pbe-n-kjpaw_psl.0.1, O.pbe-n-kjpaw_psl.1.0.0; 90/900 Ry). The pseudopotential files themselves are not redistributed; they are available from PSlibrary.

## Layout

| Folder | Content |
|---|---|
| `inputs/strain_<eps>/{FM,AFM,FiM}/` | pw.x, dos.x and projwfc.x input files for all eleven strain states (ε = −0.05 … +0.05) and the FM, C-type AFM and FiM configurations (FiM at −0.05, 0.00, +0.05 only); `strain_0.00/FM/K_Test/` contains the k-point convergence tests |
| `outputs/strain_<eps>/...` | the corresponding output files (gzip-compressed) |
| `dos/` | density-of-states and projected-DOS files used for Fig. 5 |
| `data/` | derived data files from which the tables and figures are generated |
| `scripts/` | extraction and analysis scripts (Python) |
| `scripts/figures/` | plotting scripts |
| `MANIFEST.txt` | list of all files |

## From files to results

| Manuscript item | Data file | Produced by |
|---|---|---|
| Sec. 3.1–3.2, Figs. 2–3: total energies, ΔE(AFM−FM), ΔE(FiM−FM), M_tot | `data/energies.dat` | `scripts/extract_p1.py` (from `outputs/strain_*/{FM,AFM,FiM}/scf*.out.gz`) |
| Sec. 3.3, Fig. 4: c⊥ = V/\|a×b\| | `data/cperp.dat` | `scripts/extract_p1.py` (from `inputs/strain_*/FM/scf_FM.in`) |
| Table 1: Löwdin moments and populations | `data/lowdin.dat` | `scripts/extract_p1.py` (from the projwfc.x outputs) |
| Table 2, Q₃ amplitudes of NiO₆ and MnO₆ | `data/q3_geo.dat` | `scripts/q3_geo.py` (uses `scripts/q3lib.py`) |
| Table 3: spin-resolved gaps | `data/gap.dat` | `scripts/gap3.py`; independent re-extraction with fixed band indices (n↑ = 88, n↓ = 78): `data/gap_check.dat` from `scripts/extract_p1.py` |
| Sec. 3.4: band-edge k-points; Sec. 2: k-mesh test | outputs in `outputs/strain_*/FM/` and `outputs/strain_0.00/FM/K_Test/` | `scripts/gk.py`, `scripts/gk2.py`, `scripts/edges.py`, `scripts/kedge.py` |
| Fig. 5: DOS and projected DOS | `dos/` | `scripts/dos_cizim.py` |
| Fig. 6: E_g against \|Q₃\| and the linear fit | `data/gap.dat`, `data/q3_geo.dat` | `scripts/gap.py` |
| Table 4, Fig. 7, peak positions of Sec. 3.5 | `data/gap.dat` | `scripts/figures/fig7_gap_parity.py` |
| Sec. 2: C-type spin pattern of the AFM configuration | `inputs/strain_*/AFM/` | `scripts/afm_id.py` |
| Table 5: substrate mismatch | a_pc = 3.935 Å and literature substrate parameters (manuscript Sec. 3.6) | — |

`scripts/extract_p1.py` and `scripts/figures/fig7_gap_parity.py` run directly on this repository (from its root) and print a comparison with the values quoted in the manuscript; at the time of submission all 17 checks of `extract_p1.py` and all 9 checks of `fig7_gap_parity.py` pass. The remaining scripts were run in the original working layout (`LNMO_Strain_Production/strain_*/{FM,AFM,FiM}/`, uncompressed outputs), which corresponds one-to-one to `inputs/` and `outputs/` here; some of their comments are in Turkish.

Rounding note: the band gaps are given to 0.1 meV, so several even/odd components of Table 4 fall exactly on a half-way value (e.g. T↑ at |ε| = 1 %: (1.2953 − 1.2764)/2 = 9.45 meV). The manuscript rounds these half-up (9.5 meV); Python's binary floating-point formatting may print the value 0.1 meV lower (9.4 meV). The same applies to the majority-gap maximum, +0.035 % on the 0.005 % grid of the script, quoted as +0.04 % in the manuscript.

## Requirements

Python ≥ 3.8 with NumPy and Matplotlib.

## Licence

Code: MIT (see `LICENSE`). Data: CC BY 4.0.

## Contact

Nejdet Paran, Department of Physics, Faculty of Science, Erciyes University, Kayseri, Türkiye — nejdetparan@erciyes.edu.tr — ORCID 0000-0003-0317-0438
