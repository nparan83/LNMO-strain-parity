# LNMO-strain-parity

Data and code for

> N. Paran, *Channel-selective parity of the strain response in La₂NiMnO₆: a sign-blind majority gap and a displaced minority gap*, Computational Materials Science, manuscript COMMAT-D-26-04486 (revised version R2, November 2026).

All first-principles calculations were performed with Quantum ESPRESSO 7.2 (PBE+U, Dudarev, ortho-atomic projectors, U(Mn) = 4.0 eV, U(Ni) = 5.1 eV; PSlibrary PAW datasets La.pbe-spfn-kjpaw_psl.1.0.0, Mn.pbe-spn-kjpaw_psl.0.3.1, Ni.pbe-n-kjpaw_psl.0.1, O.pbe-n-kjpaw_psl.1.0.0; 90/900 Ry). The pseudopotential files themselves are not redistributed; they are available from PSlibrary.

## Layout

| Folder | Content |
|---|---|
| `inputs/strain_<eps>/{FM,AFM,FiM}/` | pw.x, dos.x and projwfc.x input files for all eleven strain states (ε = −0.05 … +0.05) and the FM, C-type AFM and FiM configurations (FiM at −0.05, 0.00, +0.05 only); `strain_0.00/FM/K_Test/` contains the k-point convergence tests |
| `outputs/strain_<eps>/...` | the corresponding output files (gzip-compressed) |
| `dos/` | density-of-states and projected-DOS files used for Fig. 4 |
| `data/` | derived data files of the original submission |
| `scripts/`, `scripts/figures/` | extraction, analysis and plotting scripts of the original submission (Python) |
| `revision1/` | calculations added in the revision: Γ-point phonons, elastic constants, band paths, e_g/t_2g projections, a⁻a⁻a⁻ tilt pattern, and the R2 figure script (see `revision1/README.md`) |
| `MANIFEST.txt` | list of the files of the original submission |

## From files to results (revised manuscript, R2 numbering)

| Manuscript item | Data file | Produced by |
|---|---|---|
| Sec. 3.1–3.2, Fig. 2 and inset: total energies, ΔE(AFM−FM), ΔE(FiM−FM), M_tot | `data/energies.dat` | `scripts/extract_p1.py`; figure `revision1/scripts/make_figs.py` |
| Sec. 3.3, Fig. 3: c⊥ = V/\|a×b\| = c sin β | `data/cperp.dat` | `scripts/extract_p1.py`; figure `revision1/scripts/colab_R2_figs.py` |
| Sec. 3.3, Table A2: elastic constants | `revision1/R1_final.tar.gz` → `el/Cij.dat` | `revision1/el_gen.py`, `revision1/el_post.py` |
| Sec. 3.4, Table A1: Γ-point phonons | `revision1/R1_final.tar.gz` → `ph/eps_*/gamma_freq.dat`, `gamma_irreps.txt` | `revision1/ph_gen.py`, `revision1/ph_post.py`, `revision1/ph_irr.py` |
| Table 1: Löwdin moments and populations | `data/lowdin.dat` | `scripts/extract_p1.py` |
| Table 2, Q₃ amplitudes of NiO₆ and MnO₆ | `data/q3_geo.dat` | `scripts/q3_geo.py` (uses `scripts/q3lib.py`) |
| Fig. 4: DOS and projected DOS | `dos/` | `scripts/dos_cizim.py` |
| Fig. 5, Sec. 3.5: e_g/t_2g occupations | `revision1/data/eg_t2g_occ.dat` | see `revision1/README.md` |
| Fig. 6, Sec. 2.3, Sec. 3.5: band paths and band edges | `revision1/data/band_edges_R1.dat`; inputs `revision1/bands/bands_inputs.tar.gz` | `revision1/r1_bands2.sh`, `revision1/r1_bands3.sh` |
| **Table 3: spin-resolved gaps** | **`revision1/data/gap_R2.dat`** (smaller of mesh and path value where the path was computed) | from `band_edges_R1.dat` and `data/gap.dat` |
| Fig. 7: E_g against \|Q₃\| and the linear fit | `revision1/data/gap_R2.dat`, `data/q3_geo.dat` | `revision1/scripts/make_figs.py` |
| Table 4, Fig. 8, gap maxima of Sec. 3.6 | `revision1/data/gap_R2.dat` | `revision1/scripts/make_figs.py` |
| Sec. 3.7, Fig. 9, Table A3: a⁻a⁻a⁻ versus a⁻a⁻c⁺ | `revision1/data/tilt_gap_data.dat`; `revision1/R1_final.tar.gz` → `tilt/` | `revision1/tilt_gen.py`, `revision1/tilt_post.py`, `revision1/gap_aaa.py` |
| Sec. 2: C-type spin pattern of the AFM configuration | `inputs/strain_*/AFM/` | `scripts/afm_id.py` |
| Table 5: substrate mismatch | a_pc = 3.935 Å and literature substrate parameters (manuscript Sec. 3.8) | — |

`data/gap.dat` holds the 6×6×4 mesh gaps of the original submission. In the revision the band paths showed that under tension the minority valence-band maximum lies off the high-symmetry points, so Table 3 now uses the smaller of the mesh and path values (`revision1/data/gap_R2.dat`). The checks printed by `scripts/extract_p1.py` and `scripts/figures/fig7_gap_parity.py` refer to the values of the original submission. `revision1/scripts/make_figs.py` reproduces the R2 numbers (fit E_g = 1.327 − 1.370|Q₃|, r = −0.983, RMSE = 14 meV; gap maxima at +0.03 % and −1.60 %; Table 4).

The scripts in `scripts/` were run in the original working layout (`LNMO_Strain_Production/strain_*/{FM,AFM,FiM}/`, uncompressed outputs), which corresponds one-to-one to `inputs/` and `outputs/` here. Some comments, column headers and summary files are in Turkish (yol = path, ızgara = mesh).

## Requirements

Python ≥ 3.8 with NumPy and Matplotlib; Phonopy 2.26 for the phonon scripts.

## Licence

Code: MIT (see `LICENSE`). Data: CC BY 4.0.

## Contact

Nejdet Paran, Department of Physics, Faculty of Science, Erciyes University, Kayseri, Türkiye — nejdetparan@erciyes.edu.tr — ORCID 0000-0003-0317-0438
