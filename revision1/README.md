# revision1 — calculations added in the revision (R2)

Inputs, summary outputs and scripts of the calculations requested by Reviewer #3. All calculations use the production settings of the main data set (Quantum ESPRESSO 7.2, PBE+U, ortho-atomic projectors, U(Mn) = 4.0 eV, U(Ni) = 5.1 eV, 90/900 Ry, 6×6×4 mesh) and the relaxed structures in `../inputs/strain_*/FM/`. The shell scripts were run on the TRUBA cluster (SLURM); paths inside them refer to that file system. Some comments and summary labels are in Turkish (yol = path, ızgara = mesh, bitmemiş = unfinished, kararlı = stable).

## Files

| File | Content |
|---|---|
| `R1_ozet.txt` | one-page summary printed on TRUBA after all runs: elastic runs, C_ij, tilt comparison, Γ-point frequencies, a⁻a⁻a⁻ gaps |
| `R1_final.tar.gz` | `el/Cij.dat`; `ph/eps_{-0.05,0.00,0.05}/gamma_freq.dat`, `gamma_irreps.txt`; `tilt/eps_*/scf.in` (relaxed a⁻a⁻a⁻ structures); `tilt/eps_{-0.03,-0.04,-0.05}/gap/` (bands.in, bands.x outputs, band files) |
| `inputs.tar.gz` | input templates of the added calculations |
| `bands/bands_inputs.tar.gz` | bands.x / pw.x inputs of the high-symmetry paths of the a⁻a⁻c⁺ structure |
| `../data/gap_R2.dat` | **Table 3**: spin-resolved gaps; smaller of the mesh and path values where the path was computed (ε = −5, 0, +1 … +5 %), mesh values otherwise |
| `../data/band_edges_R1.dat` | path, mesh and minimum gaps and the band-edge locations (Sec. 2.3, Sec. 3.5, Fig. 6) |
| `../data/eg_t2g_occ.dat` | e_g/t_2g occupations of Ni1 and Mn1 integrated to mid-gap, and the e_g–t_2g mixing ("leak") caused by the octahedral tilts (Sec. 2.3, Sec. 3.5, Fig. 5) |
| `../data/tilt_gap_data.dat` | ΔE(a⁻a⁻a⁻ − a⁻a⁻c⁺) and the gaps of both tilt patterns (Sec. 3.7, Fig. 9, Table A3) |

## Scripts

| Calculation | Manuscript item | Generation / submission | Post-processing |
|---|---|---|---|
| Γ-point phonons, finite displacements (Phonopy 2.26.2, 0.01 Å, 30 displacements per strain) | Sec. 2.4, Sec. 3.4, Table A1 | `ph_gen.py`, `r1_ph.sh` | `ph_post.py`, `ph_irr.py` |
| Elastic constants, ±0.5 % Voigt deformations with relaxed ions | Sec. 2.4, Sec. 3.3, Table A2 | `el_gen.py`, `r1_el.sh` | `el_post.py` |
| a⁻a⁻a⁻ tilt pattern under the epitaxial constraint | Sec. 2.5, Sec. 3.7, Fig. 9, Table A3 | `tilt_gen.py`, `mkin.py`, `r1_tilt.sh`, `r1_tilt_rs.sh`; gaps `r1_aaa_gap.sh` | `tilt_post.py`, `gap_aaa.py` |
| High-symmetry band paths of the a⁻a⁻c⁺ structure | Sec. 2.3, Sec. 3.5, Fig. 6 | `r1_bands2.sh`, `r1_bands3.sh` | `../scripts/cell_bands2.py` |
| e_g/t_2g projections in the pseudocubic frame | Sec. 2.3, Sec. 3.5, Fig. 5 | projwfc.x outputs of the main data set (`r1_pack2.sh`) | `../scripts/cell_egt2g.py` |
| Figures 2, 7, 8 and the graphical abstract | Fig. 2, Fig. 7, Fig. 8, Table 4 | — | `../scripts/make_figs.py` (reads `../data/` and `../data/gap_R2.dat`) |
| Figures 2, 3, 7, 8 (Google Colab version) | Fig. 2, Fig. 3, Fig. 7, Fig. 8 | — | `../scripts/colab_R2_figs.py` (reads the data files from this repository) |

`python3 scripts/make_figs.py`, run from any folder, writes the figures to `./figs/` and prints the fit E_g = 1.327 − 1.370|Q₃| (r = −0.983, RMSE = 14 meV), the gap maxima (+0.03 % and −1.60 %) and the odd parts of Table 4.
