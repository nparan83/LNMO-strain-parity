#!/bin/bash
# r1_pack3.sh - bant kosulari kontrol + paket; O 2p PDOS paketi -> /arf/scratch/nparan/R1/ (hesap yok)
cd /arf/scratch/nparan/lnmo/LNMO_Strain_Production
for e in -0.05 0.00 0.05; do B=strain_$e/FM/R1_bands
 echo "== $e: k=$(grep -c 'Computing kpt' $B/bands.out) done=$(grep -c 'JOB DONE' $B/bands.out) $(grep 'PWSCF.*WALL' $B/bands.out | tail -1)"
 echo "   hs=$(grep -c 'high-symmetry' $B/bandsx_s1.out) gnu: $(ls $B/*.gnu 2>/dev/null | xargs -n1 basename | tr '\n' ' ')"; done
tar czf /arf/scratch/nparan/R1/R1_bands.tar.gz strain_*/FM/R1_bands/{bands.in,bandsx_s1.out,bandsx_s2.out,bands_s1.dat.gnu,bands_s2.dat.gnu}
tar czf /arf/scratch/nparan/R1/R1_pdos_O.tar.gz strain_*/FM/*'(O)_wfc#2(p)'
ls -l /arf/scratch/nparan/R1/R1_bands.tar.gz /arf/scratch/nparan/R1/R1_pdos_O.tar.gz
tar tzf /arf/scratch/nparan/R1/R1_bands.tar.gz | wc -l; tar tzf /arf/scratch/nparan/R1/R1_pdos_O.tar.gz | wc -l
