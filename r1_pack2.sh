#!/bin/bash
# r1_pack2.sh - PDOS paketi (uretim dizinlerinden) -> /arf/scratch/nparan/R1/R1_pdos_d.tar.gz
L=/arf/scratch/nparan/lnmo; mkdir -p /arf/scratch/nparan/R1
for D in $L/LNMO_DOS_Data_11_Full $L/LNMO_Strain_Production; do echo "$D"
 for s in Ni1 Mn1; do echo "  $s: $(find $D -name "*($s)_wfc#*(d)" | wc -l) dosya"; done; done
tar tzf ~/R1_pdos_d.tar.gz | sed 's|/[^/]*$||' | sort | uniq -c | head -30
cd $L/LNMO_Strain_Production && tar czf /arf/scratch/nparan/R1/R1_pdos_d.tar.gz strain_*/FM/*'(Ni1)_wfc#3(d)' \
  strain_*/FM/*'(Mn1)_wfc#4(d)' strain_*/FM/*pdos_tot strain_*/FM/scf_FM.in strain_*/FM/nscf_FM.out
ls -l /arf/scratch/nparan/R1/R1_pdos_d.tar.gz; tar tzf /arf/scratch/nparan/R1/R1_pdos_d.tar.gz | wc -l
