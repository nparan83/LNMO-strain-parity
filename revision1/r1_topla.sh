#!/bin/bash
# r1_topla.sh - CMS R1 sonuclarini topla: elastik kontrol + C_ij, tilt tablosu, fonon/irrep, a-a-a- E_g; ozet + Colab paketi
R=/arf/scratch/nparan/R1; cd $R || exit 1; O=$R/R1_ozet.txt
{
echo "### $(date '+%F %T')  R1 ozet"
echo "== 1. Elastik: JOB DONE / BFGS bitti / son adim sayisi"
for i in $(seq 1 12); do f=el/$i/pw.out
 printf ' %2d tag=%-6s done=%s bfgs_end=%s adim=%s\n' $i "$(cat el/$i/tag 2>/dev/null)" "$(grep -c 'JOB DONE' $f 2>/dev/null)" \
  "$(grep -c 'End of BFGS' $f 2>/dev/null)" "$(grep -c '^!' $f 2>/dev/null)"; done
echo "== 2. C_ij"; python3 el_post.py $R/el 2>&1 | grep -v Deprecation
echo "== 3. Tilt (aaa - aac)"; python3 tilt_post.py 2>&1 | grep -v Deprecation
echo "== 4. Fonon (en dusuk 6 / en yuksek 3, cm-1) + irrep"
for e in -0.05 0.00 0.05; do W=ph/eps_$e; [ -f $W/gamma_irreps.txt ] || python3 ph_irr.py $W > /dev/null 2>&1
 echo " $e: $(awk '!/^#/{print $3}' $W/gamma_freq.dat | head -6 | tr '\n' ' ') | $(awk '!/^#/{print $3}' $W/gamma_freq.dat | tail -3 | tr '\n' ' ')"; done
echo "== 5. a-a-a- E_g"; python3 gap_aaa.py -0.03 -0.04 -0.05 2>&1
} | tee $O
tar czf $R/R1_final.tar.gz R1_ozet.txt el/Cij.dat ph/eps_*/gamma_freq.dat ph/eps_*/gamma_irreps.txt \
 tilt/eps_*/scf.in tilt/eps_*/gap/{bands.in,bandsx_s1.out,bandsx_s2.out,bands_s1.dat.gnu,bands_s2.dat.gnu} 2>/dev/null
ls -l $R/R1_final.tar.gz; tar tzf $R/R1_final.tar.gz | wc -l
