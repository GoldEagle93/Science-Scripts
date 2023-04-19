#deactivate; conda deactivate; module purge; ml cuda/11.0.207 gcc/9.3.0 openmpi/4.0.4 mkl/2020.0.166 meld/0.4.19; source /home/alberto.perezant/Source/amber/amber.sh
#python /blue/alberto.perezant/reza/Scripts/prepare-RMSD.py
#deactivate; conda deactivate; module purge; module restore Amber20
# pdb4amber -i ref.pdb -o ref-pdb4amber.pdb
# tleap -f leap.in
#python /blue/alberto.perezant/reza/Scripts/prepare-RMSD.py
# module load python/3.8
# cpptraj -i strip.cpptraj
for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
do
cpptraj -i rmsd$i.cpptraj
done
paste rmsd??.dat | awk '{print $1 " " $2 " " $4 " " $6 " " $8 " " $10 " " $10 " " $12 " " $14 " " $16 " " $18 " " $20 " " $22 " " $24 " " $26 " " $28 " " $30 " " $32 " " $34 " " $36 " " $38 " " $40 " " $42 " " $44 " " $46 " " $48 " " $50 " " $52 " " $54 " " $56 " " $58}' > rmsd.dat
python /blue/alberto.perezant/reza/Scripts/plot-RMSD.py
# rm rmsd?.dat
# rm rmsd?.cpptraj
# rm strip.cpptraj
