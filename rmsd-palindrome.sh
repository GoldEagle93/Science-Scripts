for i in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
do
cpptraj -i rmsd$i-p.cpptraj
done
paste rmsd??-p.dat | awk '{print $1 " " $2 " " $4 " " $6 " " $8 " " $10 " " $10 " " $12 " " $14 " " $16 " " $18 " " $20 " " $22 " " $24 " " $26 " " $28 " " $30 " " $32 " " $34 " " $36 " " $38 " " $40 " " $42 " " $44 " " $46 " " $48 " " $50 " " $52 " " $54 " " $56 " " $58}' > rmsd.dat

