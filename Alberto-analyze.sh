rm -f all.dat
c='/orange/alberto.perezant/Reza/TF/NAR/'
for a in  1a74  1azp  1by4  1cdw  1dh3  1zme 2r1j 3cro
do
    rm -f ${a}.dat
    touch ${a}.dat
 for b in  11 12 13 14 15 16
 do
     #head $c/$a/${b}-*/contact-clusters/DNA-PROT/summary |awk '{print $3}' > dummy
     #wc -l $c/$a/${b}-*/contact-clusters/DNA-PROT/summary | awk '{print $1}' >> dummy
     head $c/$a/${b}-*/contact-clusters/PROT/summary |awk '{print $3}' > dummy
     wc -l $c/$a/${b}-*/contact-clusters/PROT/summary | awk '{print $1}' >> dummy
     paste ${a}.dat dummy > tmp
     mv tmp ${a}.dat
 done
 echo $a >> all.dat
 cat ${a}.dat >> all.dat
 done


     
     
