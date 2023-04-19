#!/bin/bash

# deactivate
# meld

pdbid=$1
export pdbid

cat <<- EOF > amberparam
 source leaprc.protein.ff14SBonlysc
 aa = loadpdb $pdbid-stripped.pdb
 saveamberparm aa $pdbid.top $pdbid.crd
 quit
EOF

tleap -f amberparam
rm amberparam

# cat <<- EOF > refgen
# cpptraj $pdbid.top<<EOF
#  trajin $pdbid.crd
#  trajout $pdbid-ref.pdb
#  go
#  quit
# EOF

cat <<- EOF > refgen
 parm $pdbid.top
 trajin $pdbid.crd
 trajout $pdbid-ref.pdb
 go
 quit
EOF

cpptraj -i refgen
rm refgen