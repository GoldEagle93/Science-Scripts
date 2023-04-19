#!/bin/bash

#cat<<EFO>get_top.py

cpptraj topol.top<<EOF
trajin trajectory.00.dcd  10000 10001
strip :SDM parmout topol-no-dummy.top
go
EOF

# PD

cpptraj<<EOF
parm topol-no-dummy.top [system]
parm reference.prmtop [reftop]
reference reference.inpcrd parm [reftop] [ref]
trajin contact-clusters-single-traj-pd/DNA-PROT/unique.c?.pdb parm [system]

align :XXX@N,C,CA,O5',C4',C1' ref [ref]
rms top10-pd ref [ref] :XXX@N,C,CA,O5',C4',C1' nofit out rmsd-top10-pd.dat
go
EOF

# Normal

cpptraj<<EOF
parm topol-no-dummy.top [system]
parm reference.prmtop [reftop]
reference reference.inpcrd parm [reftop] [ref]
trajin contact-clusters-single-traj/DNA-PROT/unique.c?.pdb parm [system]

align @O5',C4',C1' ref [ref]
rms fit @O5',C4',C1'
rms top10 ref [ref] :XXX@N,C,CA nofit out rmsd-top10.dat
go
EOF

# LI

cpptraj<<EOF
parm topol-no-dummy.top [system]
parm reference.prmtop [reftop]
reference reference.inpcrd parm [reftop] [ref]
trajin contact-clusters-single-traj-li/DNA/unique.c?.pdb parm [system]

align :XXX@N,C,CA ref [ref]
rms fit ref [ref] :XXX@N,C,CA
rms top10-li ref [ref] :XXX@O5',C4',C1' nofit out rmsd-top10-li.dat
go
EOF


