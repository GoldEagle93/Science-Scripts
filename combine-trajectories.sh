#!/bin/bash

cp ../14-two-sided-dummy-contacts/topol.top .
cp ../14-two-sided-dummy-contacts/reference.prmtop .
cp ../14-two-sided-dummy-contacts/reference.inpcrd .
cp topol.top topol.prmtop

cpptraj topol.top<<EOF
trajin ../11-dummy-contacts-70%/trajectory.00.dcd  10000 20000
trajin ../12-dummy-contacts-70%-trust-1/trajectory.00.dcd  10000 20000
trajin ../14-two-sided-dummy-contacts/trajectory.00.dcd  10000 20000
trajout combination-11-12-14.dcd
go
EOF


