#!/bin/bash

cpptraj topol.top<<EOF
# trajin trajectory.00.dcd  0 20000 100
trajin trajectory.00.dcd  0 20000 200
trajin trajectory.01.dcd  0 20000 3000
trajin trajectory.02.dcd  0 20000 3000
trajin trajectory.03.dcd  0 20000 3000
trajin trajectory.04.dcd  0 20000 3000
trajin trajectory.05.dcd  0 20000 3000
trajin trajectory.06.dcd  0 20000 3000
trajin trajectory.07.dcd  0 20000 3000
trajin trajectory.08.dcd  0 20000 3000
trajin trajectory.09.dcd  0 20000 3000
trajin trajectory.10.dcd  0 20000 3000
trajin trajectory.11.dcd  0 20000 3000
trajin trajectory.12.dcd  0 20000 3000
trajin trajectory.13.dcd  0 20000 3000
trajin trajectory.14.dcd  0 20000 3000
trajin trajectory.15.dcd  0 20000 3000
trajin trajectory.16.dcd  0 20000 3000
trajin trajectory.17.dcd  0 20000 3000
trajin trajectory.18.dcd  0 20000 3000
trajin trajectory.19.dcd  0 20000 3000
trajin trajectory.20.dcd  0 20000 3000
trajin trajectory.21.dcd  0 20000 3000
trajin trajectory.22.dcd  0 20000 3000
trajin trajectory.23.dcd  0 20000 3000
trajin trajectory.24.dcd  0 20000 3000
trajin trajectory.25.dcd  0 20000 3000
trajin trajectory.26.dcd  0 20000 3000
trajin trajectory.27.dcd  0 20000 3000
trajin trajectory.28.dcd  0 20000 3000
trajin trajectory.29.dcd  0 20000 3000
trajout skip.dcd
go
EOF


