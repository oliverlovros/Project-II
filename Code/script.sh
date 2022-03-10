#!/bin/sh
gfortran -o P2.out Modules/Constants_module.f90 Modules/Pbc_module.f90 Modules/Initial_state_module.f90 Modules/Forces_module.f90 Modules/Integrator_module.f90 P2.f90
./P2.out
rm P2.out
rm *.mod