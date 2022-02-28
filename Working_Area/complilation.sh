#!/bin/bash

# Main.f90 (parece una locura no poner ../Code/Modules/*.f90 pero el problema es que unos modulos llaman a otros y es importante el orde de compilacion)
gfortran ../Code/Modules/Constants_module.f90 ../Code/Modules/Pbc_module.f90 ../Code/Modules/Initial_state_module.f90 ../Code/Modules/Forces_module.f90 ../Code/Modules/Integrator_module.f90 ../Code/Main.f90 -o Main.out
# Gofr.f90
gfortran ../Code/Gofr.f90 -o Gofr.out
# borrar la basura
rm *.mod