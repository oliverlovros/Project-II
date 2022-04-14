# Parallel version of the Makefile for the Molecular Dynamics simulation

# MPI Compiler

mpi_comp = mpif90

# MPI run

mpi_exe = mpirun

# Some path and object variables are defined

OTHERDIR_MAIN=../Code

OTHERDIR_MOD=../Code/Modules

OTHERDIR_RESULTS=../Results/Data

OTHERDIR_FIG=../Results/Figures

MOD_ORDENADOS=../Code/Modules/Parallel_module.f90 ../Code/Modules/Constants_module.f90 ../Code/Modules/Pbc_module.f90 ../Code/Modules/Initial_state_module.f90 ../Code/Modules/Forces_module.f90 ../Code/Modules/Integrator_module.f90 ../Code/Modules/Gofr_module.f90 ../Code/Main.f90


all: compile exe statistics move clean
run: compile exe

compile: $(MOD_ORDENADOS)
	$(mpi_comp) -g $(MOD_ORDENADOS) -o program.exe

exe:
	source execution.sh

# statistics: Statistical analysis is performed and the figures are generated. You will be asked some questions.
.PHONY: statistics
statistics:
	@python $(OTHERDIR_MAIN)/generate_figures.py parameters.txt
	@python $(OTHERDIR_MAIN)/Block_avarage.py parameters.txt

# move  : Results files (extensions .dat and .png) are moved to Results folders (Data and Figures).
.PHONY: move	
move: 
	@mv *.dat $(OTHERDIR_RESULTS)
	@mv *.png $(OTHERDIR_FIG)
	@mv *.xyz $(OTHERDIR_RESULTS)

# clean : Objects, executables and modules are removed.
.PHONY: clean
clean:
	@rm -f $(OTHERDIR_MAIN)/*.o
	@rm -f $(OTHERDIR_MAIN)/*.mod
	@rm -f $(OTHERDIR_MOD)/*.mod
	@rm -f *.exe
	@rm -f *.mod

# clean_all : Same as "clean" and, additionally, the results and figures are also deleted (the repository returns to its original state).
.PHONY: clean_all
clean_all:
	@rm -f $(OTHERDIR_RESULTS)/*
	@rm -f $(OTHERDIR_FIG)/*
	@rm -f $(OTHERDIR_MAIN)/*.o
	@rm -f $(OTHERDIR_MAIN)/*.mod
	@rm -f $(OTHERDIR_MOD)/*.mod
	@rm -f *.dat
	@rm -f *.xyz
	@rm -f *.png
	@rm -f *.exe
	@rm -f *.mod


# The commands are explained in the order in which it is recommended to use them (if you want to use option 2).

.PHONY : help
help :
	@sed -n 's/^##//p' Makefile