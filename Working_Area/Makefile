# define variables
OTHERDIR_MAIN=../Code

OTHERDIR_MOD=../Code/Modules

OTHERDIR_RESULTS=../Results/Data

OTHERDIR_FIG=../Results/Figures

MOD_ORDENADOS=../Code/Modules/Constants_module.f90 ../Code/Modules/Pbc_module.f90 ../Code/Modules/Initial_state_module.f90 ../Code/Modules/Forces_module.f90 ../Code/Modules/Integrator_module.f90

OBJS= $(OTHERDIR_MAIN)/constantes.o $(OTHERDIR_MAIN)/initial_state.o $(OTHERDIR_MAIN)/pbc.o $(OTHERDIR_MAIN)/forces.o $(OTHERDIR_MAIN)/integrator.o $(OTHERDIR_MAIN)/main.o

# compile -> compile the (Fortran) code (main.f90 and gofr.f90).
# run -> execute the (Fortran) code (main.f90 and gofr.f90).
# all -> compile, execute the code (all) (main.f90, gofr.f90 and python programs) and remove the useless files.
# statistics -> realize the statistical analysis and the figures (python programs).
# clean -> remove objects, executables and modules.
# clean_all -> remove the same than clean and the results and figures.

all: compile run statistics move clean
compile: constantes initial_state pbc forces integrator main gofr program program_gofr

print:
	$(info Pasos: make modulos y luego make all: corremos el programa en serie, se recomienda realizar un make clean tras usar el programa para borrar archivos temporales. Se está tarbajando en una versión en paralelo del programa, si funciona sería con make paralel.)

constantes: $(OTHERDIR_MOD)/Constants_module.f90
	gfortran -c $(OTHERDIR_MOD)/Constants_module.f90 -o $(OTHERDIR_MAIN)/constantes.o

initial_state: $(OTHERDIR_MOD)/Initial_state_module.f90 $(OTHERDIR_MOD)/Constants_module.f90
	gfortran -c $(OTHERDIR_MOD)/Initial_state_module.f90 -o $(OTHERDIR_MAIN)/initial_state.o

pbc: $(OTHERDIR_MOD)/Pbc_module.f90 $(OTHERDIR_MOD)/Constants_module.f90
	gfortran -c $(OTHERDIR_MOD)/Pbc_module.f90 -o $(OTHERDIR_MAIN)/pbc.o

forces: $(OTHERDIR_MOD)/Forces_module.f90 $(OTHERDIR_MOD)/Pbc_module.f90 $(OTHERDIR_MOD)/Constants_module.f90
	gfortran -c $(OTHERDIR_MOD)/Forces_module.f90 -o $(OTHERDIR_MAIN)/forces.o	

integrator: $(OTHERDIR_MOD)/Integrator_module.f90 $(OTHERDIR_MOD)/Forces_module.f90 $(OTHERDIR_MOD)/Pbc_module.f90 $(OTHERDIR_MOD)/Constants_module.f90
	gfortran -c $(OTHERDIR_MOD)/Integrator_module.f90 -o $(OTHERDIR_MAIN)/integrator.o

main: $(OTHERDIR_MAIN)/Main.f90 $(MOD_ORDENADOS)
	gfortran -c $(OTHERDIR_MAIN)/Main.f90 -o $(OTHERDIR_MAIN)/main.o

gofr: $(OTHERDIR_MAIN)/Gofr.f90
	gfortran -c $(OTHERDIR_MAIN)/Gofr.f90 -o $(OTHERDIR_MAIN)/gofr.o
	
program: $(OBJS)
	gfortran -o program.exe $(OBJS)

program_gofr: $(OTHERDIR_MAIN)/gofr.o
	gfortran -o program_gofr.exe $(OTHERDIR_MAIN)/gofr.o
run: 
	time ./program.exe
	./program_gofr.exe

.PHONY: move	
move: 
	@mv *.dat $(OTHERDIR_RESULTS)
	@mv *.png $(OTHERDIR_FIG)
	@mv *.xyz $(OTHERDIR_RESULTS)


.PHONY: clean
clean:
	@rm -f $(OTHERDIR_MAIN)/*.o
	@rm -f *.exe
	@rm -f *.mod
	@rm -f gofr_params.txt


.PHONY: clean_all
clean_all:
	@rm -f $(OTHERDIR_RESULTS)/*
	@rm -f $(OTHERDIR_FIG)/*
	@rm -f gofr_params.txt
	@rm -f $(OTHERDIR_MAIN)/*.o
	@rm -f *.exe
	@rm -f *mod

.PHONY: statistics
statistics:
	python $(OTHERDIR_MAIN)/generate_figures.py parameters.txt
	python $(OTHERDIR_MAIN)/Block_avarage.py parameters.txt

.PHONY: vmd
vmd:
	vmd $(OTHERDIR_RESULTS)/vmd.xyz