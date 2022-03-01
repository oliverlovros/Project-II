OTHERDIR_MAIN=../Code

OTHERDIR_MOD=../Code/Modules

OTHERDIR_RESULTS=../Results/Data

MOD_ORDENADOS=../Code/Modules/Constants_module.f90 ../Code/Modules/Pbc_module.f90 ../Code/Modules/Initial_state_module.f90 ../Code/Modules/Forces_module.f90 ../Code/Modules/Integrator_module.f90

OBJS= $(OTHERDIR_MAIN)/constantes.o $(OTHERDIR_MAIN)/initial_state.o $(OTHERDIR_MAIN)/pbc.o $(OTHERDIR_MAIN)/forces.o $(OTHERDIR_MAIN)/integrator.o $(OTHERDIR_MAIN)/main.o


modulos: constantes initial_state pbc forces integrator 
all: main program run mover
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
	
program: $(OBJS)
	gfortran -o program.exe $(OBJS)
run: 
	time ./program.exe
mover: 
	mv MSD.dat $(OTHERDIR_RESULTS)
	mv temp_evol.dat $(OTHERDIR_RESULTS)	
	mv positions.dat  $(OTHERDIR_RESULTS)     


.PHONY: clean
clean:
	rm $(OTHERDIR_MAIN)/*o
	rm *mod


.PHONY: clean_all
clean_all:
	rm $(OTHERDIR_MAIN)/*o
	rm $(OTHERDIR_RESULTS)/*dat 
	rm *mod

estadistica:
	python $(OTHERDIR_MAIN)/Block_avarage.py $(OTHERDIR_RESULTS)/temp_evol.dat
vmd:
	vmd vmd.xyz

