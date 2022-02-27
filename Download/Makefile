OBJS= mod.o main.o
all: program run estadistica vmd
print:
	$(info Con make all corremos el programa en serie, se recomienda realizar un make clean tras usar el programa para borrar archivos temporales. Se está tarbajando en una versión en paralelo del programa, si funciona sería con make paralel.)




ejecutable: modulo.f90 P2.f90
	gfortran modulo.f90 P2.f90 -o ejecutable.exe
	time ./ejecutable.exe
	python Block_avarage.py P2_rho08.dat


paralel: P2.f90 subroutines.f90
	export OMP_NUM_THREADS=8
	gfortran -fopenmp P2.f90 subroutines.f90
	time ./a.out
	python Block_avarage.py P2_rho08.dat
serie:  P2.f90 subroutines.f90 
	gfortran P2.f90 subroutines.f90
	time ./a.out
	python Block_avarage.py P2_rho08.dat








program: $(OBJS)
	gfortran -o program.exe $(OBJS)

mod.o: modulo.f90
	gfortran -c modulo.f90 -o mod.o
main.o: P2.f90 modulo.f90
	gfortran -c P2.f90 -o main.o 
run:
	time ./program.exe

clean:
	rm *o *gz *png *dat *xyz program a.out

estadistica:
	python Block_avarage.py P2_rho08.dat
vmd:
	vmd vmd.xyz

