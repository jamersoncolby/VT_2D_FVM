FC=gfortran

SRC=set_precision.f90 MMS_functions_2D.f90 boundary.f90 muscl.f90 vanleer.f90 roe.f90 new_residual.f90 extract.f90 fvm_main.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(FC) -o $@ -c $<

fvm_main: $(OBJ)
	$(FC) -o $@ $(OBJ)

clean:
	@rm *.o *.mod fvm_main
