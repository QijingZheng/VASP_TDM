#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
FC= ifort -assume byterecl
FC= gfortran
FFLAGS= -g -O2
MAKE = make

#-------------------------------------------------------------------------------
# Src
#-------------------------------------------------------------------------------

SRC= prec.f90 lattice.f90 info.f90 wave.f90 option.f90 gvector.f90 \
	 tdm.f90 main.f90
SRC_PARGAMMA= prec.f90 lattice.f90 info.f90 wave.f90 option.f90 gvector_gam.f90 \
	 tdm_gam.f90 main.f90

OBJ = $(SRC:.f90=.o)
OBJ_PARGAMMA = $(SRC_PARGAMMA:.f90=.o)
EXE = vasptdm

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(FC) $(FFLAGS) -c $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
tdm:	$(OBJ)
	$(FC) $(FFLAGS) -o $(EXE) $(OBJ) $(SPGLIB)  
gam: $(OBJ_PARGAMMA)
	$(FC) $(FFLAGS) -o vasptdm_gam $(OBJ) $(SPGLIB)  

clean:
	rm -f *.mod *.a
	rm -f $(OBJ) $(EXE) vasptdm_gam $(OBJ_PARGAMMA)
tar:
	tar -czvf tdm.tgz *.f90 Makefile *.py
tag:
	ctags *.f90
