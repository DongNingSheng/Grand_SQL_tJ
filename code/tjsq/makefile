.SUFFIX:
.SUFFIXES: .o .mod .F .f .F90 .f90
F90FLAGS= -g -O2 -xhost -qopenmp -mcmodel=medium -qmkl=parallel -traceback
FC= ifort 
DIR= $1
TAG= tjsq.out
OBJS=common.o  block.o get_sys.o den_mat.o  eigen.o dmrg_sub.o  measure.o main.o
TAG = tjsq.out

all : 	$(TAG)

$(TAG):	$(OBJS) 
	$(FC)  $(OBJS) -qopenmp -mcmodel=medium -qmkl=parallel -o $(TAG)  
.f90.o:
	$(FC) -c $< $(F90FLAGS)
clean: 
	rm -f *.o *.mod core 
