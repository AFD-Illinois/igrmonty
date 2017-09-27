CC = h5pcc
CCFLAGS = -std=c99 -Ofast -Wall -Werror -fdiagnostics-color -fopenmp -L/home/brryan/Software/gsl/lib/ -I/home/brryan/Software/gsl/include/
LIB_FLAGS = -lm -ldl -lgsl -lgslcblas
CC_COMPILE = $(CC) $(CCFLAGS) -c
CC_LOAD = $(CC) $(CCFLAGS)
.c.o:
	$(CC_COMPILE) $*.c
EXE = grmonty
all: $(EXE)
SRC = compton.c geodesics.c hotcross.c init_geometry.c interp.c jnu_mixed.c main.c malloc.c model_bhlight3d.c radiation.c random.c scatter_super_photon.c tetrads.c track_super_photon.c utils.c 
OBJ = compton.o geodesics.o hotcross.o init_geometry.o interp.o jnu_mixed.o main.o malloc.o model_bhlight3d.o radiation.o random.o scatter_super_photon.o tetrads.o track_super_photon.o utils.o 
INC = constants.h decs.h 
$(OBJ) : $(INC) makefile
$(EXE): $(OBJ) $(INC) makefile
	$(CC_LOAD) $(OBJ) $(LIB_FLAGS) -o $(EXE)