
# this bakes the current version into the executable
GIT_VERSION = $(shell git describe --dirty --always --tags)

CC = h5pcc  # or h5cc.  enables compilation with hdf5

# change DMODEL if you use a non-iharm model
CCFLAGS = -shlib -DMODEL=iharm -DNOTES= -DVERSION=$(GIT_VERSION) -O3 -Wall -fopenmp 

LIB_FLAGS = -lm -lgsl -lgslcblas  # gsl is used for special functions, integrations, and
                                  # random number generation

CC_COMPILE = $(CC) $(CCFLAGS) -c
CC_LOAD = $(CC) $(CCFLAGS)

.c.o:
	$(CC_COMPILE) $*.c

EXE = grmonty

all: $(EXE)

# to change the model, copy the appropriate file out of model/ into model.c
SRC = compton.c geodesics.c h5io.c hdf5_utils.c hotcross.c init_geometry.c interp.c jnu_mixed.c main.c malloc.c model.c par.c radiation.c random.c scatter_super_photon.c tetrads.c track_super_photon.c utils.c 
OBJ = compton.o geodesics.o h5io.o hdf5_utils.o hotcross.o init_geometry.o interp.o jnu_mixed.o main.o malloc.o model.o par.o radiation.o random.o scatter_super_photon.o tetrads.o track_super_photon.o utils.o 
INC = constants.h custom.h decs.h h5io.h hdf5_utils.h model.h par.h 

$(OBJ) : $(INC) makefile
$(EXE): $(OBJ) $(INC) makefile
	$(CC_LOAD) $(OBJ) $(LIB_FLAGS) -o $(EXE)

clean:
	rm *.o grmonty
