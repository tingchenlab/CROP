# Project: CROPLinux

CPP  = g++ -D_DEBUG_
CC   = gcc -D_DEBUG_
RES  = 
OBJ  = alignment.o Unique.o bayesianclustering.o CROP.o Main.o $(RES)
LINKOBJ  = alignment.o Unique.o bayesianclustering.o CROP.o Main.o $(RES)
LIBS =  -g3 -lgsl -lgslcblas -fopenmp -m64
INCS =  
CXXINCS =  
BIN  = CROPLinux
CXXFLAGS = $(CXXINCS)   -g3 -fopenmp -m64
CFLAGS = $(INCS)   -g3 -fopenmp -m64
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before CROPLinux all-after


clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o CROPLinux $(LIBS)

alignment.o: alignment.c alignment.h common.h
	$(CPP) -c alignment.c -o alignment.o $(CXXFLAGS)

Unique.o: Unique.cpp Unique.h common.h
	$(CPP) -c Unique.cpp -o Unique.o $(CXXFLAGS)

bayesianclustering.o: bayesianclustering.cpp bayesianclustering.h alignment.h common.h
	$(CPP) -c bayesianclustering.cpp -o bayesianclustering.o $(CXXFLAGS)

CROP.o: CROP.cpp CROP.h bayesianclustering.h common.h
	$(CPP) -c CROP.cpp -o CROP.o $(CXXFLAGS)

Main.o: Main.cpp CROP.h Unique.h
	$(CPP) -c Main.cpp -o Main.o $(CXXFLAGS)