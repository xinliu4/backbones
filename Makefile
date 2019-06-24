#CPP = /usr/local/mpich-install/bin/mpicxx
CPP = mpicxx
CPPFLAGS = -fopenmp 
LIBS = /projects/DarkUniverse_esp/childh/hacc2/cooley.thrust/mpi/lib/libGenericIOMPI.a -L/usr/lib64 -lhdf5 -lhdf5_cpp -lgsl -lgslcblas
INCLUDES = -I /projects/DarkUniverse_esp/childh/hacc2/genericio -I/usr/include

all: backbones.exe clean

backbones.exe : backbones.o
	$(CPP) $(CPPFLAGS) -o backbones.exe backbones.o $(LIBS)

backbones.o: backbones.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c backbones.cpp -o backbones.o

clean: 
	rm *.o 
