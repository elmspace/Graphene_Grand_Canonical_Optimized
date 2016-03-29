LFLAGS = -lm -lfftw3  -O3 -lfftw3_threads -fopenmp -llapack -lblas -lgfortran
DEBUG = -g
LIBS = -lm -lstdc++ -lfftw3 -llapack

main: main.cpp
	g++ -std=c++0x $(LFLAGS) -o $@ $(MOBLIB) $(SUBFILES) main.cpp $(LIBS)
