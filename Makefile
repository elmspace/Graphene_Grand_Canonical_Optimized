LFLAGS = -lm -lfftw3  -O3 -lfftw3_threads -fopenmp -lpthread
DEBUG = -g
LIBS = -lm -lstdc++ -lfftw3

main: main.cpp
	g++ -std=c++11 $(LFLAGS) -o $@ $(MOBLIB) $(SUBFILES) main.cpp $(LIBS)
