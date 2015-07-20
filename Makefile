LFLAGS = -lm -lfftw3  -O3
DEBUG = -g
LIBS = -lm -lstdc++ -lfftw3

main: main.cpp
	g++ -std=c++11 $(LFLAGS) -o $@ $(MOBLIB) $(SUBFILES) main.cpp $(LIBS)

