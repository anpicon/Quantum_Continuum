CXX = g++ #gcc-mp-8
CXXFLAGS = -ansi -std=c++14 -pg #-W -g -Wall  -Wwrite-strings -Wshadow -lstdc++ 
OPTFLAGS = -O3 -ffast-math -fopenmp # -O3 -framework Accelerate
RPTFLAGS = -fopt-info-vec-all #-fopt-report=5 -fopt-report-phase=vec -fopt-report-phase=par
#/usr/local/lib: fftw, /urs/lib: armadillo
LFLAGS = -L/home/emilioerc/gsl/lib -L/usr/lib
#/usr/local/include: fftw
#Users/apicon/Documents/ProgramasC/Cpp/Libraries: Laser.h, iodata.h, Constants.h
INCLUDES = -I./headers -I/usr/local/include #-I/opt/local/include/gcc8/c++/
#-lm -fftw3: fftw, -larmadillo -framework Accelerate: armadillo
LIBS = -lm -L/opt/local/lib/ -larmadillo
NAMEXE = QuantumC.x
OBJS = QuantumC.o

.PHONY: exec clean #clean all 

all: ending clean
	@echo All programs compiled

ending: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LFLAGS) $(OBJS) $(INCLUDES) $(LIBS) -o $(NAMEXE)
	@echo Linked

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(RPTFLAGS) $(INCLUDES) -c $<
	@echo Compiled

clean: 
	$(RM) $(OBJS)
# DO NOT DELETE
