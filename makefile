CXX = gcc-mp-8 #gcc-mp-8
CXXFLAGS = -ansi -lstdc++ #-W -g -Wall  -Wwrite-strings -Wshadow
OPTFLAGS = -O3 -ffast-math -framework Accelerate -fopenmp
RPTFLAGS = -fopt-info-vec-all #-fopt-report=5 -fopt-report-phase=vec -fopt-report-phase=par
#/usr/local/lib: fftw, /urs/lib: armadillo
LFLAGS =
#/usr/local/include: fftw
#Users/apicon/Documents/ProgramasC/Cpp/Libraries: Laser.h, iodata.h, Constants.h
INCLUDES = -I/opt/local/include/gcc8/c++/ -I./headers -I/usr/local/include
#-lm -fftw3: fftw, -larmadillo -framework Accelerate: armadillo
LIBS = -lm -L/opt/local/lib/ -larmadillo
NAMEXE = QuantumC.out
OBJS = QuantumC.o

.PHONY: clean all 

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
