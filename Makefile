NVCC = nvc++

CFLAGS = -acc -Wall -O0 -g

MATHLIB = -cudalib=cufft,curand

SOURCES := $(wildcard src/*.cpp)
HEADERS := $(wildcard src/*.h)
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCES))
MAIN    := src/main.o

OUT := bin/3DMultiPhysics

$(OUT): $(OBJECTS) $(MAIN)
	$(NVCC) $(CFLAGS) $^ -o $(OUT) $(MATHLIB)

%.o: %.cpp
	$(NVCC) $(CFLAGS) -c $< -o $@ $(MATHLIB)

$(MAIN): src/main.cu
	$(NVCC) $(CFLAGS) -c src/main.cu -o src/main.o $(MATHLIB)

.phony: clean

clear: 
	rm bin/*.vtk bin/*.dat bin/*log*

clean:
	rm $(OUT) $(OBJECTS)

