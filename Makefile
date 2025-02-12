NVCC = nvc++

INCFLAGS =  -I./include
CFLAGS = -acc -Wall -fast -gpu=cc80,cuda11.5

MATHLIB = -cudalib=cufft,curand

SOURCES := $(wildcard src/*.cpp)
HEADERS := $(wildcard include/*.h)
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCES))
MAIN    := src/main.o

OUT := bin/3DMultiPhysics

$(OUT): $(OBJECTS) $(MAIN)
	$(NVCC) $(CFLAGS) $(INCFLAGS) $^ -o $(OUT) $(MATHLIB)

%.o: %.cpp
	$(NVCC) $(CFLAGS) $(INCFLAGS) -c $< -o $@ $(MATHLIB)

$(MAIN): src/main.cu
	$(NVCC) $(CFLAGS) $(INCFLAGS) -c src/main.cu -o src/main.o

.phony: clean

clear: 
	rm bin/*.vtk bin/*.dat bin/*log*

clean:
	rm $(OUT) $(OBJECTS) $(MAIN)

