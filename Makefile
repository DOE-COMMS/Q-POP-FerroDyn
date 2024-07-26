NVCC = nvc++
H5CC = h5c++

H5FLAGS = -I./external/HighFive/build/install/include
INCFLAGS =  -I./include
CFLAGS = -acc -Wall -O0 -g -gpu=cc86,cuda11.5

MATHLIB = -cudalib=cufft,curand

SOURCES := $(wildcard src/*.cpp)
HEADERS := $(wildcard src/*.h)
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCES))
MAIN    := src/main.o
H5      := src/write_h5.o

OUT := bin/3DMultiPhysics

$(OUT): $(OBJECTS)
	$(NVCC) $(CFLAGS) $(H5FLAGS) $(INCFLAGS) $^ -o $(OUT) $(MATHLIB)

%.o: %.cpp
	$(NVCC) $(CFLAGS) $(H5FLAGS) $(INCFLAGS) -c $< -o $@ $(MATHLIB)

H5: src/io/write_h5.cpp
	$(H5CC) $(H5FLAGS) $(INCFLAGS) -c src/io/write_h5.cpp -o src/write_h5.o

# $(MAIN): src/main.cu
# 	$(NVCC) $(CFLAGS) $(INCFLAGS) -c src/main.cu -o src/main.o

.phony: clean

clear: 
	rm bin/*.vtk bin/*.dat bin/*log*

clean:
	rm $(OUT) $(OBJECTS)

