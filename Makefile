NVCC = nvc++

CFLAGS = -acc -fast -Minfo=accel -Mcudalib=cufft,curand -gpu=cuda11.7

MATHLIB = -I/opt/nvidia/hpc_sdk/Linux_x86_64/22.7/math_libs/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/math_libs/lib64

SOURCES := $(wildcard *.cpp)
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCES))

OUT := 3DMultiPhysics

all: $(OBJECTS) main.o
        $(NVCC) $(CFLAGS) $^ -o $(OUT)

%.o: %.cpp
        $(NVCC) $(CFLAGS) -c $< -o $@
        $(NVCC) $(CFLAGS) -c main.cu -o main.o

main.o: main.cu
        $(NVCC) $(CFLAGS) -c main.cu -o main.o

.phony: clean

clean:
        rm $(OUT) *.o