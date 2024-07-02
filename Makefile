NVCC = nvc++

CFLAGS = -acc -gpu=cc86,cuda11.8 -Wall -O0 -g -stdpar

MATHLIB = -lcufft -lcurand

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

clean:
	rm $(OUT) $(OBJECTS)
