# Q-POP-FerroDyn

## Installation

The single-GPU version of this module uses OpenACC to offload computations to the GPU. The upcoming exascale version will use Kokkos to offload computing, and MPI to scale across nodes.

Required packages:
1. Nvidia HPC-SDK (22.1+, but 22.1 is preferred)
   - Provides OpenACC via nvc++
   - Required CUDA libraries: cuFFT, cuRand
2. HDF5
   - Present on most supercomputers as a module

CMake will automatically detect packages if they have been installed correctly.
```
cmake -B build .
cd build
make
```

It may be more convenient for some Perlmutter users to use a container to compile and run Q-POP-FerroDyn.
Navigate to the root directory of the repository. The Dockerfile located here can be used to build a minimal container.
```
podman-hpc build -t q-pop-ferrodyn:latest .
podman-hpc migrate q-pop-ferrodyn:latest
podman-hpc run --rm -it -gpu -v $PWD:/ferrodyn q-pop-ferrodyn:latest /bin/bash
cd /ferrodyn
make
cd bin
./Q-POP-FerroDyn
```

For documentation on podman-hpc and containerization, please refer https://docs.nersc.gov/development/containers/podman-hpc/overview/.
