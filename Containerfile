FROM nvcr.io/nvidia/nvhpc:22.1-devel-cuda11.5-ubuntu20.04

RUN apt-get update && \
    apt-get install apt-transport-https ca-certificates gnupg software-properties-common wget & \
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add - & \
    apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main' & \
    apt-get update & \
    apt-get install cmake

WORKDIR /Q-POP-FerroDyn

COPY ./src /Q-POP-FerroDyn/src
COPY ./include /Q-POP-FerroDyn/include
COPY ./Makefile /Q-POP-FerroDyn/Makefile
COPY ./input /Q-POP-FerroDyn/input
RUN mkdir bin

CMD ["bash"]