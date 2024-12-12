cmake_minimum_required(VERSION 3.20)
set(CMAKE_BUILD_TYPE Release)
set (CMAKE_CXX_STANDARD 11)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 80)
endif()

project(Q-POP-FerroDyn VERSION 1.0
                       DESCRIPTION "Q-POP-FerroDyn OpenACC"
                       LANGUAGES CXX CUDA)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  set(CMAKE_CXX_EXTENSIONS OFF)

  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)

file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "src/*.h")

set(MAIN "src/main.cu")
set_source_files_properties(src/main.cu PROPERTIES LANGUAGE CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(3DMultiPhysics ${SOURCES} ${MAIN})

set_target_properties(3DMultiPhysics PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

target_compile_options(3DMultiPhysics PRIVATE $<$<COMPILE_LANGUAGE:CXX>:
  -acc
  -Wall
  -O0
  -g
  -gpu=cc80,cuda11.5
  -cudalib=cufft,curand
>)

target_link_options(3DMultiPhysics PRIVATE $<$<COMPILE_LANGUAGE:CXX>:
  -acc
  -Wall
  -O0
  -g
  -gpu=cc80,cuda11.5
  -cudalib=cufft,curand
>)

add_custom_target(clear
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/*.vtk
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/*.dat
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/*.log
)


