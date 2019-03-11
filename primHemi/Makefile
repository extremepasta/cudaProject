# operating system
HOST_OS := $(shell uname -s 2>/dev/null | tr "[:upper:]" "[:lower:]")

# architecture 
ARCH := $(shell getconf LONG_BIT)

CUDA_HOME = /usr/local/cuda-9.1

NVCC := $(CUDA_HOME)/bin/nvcc

ifeq ($(HOST_OS),darwin)
	CXX := clang++
else
	CXX := g++
endif

STD := --std=c++11

CXX_FLAGS  := $(STD) -I../ -I../../

CUDA_INC = $(CUDA_HOME)/include 

HEMI_INC =  /home/pesco/Scrivania/HEMIGPU/primHemi

BOOST_INC = /home/pesco/Scrivania/HEMIGPU/primHemi/boost_1_68_0

HEMI_HEADERS = $(HEMI_INC)/hemi/*.h

HOST_ONLY_FLAGS := -DHEMI_CUDA_DISABLE

CXX_FLAGS += -O3

CXX_FLAGS += $(DEBUG_FLAGS)

NVCCFLAGS += -g -arch=sm_30 --std=c++11 --expt-extended-lambda

cpp_sources := quantiliMain.cpp

cu_sharedlibs := $(patsubst %.cu, %.o,$(cu_sources))

cpp_sharedlibs := $(patsubst %.cpp, %.o, $(cpp_sources))

OBJECTS = $(cpp_sharedlibs) 

all : quantiliMain.o

%.o: %.cpp $(cpp_sources)
		nvcc -x cu $(NVCCFLAGS) -I$(HEMI_INC) -I$(BOOST_INC) -I$(CUDA_HOME)/include $^ -c -o $@
