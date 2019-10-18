CXX = mpic++
CXXFLAGS = -Wall -std=c++11 -O3 -MMD -fopenmp
LDFLAGS =

nbody-with-cuda: CXXFLAGS += -DWITH_CUDA

ifeq ($(shell uname), Linux)
	CXXFLAGS += -Wno-unused-result
	LDFLAGS += -lm
endif

RM = /bin/rm

SRCS = $(wildcard ./src/*.cpp)
OBJS = $(subst .cpp,.o,$(SRCS))

all: nbody

nbody: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LM)

nbody-with-cuda: $(OBJS) cuda
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) ./src/cuda.o $(LM) -lcudart -L/usr/local/cuda-10.1/lib64/

cuda: ./src/cuda.cu
	nvcc -c ./src/cuda.cu -o ./src/cuda.o -I/usr/lib/mpich/include -L/usr/lib/mpich/lib

clean:
	$(RM) -f ./src/*.o ./src/*.d ./nbody ./nbody-with-cuda

-include $(OBJS:.o=.d)
