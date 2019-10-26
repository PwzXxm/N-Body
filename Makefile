#Authors:
#    - Weizhi Xu   (weizhix)  752454
#    - Zijun Chen  (zijunc3)  813190

CXX = mpic++
CXXFLAGS = -Wall -std=c++11 -O3 -MMD
LDFLAGS =

nbody-with-cuda: CXXFLAGS += -DWITH_CUDA

ifeq ($(shell uname), Linux)
	CXXFLAGS += -Wno-unused-result -fopenmp
	LDFLAGS += -lm
endif

ifeq ($(shell uname), Darwin)
	CXXFLAGS += -Xpreprocessor -fopenmp
	LDFLAGS += -lomp
endif

RM = /bin/rm

SRCS = $(wildcard ./src/*.cpp)
OBJS = $(subst .cpp,.o,$(SRCS))

all: nbody

nbody: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

nbody-with-cuda: $(OBJS) cuda
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) ./src/cuda.o $(LDFLAGS) -lcudart -L/usr/local/cuda-10.1/lib64/

cuda: ./src/cuda.cu
	nvcc -c ./src/cuda.cu -DHIDE_MPI -o ./src/cuda.o

clean:
	$(RM) -f ./src/*.o ./src/*.d ./nbody ./nbody-with-cuda

-include $(OBJS:.o=.d)
