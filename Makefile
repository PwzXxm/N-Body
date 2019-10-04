CC = mpicc
CFLAGS = -Wall -O3 -std=c99 -MMD
OPENMPFLAGS = -fopenmp
OPENMPL = 
LM =

nbody-with-cuda: CFLAGS += -DWITH_CUDA

ifeq ($(shell uname), Darwin)
	# for apple's clang
	OPENMPFLAGS = -Xpreprocessor -fopenmp
	OPENMPL = -lomp
	CFLAGS += -Wno-format
endif

ifeq ($(shell uname), Linux)
	CFLAGS += -Wno-unused-result -Wno-format-security
	LM = -lm
endif

CFLAGS += $(OPENMPFLAGS)

RM = /bin/rm

SRCS = $(wildcard ./src/*.c)
OBJS = $(subst .c,.o,$(SRCS))

all: nbody

nbody: $(OBJS)
	$(CC) $(CFLAGS) $(OPENMPL) -o $@ $(OBJS) $(LM)

nbody-with-cuda: $(OBJS) cuda
	mpic++ $(CFLAGS) $(OPENMPL) -o $@ $(OBJS) ./src/cuda.o $(LM) -lcudart -L/usr/local/cuda-10.1/lib64/

cuda: ./src/cuda.cu
	nvcc -c ./src/cuda.cu -o ./src/cuda.o -I/usr/lib/mpich/include -L/usr/lib/mpich/lib

clean:
	$(RM) -f ./src/*.o ./src/*.d ./nbody ./nbody-with-cuda

-include $(OBJS:.o=.d)
