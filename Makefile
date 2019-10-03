CC = mpicc
CFLAGS = -Wall -O3 -std=c99 -MMD
OPENMPFLAGS = -fopenmp

ifeq ($(shell uname), Darwin)
	# for apple's clang
	OPENMPFLAGS = -Xpreprocessor -fopenmp
	CFLAGS += -Wno-format
endif

CFLAGS += $(OPENMPFLAGS)

RM = /bin/rm

SRCS = $(wildcard ./src/*.c)
OBJS = $(subst .c,.o,$(SRCS))

all: nbody

nbody: $(OBJS)
	$(CC) $(CFLAGS) -lomp -o $@ $(OBJS)

clean:
	$(RM) -f ./src/*.o ./src/*.d ./nbody

-include $(OBJS:.o=.d)
