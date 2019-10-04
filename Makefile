CC = mpicc
CFLAGS = -Wall -O3 -std=c99 -MMD
OPENMPFLAGS = -fopenmp
OPENMPL = 
LM =

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

clean:
	$(RM) -f ./src/*.o ./src/*.d ./nbody

-include $(OBJS:.o=.d)
