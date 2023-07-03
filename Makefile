CC = gcc
CFLAGS = -g -Wall -c
LDFLAGS := -lm
VERSIONS = mpi

SRC_FOLDER = src
BUILD_FOLDER = build

EXECUTABLES = $(foreach version, $(VERSIONS), $(BUILD_FOLDER)/$(version).out)
SRCS = $(foreach dir, $(VERSIONS), $(wildcard $(SRC_FOLDER)/$(dir)/*.c))
OBJS = $(patsubst $(SRC_FOLDER)/%.c, $(BUILD_FOLDER)/%.o, $(SRCS))

.PHONY: all clean

all: $(EXECUTABLES)

mpi: CC = mpicc
mpi: $(EXECUTABLES)

omp: CFLAGS += -fopenmp
omp: LDFLAGS += -fopenmp
omp: $(EXECUTABLES)

$(EXECUTABLES): $(OBJS)
	$(CC) $^ $(LDFLAGS) -o $@

$(OBJS): $(BUILD_FOLDER)/%.o: $(SRC_FOLDER)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

profileserial: $(EXECUTABLES)
	@valgrind --tool=callgrind --dump-instr=yes --simulate-cache=yes --collect-jumps=yes ./build/serial.out

clean:
	rm -rf build
