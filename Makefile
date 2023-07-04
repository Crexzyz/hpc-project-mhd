CC = gcc
CFLAGS = -g -Wall -c
LDFLAGS := -lm
VERSIONS = mpi

SRC_FOLDER = src
BUILD_FOLDER = build

STRATEGY :=

ifneq ($(filter serial,$(MAKECMDGOALS)),)
	STRATEGY = serial
else ifneq ($(filter mpi,$(MAKECMDGOALS)),)
	STRATEGY = mpi
else ifneq ($(filter omp,$(MAKECMDGOALS)),)
	STRATEGY = omp
else ifneq ($(filter vectorization,$(MAKECMDGOALS)),)
	STRATEGY = vectorization
endif

EXECUTABLE = $(BUILD_FOLDER)/$(STRATEGY).out

SRCS := $(wildcard $(SRC_FOLDER)/$(STRATEGY)/*.c)
OBJS := $(patsubst $(SRC_FOLDER)/%.c, $(BUILD_FOLDER)/%.o, $(SRCS))

$(EXECUTABLE): $(OBJS)
	$(CC) $^ $(LDFLAGS) -o $@

$(OBJS): $(BUILD_FOLDER)/%.o: $(SRC_FOLDER)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

serial: $(EXECUTABLE)

omp: CFLAGS += -fopenmp
omp: LDFLAGS += -fopenmp
omp: $(EXECUTABLE)

mpi: CC = mpicc
mpi: $(EXECUTABLE)

vectorization: CFLAGS += -std=gnu11 -march=native -mtune=native -fopenmp
vectorization: LDFLAGS += -O3 -fstrict-aliasing -ftree-vectorize -fopt-info-vec-optimized -fopenmp
vectorization: $(EXECUTABLE)

profileserial: $(STRATEGY)
	@valgrind --tool=callgrind --dump-instr=yes --simulate-cache=yes --collect-jumps=yes ./build/serial.out

.PHONY: clean

clean:
	rm -rf build
