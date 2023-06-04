CC = gcc
CFLAGS = -g -c
LDFLAGS = -lm
VERSIONS = serial

SRC_FOLDER = src
BUILD_FOLDER = build

EXECUTABLES = $(foreach version, $(VERSIONS), $(BUILD_FOLDER)/$(version).out)
SRCS = $(foreach dir, $(VERSIONS), $(wildcard $(SRC_FOLDER)/$(dir)/*.c))
OBJS = $(patsubst $(SRC_FOLDER)/%.c, $(BUILD_FOLDER)/%.o, $(SRCS))

.PHONY: all clean

all: $(EXECUTABLES)

$(EXECUTABLES): $(OBJS)
	$(CC) $^ $(LDFLAGS) -o $@

$(OBJS): $(BUILD_FOLDER)/%.o: $(SRC_FOLDER)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf build
