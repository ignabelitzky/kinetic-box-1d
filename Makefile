CC := gcc
CFLAGS := -Wall -Werror -Wextra -pedantic -std=c99 -O3 -march=native
CPPFLAGS := -D_POSIX_C_SOURCE=200809L
CLIBS := -lm
OMP := -fopenmp

TARGET := main

SRC := src/g1D-sp-075d.c src/utils.c
OBJ := $(SRC:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(OMP) -o $@ $^ $(CLIBS)

%.o: %.c include/utils.h include/constants.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(OMP) -c -o $@ $<

format:
	clang-format -style=Microsoft -i src/*.c include/*.h

clean:
	rm -f $(TARGET) $(OBJ)

.PHONY: all clean format test

test: $(TARGET)
	python3 tests/test_simulation.py ./$(TARGET)
