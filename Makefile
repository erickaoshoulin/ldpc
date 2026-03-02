CC ?= gcc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -pedantic
LDFLAGS ?= -lm

SRC := src/main.c src/ldpc.c
OBJ := $(SRC:.c=.o)
TARGET := ldpc_eval

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -Iinclude -o $@ $(OBJ) $(LDFLAGS)

src/%.o: src/%.c include/ldpc.h
	$(CC) $(CFLAGS) -Iinclude -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJ) $(TARGET) results/performance.csv results/performance.png results/performance.svg

.PHONY: all run clean
