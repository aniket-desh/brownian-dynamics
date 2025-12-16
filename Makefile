CC=gcc
CFLAGS=-std=c11 -O2 -Wall -Wextra -pedantic
INCLUDES=-Iinclude
SRC=$(wildcard src/*.c)
OBJ=$(SRC:.c=.o)
TARGET=bd_sim
LDLIBS=-lm

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
