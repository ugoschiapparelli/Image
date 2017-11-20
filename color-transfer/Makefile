ROOT=../bcl

CFLAGS=-Wall -Wextra -Werror -std=c99 -I$(ROOT)/include
LDFLAGS=-L$(ROOT)/lib -lbcl
LDLIBS=-lm -lfftw3 -lbcl

BIN=\
	color-transfer

all: $(BIN)

clean:
	$(RM) $(BIN) *.ppm

.PHONY: all clean
