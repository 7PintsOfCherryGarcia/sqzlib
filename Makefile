CC=gcc
CFLAGS= -Wall -O3 -g -std=c11 -pedantic

SRC=src
INC=-Ilibs
LIB=-lz

PROG=sqz
OBJS= squeezma.o

.PHONY:all clean
.SUFFIXES:.c .o

%.o:$(SRC)/%.c
	$(CC) -c $(CFLAGS) $(INC) $< -o $@

all:$(PROG)

sqz:$(OBJS)
	$(CC) $(CFLAGS) $(INC) -o $@ $(SRC)/sqz.c $(OBJS) $(LIB)

clean:
	rm -rf zqz *.o


# makedepend line not in use in current compilation enviroanment
