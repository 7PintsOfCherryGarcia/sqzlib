CC=gcc
CFLAGS= -Wall -g -std=c11 -pedantic
CSHFLAG= -shared

SRC= src
INC= -Ilibs
LIB= -L.
LIBS= -lz

PROG=libsqz sqz
OBJS=sqz_init.o
SOBJS=sqz_initS.o

.PHONY:all clean
.SUFFIXES:.c .o

%.o:$(SRC)/%.c
	$(CC) -c $(CFLAGS) $(INC) $< -o $@

%S.o:$(SRC)/%.c
	$(CC) -c -fPIC $(CFLAGS) $(INC) $< -o $@

all:$(PROG)

libsqz:$(OBJS) $(SOBJS)
	ar rcs $@.a $(OBJS)
	$(CC) $(CSHFLAG) $(SOBJS) -o $@.so

sqz:
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ $(SRC)/sqz.c libsqz.a $(LIBS)

clean:
	rm -rf sqz *.o


# makedepend line not in use in current compilation enviroanment
