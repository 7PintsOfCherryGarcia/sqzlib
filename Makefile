CC=gcc
CFLAGS= -Wall -g -std=c11 -pedantic
CSHFLAG= -shared

SRCDIR= src
INC= -Ilibs
LIB= -L.
LIBS= -lz

PROG=libsqz sqz
SRC=sqz_kseq.c sqz_init.c sqz_coding.c
OBJS=$(SRC:%.c=$(SRCDIR)/%.o)
SOBJS=$(SRC:%.c=$(SRCDIR)/%S.o)
#SOBJS=sqz_initS.o sqz_kseqS.o sqz_codingS.o

.PHONY:all clean
.SUFFIXES:.c .o

%.o:%.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

%S.o:%.c
	$(CC) -c -fPIC $(CFLAGS) $(INC) $< -o $@

all:$(PROG)

libsqz:$(OBJS) $(SOBJS)
	ar rcs $@.a $(OBJS)
	$(CC) $(CSHFLAG) $(SOBJS) -o $@.so

sqz:
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ $(SRCDIR)/sqz.c libsqz.a $(LIBS)

clean:
	rm -rf sqz $(SRCDIR)/*.o


# makedepend line not in use in current compilation enviroanment
