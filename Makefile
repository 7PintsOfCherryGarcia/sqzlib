
CC=gcc
CFLAGS= -Wall -g -std=c11 -pedantic
CDPLYFLAGS= -Wall -std=c11 -pedantic -O2
CSHFLAG= -shared

SRCDIR= src
INC= -Ilibs
LIB= -L.
LIBS= -lz

SRC=sqz_kseq.c sqz_init.c sqz_coding.c sqz_zlib.c sqz_cmp.c sqz_dcp.c

OBJS=$(SRC:%.c=$(SRCDIR)/%.o)
SOBJS=$(SRC:%.c=$(SRCDIR)/%S.o)
DPLYOBJS=$(SRC:%.c=$(SRCDIR)/%DPLY.o)
DPLYSOBJS=$(SRC:%.c=$(SRCDIR)/%DPLYS.o)


.PHONY:all clean
.SUFFIXES:.c .o


all:libsqzDPLY sqzDPLY

%.o:%.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

%S.o:%.c
	$(CC) -c -fPIC $(CFLAGS) $(INC) $< -o $@

%DPLY.o:%.c
	$(CC) $(CDPLYFLAGS) $(INC) -c $< -o $@

%DPLYS.o:%.c
	$(CC) -c -fPIC $(CDPLYFLAGS) $(INC) $< -o $@



libsqz:$(OBJS) $(SOBJS)
	ar rcs $@.a $(OBJS)
	$(CC) $(CSHFLAG) $(SOBJS) -o $@.so

sqz:
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ $(SRCDIR)/sqz.c libsqz.a $(LIBS)

libsqzDPLY:$(DPLYOBJS) $(DPLYSOBJS)
	ar rcs libsqz.a $(DPLYOBJS)
	$(CC) $(CSHFLAG) $(DPLYSOBJS) -o libsqz.so

sqzDPLY:
	$(CC) $(CDPLYFLAGS) $(INC) $(LIB) -o sqz $(SRCDIR)/sqz.c libsqz.a $(LIBS)



build:wipe libsqz sqz

deploy:wipe libsqzDPLY sqzDPLY

clean:
	rm -rf $(SRCDIR)/*.o

wipe:
	rm -rf sqz libsqz.so libsqz.a $(SRCDIR)/*.o

# makedepend line not in use in current compilation enviroanment
