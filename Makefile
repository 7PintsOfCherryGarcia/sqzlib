
CC=gcc
CFLAGS= -Wall -g -std=c11 -pedantic -O2
CFLAGSB= -Wall -g -std=c11 -pedantic
CSHFLAG= -shared

SRCDIR= src
INC= -Ilibs
LIB= -L.
LIBS= -lz

SRC=sqz_kseq.c sqz_init.c sqz_coding.c sqz_zlib.c sqz_cmp.c sqz_dcp.c

OBJS=$(SRC:%.c=$(SRCDIR)/%.o)
SOBJS=$(SRC:%.c=$(SRCDIR)/%S.o)

OBJSB=$(SRC:%.c=$(SRCDIR)/%B.o)
SOBJSB=$(SRC:%.c=$(SRCDIR)/%SB.o)


.PHONY:all clean wipe examples
.SUFFIXES:.c .o


all:libsqz sqz examples

%.o:%.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

%S.o:%.c
	$(CC) $(CFLAGS) $(INC) -c -fPIC $< -o $@

%B.o:%.c
	$(CC) $(CFLAGSB) $(INC) -c $< -o $@

%SB.o:%.c
	$(CC) $(CFLAGSB) $(INC) -c -fPIC  $< -o $@


libsqz:libsqzflag $(OBJS) $(SOBJS)
	ar rcs $@.a $(OBJS)
	$(CC) $(CSHFLAG) $(SOBJS) -o $@.so
	cp src/sqzlib.h .
	cp src/sqz_data.h .


sqz:sqzflag
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ $(SRCDIR)/sqz.c libsqz.a $(LIBS)


examples:exampleflag
	make -C examples


libsqzB:libsqzflag $(OBJSB) $(SOBJSB)
	ar rcs libsqz.a $(OBJSB)
	$(CC) $(CSHFLAG) $(SOBJSB) -o libsqz.so
	cp src/sqzlib.h .
	cp src/sqz_data.h .

libsqzflag:
	$(info ****Building libsqz****)


sqzB:sqzflag
	$(CC) $(CFLAGSB) $(INC) $(LIB) -o sqz $(SRCDIR)/sqz.c libsqz.a $(LIBS)

sqzflag:
	$(info ****Building sqz****)

examplesB:exampleflag
	make build -C examples

exampleflag:
	$(info ****Building examples****)


build:wipe libsqzB sqzB examplesB


clean:
	rm -rf $(SRCDIR)/*.o

wipe:
	rm -rf sqz libsqz.so libsqz.a $(SRCDIR)/*.o
	make clean -C examples

# makedepend line not in use in current compilation enviroanment
