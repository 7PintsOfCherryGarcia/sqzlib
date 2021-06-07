CC=gcc
CFLAGS=  -Wall -std=c11 -pedantic -O2
CFLAGSB= -Wall -std=c11 -pedantic -g
CSHFLAG= -shared

SRCDIR= src
INC= -Ilibs
LIB= -L.
LIBS= -lz -lpthread

SRC=sqz_kseq.c sqz_init.c sqz_coding.c\
    sqz_zlib.c sqz_cmp.c sqz_dcp.c\
    sqz_filefun.c pthread/sqz_pthread.c

OBJS=$(SRC:%.c=$(SRCDIR)/%.o)
SOBJS=$(SRC:%.c=$(SRCDIR)/%S.o)

OBJSB=$(SRC:%.c=$(SRCDIR)/%B.o)
SOBJSB=$(SRC:%.c=$(SRCDIR)/%SB.o)


.PHONY:all clean wipe examples
.SUFFIXES:.c .o


all:libsqz sqz

%.o:%.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

%S.o:%.c
	$(CC) $(CFLAGS) $(INC) -c -fPIC $< -o $@

%B.o:%.c
	$(CC) $(CFLAGSB) $(INC) -c $< -o $@

%SB.o:%.c
	$(CC) $(CFLAGSB) $(INC) -c -fPIC  $< -o $@


libsqzflag:
	$(info >>>>>Building libsqz<<<<<)

sqzflag:
	$(info >>>>>Building sqz<<<<<)

exampleflag:
	$(info >>>>>Building examples<<<<<)



libsqz:libsqzflag $(OBJS)
	ar rcs $@.a $(OBJS)
	cp src/sqzlib.h .
	cp src/sqz_data.h .

libsqz_shared:libsqzflag $(SOBJS)
	$(CC) $(CSHFLAG) $(SOBJS) -o libsqz.so
	cp src/sqzlib.h .
	cp src/sqz_data.h .


libsqzB:libsqzflag $(OBJSB)
	ar rcs libsqz.a $(OBJSB)
	cp src/sqzlib.h .
	cp src/sqz_data.h .

libsqz_sharedB:libsqzflag $(SOBJSB)
	$(CC) $(CSHFLAG) $(SOBJSB) -o libsqz.so
	cp src/sqzlib.h .
	cp src/sqz_data.h .


sqz:sqzflag
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ $(SRCDIR)/sqz.c libsqz.a $(LIBS)

sqzB:sqzflag
	$(CC) $(CFLAGSB) $(INC) $(LIB) -o sqz $(SRCDIR)/sqz.c libsqz.a $(LIBS)


examples:exampleflag
	make -C examples

examplesB:exampleflag
	make build -C examples




build:wipe libsqzB sqzB #examplesB


clean:
	rm -rf $(SRCDIR)/*.o $(SRCDIR)/pthread/*.o

wipe:
	rm -rf sqz libsqz.so libsqz.a sqzlib.h sqz_data.h $(SRCDIR)/*.o $(SRCDIR)/pthread/*.o
	#make clean -C examples

# makedepend line not in use in current compilation enviroanment
