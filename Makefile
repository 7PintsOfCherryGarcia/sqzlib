CC=gcc
CFLAGS=  -Wall -Wextra -std=c11 -pedantic -O3
CFLAGSB= -Wall -Wextra -std=c11 -pedantic -g
CSHFLAG= -shared

SRCDIR= src
INC= -Ilibs -I.
LIB= -L.
LIBS= -lz -lzstd -lpthread

SRC=sqzlib/sqz_init.c sqzlib/sqz_kseq.c sqzlib/sqz_coding.c\
    sqzlib/sqz_zlib.c sqzlib/sqz_zstd.c sqzlib/sqz_filefun.c

OBJS=$(SRC:%.c=$(SRCDIR)/%.o)
SOBJS=$(SRC:%.c=$(SRCDIR)/%S.o)

OBJSB=$(SRC:%.c=$(SRCDIR)/%B.o)
SOBJSB=$(SRC:%.c=$(SRCDIR)/%SB.o)


.PHONY:all clean wipe
.SUFFIXES:.c .o


all:libsqz pthr sqz


%.o:%.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

%S.o:%.c
	$(CC) $(CFLAGS) $(INC) -c -fPIC $< -o $@

%B.o:%.c
	$(CC) $(CFLAGSB) $(INC) -c $< -o $@

%SB.o:%.c
	$(CC) $(CFLAGSB) $(INC) -c -fPIC  $< -o $@


pthr:
	$(CC) $(CFLAGS) $(INC) -Wno-unused-function -c -o src/pthread/sqz_pthread.o src/pthread/sqz_pthread.c

pthrB:
	$(CC) $(CFLAGSB) $(INC) -Wno-unused-function -c -o src/pthread/sqz_pthreadB.o src/pthread/sqz_pthread.c


libsqzflag:
	$(info >>>>>Building libsqz<<<<<)

sqzflag:
	$(info >>>>>Building sqz<<<<<)


libsqz:libsqzflag $(OBJS)
	ar rcs $@.a $(OBJS)
	cp src/sqzlib/sqzlib.h .

libsqz_shared:libsqzflag $(SOBJS)
	$(CC) $(CSHFLAG) $(SOBJS) -o libsqz.so
	cp src/sqzlib/sqzlib.h .

libsqzB:libsqzflag $(OBJSB) pthrB
	ar rcs libsqz.a $(OBJSB)
	cp src/sqzlib/sqzlib.h .

libsqz_sharedB:libsqzflag $(SOBJSB)
	$(CC) $(CSHFLAG) $(SOBJSB) -o libsqz.so
	cp src/sqzlib/sqzlib.h .


sqz:sqzflag pthr
	$(CC) $(CFLAGS) $(INC) $(LIB) -Wno-unused-function -o $@ $(SRCDIR)/sqz.c src/pthread/sqz_pthread.o libsqz.a $(LIBS)

sqzB:sqzflag pthrB
	$(CC) $(CFLAGSB) $(INC) $(LIB) -o sqz $(SRCDIR)/sqz.c src/pthread/sqz_pthreadB.o libsqz.a $(LIBS)

build:wipe libsqzB sqzB

clean:
	rm -rf sqz libsqz.a sqzlib.h $(SRCDIR)/sqzlib/*.o $(SRCDIR)/pthread/*.o

wipe:
	rm -rf sqz libsqz.so libsqz.a sqzlib.h sqz_data.h $(SRCDIR)/sqzlib/*.o $(SRCDIR)/pthread/*.o

# makedepend line not in use in current compilation enviroanment
