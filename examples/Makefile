CC=gcc
CFLAGS= -Wall -std=c11 -pedantic -O3

INC= -I../libs -I..
LIBS= -lz -lzstd


.PHONY:all clean
.SUFFIXES:.c .o


all:countseqs

countseqsflag:
	$(info >>>>>Building countseqs<<<<<)


countseqs:countseqsflag
	$(info countseqs)
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ countseqs.c ../libsqz.a $(LIBS)


clean:
	rm  -rf countseqs

