CC=gcc
CFLAGS= -Wall -std=c11 -pedantic -O3

INC= -I../libs -I..
LIBS= -lz -lzstd


.PHONY:all clean
.SUFFIXES:.c .o


all:countseqs printseqs

countseqsflag:
	$(info >>>>>Building countseqs<<<<<)

printseqsflag:
	$(info >>>>>Building printseqs<<<<<)

countseqs:countseqsflag
	$(info countseqs)
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ countseqs.c ../libsqz.a $(LIBS)

printseqs:printseqsflag
	$(info printseqs)
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ printseqs.c ../libsqz.a $(LIBS)


clean:
	rm  -rf countseqs


# makedepend line not in use in current compilation enviroanment
