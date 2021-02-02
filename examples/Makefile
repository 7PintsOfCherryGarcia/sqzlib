
CC=gcc
CFLAGS= -Wall -std=c11 -pedantic -O2
CFLAGSB= -Wall -g -std=c11 -pedantic

INC= -I../libs -I..
LIB= -L..
LIBS= -lz


.PHONY:all clean seqstats
.SUFFIXES:.c .o


all:seqstats countseqs


countseqs:
	$(info countseqs)
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ countseqs.c ../libsqz.a $(LIBS)


seqstats:
	$(info seqstats)
	make -C seqstats


countseqsB:
	$(info countseqs)
	$(CC) $(CFLAGSB) $(INC) $(LIB) -o countseqs countseqs.c ../libsqz.a $(LIBS)


seqstatsB:
	$(info seqstats)
	make build -C seqstats


clean:
	rm  -rf countseqs
	make -C seqstats clean


build:clean countseqsB seqstatsB

# makedepend line not in use in current compilation enviroanment
