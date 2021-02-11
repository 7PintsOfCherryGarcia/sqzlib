CC=gcc
CFLAGS= -Wall -std=c11 -pedantic -O2
CFLAGSB= -Wall -g -std=c11 -pedantic

INC= -I../libs -I..
LIB= -L..
LIBS= -lsqz -lz


.PHONY:all clean
.SUFFIXES:.c .o


all:countseqs seqstats minimap2

countseqsflag:
	$(info >>>>>Building countseqs<<<<<)

countseqs:countseqsflag
	$(info countseqs)
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ countseqs.c $(LIBS)

countseqsB:countseqsflag
	$(info countseqs)
	$(CC) $(CFLAGSB) $(INC) $(LIB) -o countseqs countseqs.c $(LIBS)


seqstatsflag:
	$(info >>>>>Building seqstats<<<<<)

seqstats:seqstatsflag
	$(info seqstats)
	make -C seqstats

seqstatsB:seqstatsflag
	$(info seqstats)
	make build -C seqstats


minimap2flag:
	$(info >>>>>Building minimap2<<<<<)

minimap2:minimap2flag
	$(info minimap2)
	make -C minimap2


clean:
	rm  -rf countseqs
	make -C seqstats clean
	make -C minimap2 clean


build:clean countseqsB seqstatsB minimap2

# makedepend line not in use in current compilation enviroanment
