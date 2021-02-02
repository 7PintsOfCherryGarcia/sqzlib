
CC=gcc

.PHONY:seqstats seqstatsB

LIB= -lz -lm
INC= -I../..

seqstats: seqstats.c klib/kseq.h klib/ksort.h
	$(CC) $< -Wall -O3 $(INC)  -o $@ ../../libsqz.a $(LIB)


seqstatsB: seqstats.c klib/kseq.h klib/ksort.h
	$(CC) $< -Wall $(INC)  -o seqstats ../../libsqz.a $(LIB)

build:clean seqstatsB


clean:
	rm -rf seqstats
