CC=gcc
CFLAGS= -Wall -std=c11 -pedantic -O3

INC= -I../libs -I..
LIB= -L..
LIBS= -lsqz -lz


.PHONY:all clean
.SUFFIXES:.c .o


all:countseqs

countseqsflag:
	$(info >>>>>Building countseqs<<<<<)

countseqs:countseqsflag
	$(info countseqs)
	$(CC) $(CFLAGS) $(INC) $(LIB) -o $@ countseqs.c $(LIBS)



clean:
	rm  -rf countseqs


# makedepend line not in use in current compilation enviroanment
