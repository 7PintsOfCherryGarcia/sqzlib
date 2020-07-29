#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#define LOAD_SIZE 16*1024

/*
  "sqzfastx_t"
  libsqueezma main data loading structure. Defines the buffers and flags for
  reading sequencing data into.
*/
typedef struct {
    const char *filename;
    gzFile     fp;
    kseq_t     *seq;
    uint8_t    *seqbuffer;
    size_t     seqlen;
    uint8_t    *namebuffer;
    size_t     namelen;
    uint8_t    *qualbuffer;
    char       fmt;
    size_t     n;
    char       endflag;  //Sequece has not completely been read into a buffer flag
    size_t     rem;      //Length of sequence remaining to be read
    size_t     toread;   //Size of sequence still needed to be read
    size_t     prevlen;  //Size of sequence currently being read
} sqzfastx_t;

char sqz_getformat(const char *filename);


char sqz_kseqinit(sqzfastx_t *sqz);


size_t sqz_loadfastq(sqzfastx_t *sqz);


size_t sqz_newblock(sqzfastx_t *sqz);


size_t sqz_endblock(sqzfastx_t *sqz);


void sqz_kill(sqzfastx_t *sqz);
