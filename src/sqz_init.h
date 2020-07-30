#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#define LOAD_SIZE 16*1024*1024

typedef struct  kseq_t kseq_t;


/*
  "sqzfastx_t"
  libsqueezma main data loading structure. Defines the buffers and flags for
  reading sequencing data into.
*/
typedef struct {
    //file members
    const char *filename;
    gzFile     fp;
    //data members
    size_t     offset;
    kseq_t     *seq;
    uint8_t    *seqbuffer;
    uint8_t    *qualbuffer;
    uint8_t    *namebuffer;
    size_t     namelen;
    size_t     n;
    //flags
    char       fmt;
    char       endflag; //Sequece has not completely been read into a buffer flag
    //miscelaneous
    size_t     rem;     //Length of sequence remaining to be read
    size_t     toread;  //Size of sequence still needed to be read
    size_t     prevlen; //Size of sequence currently being read
} sqzfastx_t;


sqzfastx_t *sqz_fastxinit(const char *filename, size_t buffersize);


char sqz_getformat(const char *filename);


char sqz_kseqinit(sqzfastx_t *sqz);
