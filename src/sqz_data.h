#include <zlib.h>
#include <stdint.h>

#define LOAD_SIZE 16*1024*1024

#define NEND 128

#define CHUNK 131072

#define MAGIC 151324677

#ifndef KLIB
#define kseq_t struct kseq_t
#endif



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
    size_t     bases;
    //flags
    char       fmt;
    char       endflag; //Sequece has not completely been read into a buffer flag
    //miscelaneous
    size_t     rem;     //Length of sequence remaining to be read
    size_t     toread;  //Size of sequence still needed to be read
    size_t     prevlen; //Size of sequence currently being read
} sqzfastx_t;


typedef struct {
    //Code data buffer
    uint8_t *codebuff;
    size_t  blksize;
    char    newblk;
    //Compression members
    uint8_t *cmpbuff;
    size_t cmpsize;
} sqzblock_t;


