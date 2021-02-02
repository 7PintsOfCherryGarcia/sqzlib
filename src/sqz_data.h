#include <zlib.h>
#include <stdint.h>

#define LOAD_SIZE 8L*1024L*1024L
#define NAME_SIZE 1L*1024L*1024L
#define B64       sizeof(uint64_t)
#define HEADLEN   8
#define NEND      128
#define CHUNK     131072
#define MAGIC     151324677

#define FQH       '@'
#define FAH       '>'
#define FQS       '+'
#define NL        '\n'

#ifdef KLIB
typedef struct kseq_t kseq_t;
#endif

#ifndef SQZLIB
typedef struct sqzfastx_t sqzfastx_t;
typedef struct sqzblock_t sqzblock_t;
#endif

#ifdef SQZLIB
/*
  "sqzfastx_t"
  libsqueezma main data loading structure. Defines the buffers and flags for
  reading sequencing data into.
*/
typedef struct {
    //file members
    const char *filename;
    gzFile     fp;
    //flags
    char       fmt;
    char       endflag; //Sequece has not completely been read into a buffer flag
    //data members
    uint64_t   offset;
    kseq_t     *seq;
    uint8_t    *seqbuffer;
    uint8_t    *qualbuffer;
    uint8_t    *namebuffer;
    uint8_t    *readbuffer;
    size_t     namesize;
    size_t     namepos;

    size_t     n;
    size_t     bases;
    //miscelaneous
    size_t     rem;     //Length of sequence remaining to be read
    size_t     toread;  //Size of sequence still needed to be read
    size_t     prevlen; //Size of sequence currently being read
} sqzfastx_t;


typedef struct {
    //Code data buffer
    uint8_t   *blkbuff; //Buffer to hold encoded fastX data
    uint64_t   blksize; //Size of blkbuff
    uint64_t   blkpos;  //Position within blkbuff
    uint64_t   namepos; //Position within namebuffer TODO Rethink this memeber
    char       newblk;  //flag
    //Compression members
    uint8_t   *cmpbuff;
    uint64_t   cmpsize;
    uint64_t   cmppos;
} sqzblock_t;

#endif

typedef struct {
    FILE       *fp;
    sqzfastx_t *sqz;
    sqzblock_t *blk;
    uint64_t   size;
    uint8_t    ff;
    uint64_t   filepos;
} sqz_File;



