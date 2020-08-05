#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#define LOAD_SIZE 64*1024

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


/*
  "sqz_fastxinit"
*/
sqzfastx_t *sqz_fastxinit(const char *filename, size_t bsize);

/*
  "sqz_getformat"
  Reads the first sequence record in a fast q/a file. If a quality string is found, fastq
  format is assumed, fasta otherwise
*/
char sqz_getformat(const char *filename);

/*
  "sqz_kseqinit"
  Starts kseq objects
*/
char sqz_kseqinit(sqzfastx_t *sqz);

/*
  "sqz_kill"
  frees allocated buffers in sqzfastx_t struct
*/
void sqz_kill(sqzfastx_t *sqz);


/*
  "sqz_loadfastq"
  loads fastq data to sqzfastx_t object, returns number of sequences loaded
*/
size_t sqz_loadfastq(sqzfastx_t *sqz);


/*
    "sqz_sqzblkinit"
    Starts sqzblock_t object
*/
sqzblock_t *sqz_sqzblkinit(size_t size);


char sqz_encode(sqzfastx_t *sqz, sqzblock_t *blk);


void sqz_blkdestroy(sqzblock_t *blk);


size_t sqz_deflate(sqzblock_t *blk, int level);


char sqz_zlibcmpdump(sqzblock_t *blk, size_t size, FILE *ofp);
