#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include <stdlib.h>
#include <zlib.h>

#define LOAD_SIZE 8*1024*1024

typedef struct kseq_t kseq_t;


typedef struct {
  //file members
  const char *filename;
  gzFile     fp;
  //flags
  char       fmt;
  char       endflag; //Sequece has not completely been read into a buffer flag
  //data members
  size_t     offset;
  kseq_t     *seq;
  uint8_t    *seqbuffer;
  uint8_t    *qualbuffer;
  uint8_t    *namebuffer;
  size_t     namelen;
  size_t     maxname;
  size_t     n;
  size_t     bases;
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
    Initialize sqz struct
*/
sqzfastx_t *sqz_fastxinit(const char *filename, uint64_t bsize);

/*
    Free sqz object
*/
void sqz_kill(sqzfastx_t *sqz);

/*
    Initialize blk struct
*/
sqzblock_t *sqz_sqzblkinit(size_t size);

/*
    Write sqz file header
*/
char sqz_filehead(sqzfastx_t *sqz, FILE *ofp);


/*
    Load fastX data into sqz struct
*/
size_t sqz_loadfasta(sqzfastx_t *sqz);


/*
    Encode fastX data already loaded in sqz struct and sotore it in blk struct
*/
char sqz_fastaencode(sqzfastx_t *sqz, sqzblock_t *blk);

/*
    Compress data already loaded in blk struct
*/
size_t sqz_deflate(sqzblock_t *blk, int level);


/*
    
*/
uint64_t sqz_loadfastq(sqzfastx_t *sqz);
