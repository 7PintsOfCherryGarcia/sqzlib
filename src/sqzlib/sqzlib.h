/** \file sqzlib.h
    TODO: Add low level API interface to initialize a input files
    TODO: Do not write N block length when it's 127. 127 is implied
          if not last N blck
    TODO: Change member names of sqzblock_t to reflect what is coded data
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include <zlib.h>

#define LOAD_SIZE 4L*1024L*1024L
#define NAME_SIZE 1L*1024L*1024L
#define B64       sizeof(uint64_t)
#define HEADLEN   8
#define NEND      128
#define MAGIC     151324677

#define FQH       '@'
#define FAH       '>'
#define FQS       '+'
#define NL        '\n'


/*
  "sqzfastx_t"
  libsqueezma main data loading structure. Defines the buffers and flags for
  reading sequencing data into.
*/
typedef struct {
    uint8_t    nthread;
    uint8_t    endthread;
    //flags
    char       endflag; //Sequece has not completely been read into a buffer flag
    char       cmpflag;
    //data members
    uint64_t   offset;
    uint8_t    *seqbuffer;
    uint8_t    *qualbuffer;
    uint8_t    *namebuffer;
    uint8_t    *readbuffer;
    uint64_t    namesize;
    uint64_t    namepos;
    //Partially loaded sequence and qualities
    char       *pseq; //previous seq
    char       *pqlt; //previous quality
    uint64_t    plen;
    //return members
    uint64_t    n;
    uint64_t    bases;
    //miscelaneous
    uint64_t    seqread;
    uint64_t    rem;     //Length of sequence remaining to be read
    uint64_t    toread;  //Size of sequence still needed to be read
    uint64_t    prevlen; //Size of sequence currently being read
} sqzfastx_t;


/*
 Docstring
*/
typedef struct {
    //Code data buffer
    uint8_t   *blkbuff; //Buffer to hold encoded fastX data
    uint64_t   blksize; //Size of blkbuff
    uint64_t   mblksize;//Max size of blkbuff
    uint64_t   blkpos;  //Position within blkbuff
    uint64_t   namepos; //Position within namebuffer TODO Rethink this memeber
    uint64_t   prevlen; //How much of current sequence had been decoded
    uint8_t    newblk;  //flag
    //Compression members
    uint8_t   *cmpbuff; //Buffer to hold compressed blk data
    uint64_t   cmpsize;
    uint64_t   mcmpsize;//Max size of cmpbuff
    uint64_t   cmppos;
} sqzblock_t;



/*
  Docstring
*/
typedef struct {
    FILE       *fp;
    sqzfastx_t *sqz;
    sqzblock_t *blk;
    uint64_t   size;
    uint8_t    ff;
    uint64_t   filepos;
} sqz_File;


/*
##############################################################################
 sqz file routines
##############################################################################
*/

/*
  Write an sqz header
*/
char sqz_filehead(unsigned char fmt, FILE *ofp);


/*
  Write an sqz tail
*/
char sqz_filetail(uint64_t numseqs, FILE *ofp);


/*
  Get format of a fastX file
*/
uint8_t  sqz_getformat(const char *filename);


/*
  Get size of an sqz file
*/
uint64_t sqz_filesize(FILE *fp);


/*
  ##############################################################################
  sqz structure routines
  ##############################################################################
*/

/*
  Initialize an sqzblock structure
*/
sqzblock_t *sqz_sqzblkinit(uint64_t size);


/*
  Terminate an sqzblock structure
*/
void sqz_blkkill(sqzblock_t *blk);


/*
  ##############################################################################
  sqz data routines
  ##############################################################################
*/

/*
  Read a data block into an sqz block structure
*/
char sqz_readblksize(sqzblock_t *blk, FILE *fp);


/*
  Decode data in sqz block to memory
*/
uint64_t
sqz_fastXdecode(sqzblock_t *blk, uint8_t *buff, uint64_t size,char fqflag);
