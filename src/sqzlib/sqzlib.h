/** \file sqzlib.h
    TODO: Add low level API interface to initialize a input files
    TODO: Do not write N block length when it's 127. 127 is implied
          if not last N blck
    TODO: Change member names of sqzblock_t to reflect what is coded data
    TODO: Significant work is needed in the decoding low level API. Too much
          code redundancy and therefore inefficiencies.
    TODO: Change quality binning to mimic qualbin.c from htsbox
    TODO: Incorporate reading modes to sqzFile
    TODO: Set sqzrewind to only work in reading mode
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>


#define LOAD_SIZE 4L*1024L*1024L   //Sequence buffer size
#define NAME_SIZE 1L*1024L*1024L   //Sequence name buffer size
#define B64       8                //64bits - 8 bytes
#define HEADLEN   8
#define NEND      128
#define MAGIC     151324677        //4 bytes: 5 8 5 9

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
    uint8_t     nthread;
    uint8_t     endthread;
    //flags
    char        endflag; //Sequece has not completely been read into a buffer flag
    char        cmpflag;
    //data members
    uint64_t    offset;
    uint8_t     *seq;
    uint8_t     *qlt;
    uint8_t     *namebuffer;
    uint8_t     *readbuffer;
    uint64_t    namesize;
    uint64_t    namepos;
    //return members
    uint64_t    n;
    uint64_t    bases;
    uint64_t    blks;
    //Partially loaded sequence
    uint8_t     *pseq;
    uint8_t     *pqlt;
    uint64_t    psize;
    //Partially decoded sequences
    uint64_t    seqread; //Amount of sequence read.
    uint64_t    rem;     //Length of sequence remaining to be read
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
typedef struct sqzFile_s {
    FILE       *fp;
    gzFile      gzfp;
    sqzfastx_t *sqz;
    sqzblock_t *blk;
    uint64_t    size;
    uint8_t     ff;
    uint8_t     fmt;
    uint8_t     libfmt;
    uint64_t    filepos;
} *sqzFile;


/*
  ##############################################################################
  kseq compatibility routines
  ##############################################################################
*/
/*
  Read len decoded and decompressed bytes into buff
  Returns number of bytes written into buff. If returned value is less than
  len but greater that or equal to 0, end of file has been reached. -1 on error.
*/
int64_t sqzread(sqzFile file, void *buff, uint64_t len);


/*
  Open sqzFile
  Returns sqzFile object or NULL on failure
*/
sqzFile sqzopen(const char *filename, const char *mode);


/*
  Associate an sqzFile with the file descriptor fd
*/
sqzFile sqzdopen(int fd, const char *mode);


/*
 Rewinds sqzFile. Only suported during reading
*/
void sqzrewind(sqzFile file);


/*
  Close sqzFile
*/
void sqzclose(sqzFile file);


/*
##############################################################################
 sqz file routines
##############################################################################
*/

/*
  Write an sqz header
*/
char sqz_filehead(uint8_t fmt, uint8_t libfmt, FILE *ofp);


/*
  Write an sqz tail
*/
char sqz_filetail(uint64_t numseqs, uint64_t numblks, FILE *ofp);


/*
  Get format of a fastX file
*/
uint8_t  sqz_getformat(sqzFile sqzfp);


/*
  Get size of an sqz file
*/
uint64_t sqz_filesize(sqzFile sqzfp);

/*
  Get number of blocks in sqz file
*/
uint8_t sqz_getblocks(sqzFile fp, uint64_t *b);


/*
  ##############################################################################
  sqz structure routines
  ##############################################################################
*/

/*
  Initialize sqzfastx structure
*/
sqzfastx_t *sqz_fastxinit(uint8_t fmt, uint64_t size);


/*
  Free sqzfastx structure
*/
void sqz_fastxkill(sqzfastx_t *sqz);

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
char sqz_readblksize(sqzblock_t *blk, FILE *fp, uint8_t libfmt);


/*
  Compress encoded data block
*/
int64_t sqz_blkcompress(sqzblock_t *blk, int level, uint8_t libfmt);


/*
  Encode data from sqz into blk
*/
char sqz_fastXencode(sqzfastx_t *sqz, sqzblock_t *blk, uint8_t fqflag);


/*
  Decode data in sqz block to memory
*/
uint64_t sqz_fastXdecode(sqzblock_t *blk,
                         uint8_t *buff,
                         uint64_t size,
                         uint8_t fqflag);


uint64_t sqz_loadfastX(sqzfastx_t *sqz, uint8_t fqflag, void *seq);


char sqz_blkdump(void *cmpbuff, uint64_t *blksize, uint64_t cmpsize, FILE *ofp);
