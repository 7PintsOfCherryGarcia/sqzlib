#include <zlib.h>

#define LOAD_SIZE 4L*1024L*1024L
#define NAME_SIZE 1L*1024L*1024L
#define B64       sizeof(uint64_t)
#define HEADLEN   8
#define NEND      128
#define MAGIC     151324677

#define NBLK      255U
#define QBLK        0U
#define TNBLK      15U

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
    //Partially decoded sequences
    uint64_t    seqread; //Amount of sequence read.
    uint64_t    rem;     //Length of sequence remaining to be read
    uint64_t    toread;  //Size of sequence still needed to be read
    uint64_t    prevlen; //Size of sequence currently being read
} sqzfastx_t;


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
