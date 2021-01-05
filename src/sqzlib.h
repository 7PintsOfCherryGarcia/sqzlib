/** \file sqzlib.h
*/

#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include <stdlib.h>
#include <zlib.h>

#define LOAD_SIZE 8*1024*1024
#define B64       sizeof(uint64_t)
#define HEADLEN   8

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
    Load fasta data into sqz struct
*/
size_t sqz_loadfasta(sqzfastx_t *sqz);


/*
    Encode fasta data already loaded in sqz struct and store it in blk struct
*/
char sqz_fastaencode(sqzfastx_t *sqz, sqzblock_t *blk);

/*
    Compress data already loaded in blk struct
*/
size_t sqz_deflate(sqzblock_t *blk, int level);


/*
    Load fastq data into sqz struct
*/
uint64_t sqz_loadfastq(sqzfastx_t *sqz);


/*
    Load fastq data already loaded in sqz struct and store it in blk struct
*/
char sqz_fastqencode(sqzfastx_t *sqz, sqzblock_t *blk);


/*
    Write compressed data stored in blk to file
*/
char sqz_zlibcmpdump(sqzblock_t *blk, size_t size, FILE *ofp);


/*
    Write sqz file tail
*/
char sqz_filetail(size_t numseqs, FILE *ofp);

/*
    Free blk object
*/
void sqz_blkdestroy(sqzblock_t *blk);


/*
    Decode data stored in buff. For fastq data streams
*/
size_t sqz_fastqdecode(const uint8_t *buff, size_t size);


/*
    Decode data stored in buff. For fasta data streams
*/
size_t sqz_fastadecode(const uint8_t *buff, size_t size);


/*
    Decompress zlib data stored in blk struct
*/
size_t sqz_inflate(sqzblock_t *blk);


/*
#########################################PENE###################################
KSEQ compatibility functions

In order to use kseq.h and take advantage of the quick fastX parsing, a
"reading" function must be provided. This "reading" function loads uncompressed
data into a buffer which is subsequently parsed by the rest of kseq.h functions.

The reading function must be of the form:
see gzread(): https://github.com/madler/zlib/blob/master/zlib.h
       int reader(FILE *file, void *buff, unsigned len)
Where
    file corresponds to a filepointer
    buff is a non empty buffer where compressed data is extracted
    len number of uncompressed bytes that should be read

This function returns the number of uncompressed bytes read, less than len for
end of file, or -1 for error.

################################################################################
*/
typedef struct {
    FILE *fp;
    sqzfastx_t *sqz;
    sqzblock_t *blk;
} sqz_File;


sqz_File sqz_sqzopen(char *filename);

//! A function
/**
 * Read uncompressed data in sqz file into a buffer.
 *
 * @param file - sqz_File pointer previously opened with sqz_sqzopen
 */
size_t sqz_sqzread(sqz_File *file, void *buff, size_t len);
