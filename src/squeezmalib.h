#include <stdio.h>
#include <stdint.h>

#define LOAD_SIZE 16*1024

typedef struct  kseq_t kseq_t;

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
