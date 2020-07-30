#include <stdio.h>

#define LOAD_SIZE 16*1024
#define NEND (128)

//Table to change "ACGT" to 0123 else to 4
unsigned char seq_nt4_table[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct kseq_t kseq_t;


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


typedef struct {
    uint8_t *codebuff;
    size_t  offset;
    char    newblk;
} sqzblock_t;


char sqz_encode(sqzfastx_t *sqz, sqzblock_t *blk);


size_t sqz_seqencode(const unsigned char *str,
                     uint32_t strlen,
                     sqzblock_t *codeblk);


size_t sqz_qualencode(const unsigned char *strqual, uint8_t *codebuff);


unsigned char sqz_8binqual(char q);


const unsigned char *sqz_findn(const unsigned char *strseq);


uint64_t sqz_bit2encode(const unsigned char *str, uint32_t strlen);


sqzblock_t *sqz_sqzblkinit(size_t size);


size_t sqz_headblk(sqzfastx_t *sqz, sqzblock_t *blk);


size_t sqz_tailblk(sqzfastx_t *sqz, sqzblock_t *blk);
