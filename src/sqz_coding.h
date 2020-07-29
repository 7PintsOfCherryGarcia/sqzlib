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


typedef struct {
    uint8_t *codebuff;
    size_t offset;
} sqzcodeblock_t;


int sqz_encode(sqzfastx_t *sqz, size_t sqzsize, sqzcodeblock_t *blk);


size_t sqz_seqencode(const unsigned char *str, uint32_t strlen, sqzcodeblock_t *codeblk);


size_t sqz_qualencode(const unsigned char *strqual, uint8_t *codebuff);


unsigned char sqz_8binqual(char q);


const unsigned char *sqz_findn(const unsigned char *strseq);


uint64_t bit2encode(const unsigned char *str, uint32_t strlen);
