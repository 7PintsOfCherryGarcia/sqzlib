#include <stdio.h>

#include "sqz_data.h"


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


char sqz_encode(sqzfastx_t *sqz, sqzblock_t *blk);


void sqz_seqencode(const uint8_t *seq, size_t seqlen, sqzblock_t *blk);


size_t sqz_qualencode(const unsigned char *strqual, uint8_t *codebuff);


unsigned char sqz_8binqual(char q);


const uint8_t *sqz_findn(const uint8_t *seq);


uint64_t sqz_bit2encode(const uint8_t *seq, size_t seqlen);


sqzblock_t *sqz_sqzblkinit(size_t size);


char sqz_headblk(sqzfastx_t *sqz, sqzblock_t *blk);


char sqz_tailblk(sqzfastx_t *sqz, sqzblock_t *blk);


void sqz_blkdestroy(sqzblock_t *blk);
