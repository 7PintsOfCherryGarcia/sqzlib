#include <stdio.h>

#include "sqz_data.h"

#define TWO_BIT_MASK (3)

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


//Table to change 01234 to ACGTN
unsigned char seq_dec_table[128] = {
    'A', 'C','G', 'T',  'N', 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


unsigned char qual_val_table[8] = {33,39,48,55,60,66,70,73};


char sqz_encode(sqzfastx_t *sqz, sqzblock_t *blk);


void sqz_seqencode(const uint8_t *seq,
                   uint64_t readlen,
                   sqzblock_t *blk,
                   uint64_t seqlen);


size_t sqz_qualencode(const uint8_t *qual,
                      size_t readlen,
                      sqzblock_t *blk,
                      uint64_t seqlen);


uint8_t sqz_8binqual(uint8_t q);


const uint8_t *sqz_findn(const uint8_t *seq);


uint64_t sqz_bit2encode(const uint8_t *seq, size_t seqlen);


sqzblock_t *sqz_sqzblkinit(size_t size);


char sqz_headblk(sqzfastx_t *sqz, sqzblock_t *blk);


char sqz_tailblk(sqzfastx_t *sqz, sqzblock_t *blk);


void sqz_blkdestroy(sqzblock_t *blk);


size_t sqz_qualdecode(const uint8_t *codebuff, char *qualstr, size_t length);
//size_t sqz_qualdecode(const uint8_t *buff, char *uncode, size_t length);


unsigned char sqz_bit2decode(const uint64_t *mer, char *decoded, unsigned char len);


unsigned char sqz_writens(unsigned char numn, char *decoded);


size_t sqz_seqdecode(const uint8_t *codebuff,
                     char *seqstr,
                     char *qualstr,
                     size_t length);
//size_t sqz_loopdecode(size_t length, const uint8_t *seqbuffer, char *seqstr);
