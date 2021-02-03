#include <stdio.h>

#include "sqz_data.h"

#define TWO_BIT_MASK (3)

//Table to change "ACGT" to 0123 else to 4
unsigned char seq_nt4_table[128] = {
    128, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
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


char sqz_fastqencode(sqzfastx_t *sqz, sqzblock_t *blk);


uint64_t sqz_seqencode(const uint8_t *seq,
                       uint64_t seqlen,
                       uint8_t *blkbuff,
                       uint64_t seqlenOG);


uint64_t sqz_qualencode(const uint8_t *qual,
                        uint64_t quallen,
                        uint8_t *blkbuff,
                        uint64_t seqlen);


char sqz_fastaencode(sqzfastx_t *sqz, sqzblock_t *blk);




uint8_t sqz_8binqual(uint8_t q);


const uint8_t *sqz_findn(const uint8_t *seq);


uint64_t sqz_bit2encode(const uint8_t *seq, size_t seqlen);


sqzblock_t *sqz_sqzblkinit(size_t size);


char sqz_headblk(sqzfastx_t *sqz, sqzblock_t *blk);


char sqz_tailblk(sqzfastx_t *sqz, sqzblock_t *blk);


char sqz_fastaheadblk(sqzfastx_t *sqz, sqzblock_t *blk);


char sqz_fastatailblk(sqzfastx_t *sqz, sqzblock_t *blk);


void sqz_blkdestroy(sqzblock_t *blk);


size_t sqz_qualdecode(const uint8_t *codebuff,
                      uint8_t *qualstr,
                      size_t length,
                      uint64_t *wbytes);


unsigned char sqz_bit2decode(const uint64_t *mer,
                             uint8_t *decoded,
                             unsigned char len);


unsigned char sqz_writens(unsigned char numn, uint8_t *decoded);


size_t sqz_seqdecode(const uint8_t *codebuff,
                     char *seqstr,
                     char *qualstr,
                     size_t length,
                     char flag);


size_t sqz_fastadecode(const uint8_t *buff, size_t size);

//tmp
size_t sqz_seqdecode2(const uint8_t *codebuff,
                      uint8_t *decodebuff,
                      size_t length,
                      char qflag,
                      uint64_t *wbytes);
