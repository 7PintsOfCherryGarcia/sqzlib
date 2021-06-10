#include <stdio.h>

#include "sqz_data.h"

#define TWO_BIT_MASK (3)

//Table to change "ACGT" to 0123 else to 4
unsigned char seq_nt4_tableSQZ[128] = {
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
unsigned char seq_dec_tableSQZ[128] = {
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


uint64_t sqz_seqencode(const uint8_t *seq, uint64_t lentocode, uint8_t *blkbuff);


uint64_t sqz_qualencode(const uint8_t *qual, uint8_t *blkbuff);


char sqz_fastaencode(sqzfastx_t *sqz, sqzblock_t *blk);




uint8_t sqz_8binqual(uint8_t q);


const uint8_t *sqz_findn(const uint8_t *seq);


uint64_t sqz_bit2encode(const uint8_t *seq, size_t seqlen);


sqzblock_t *sqz_sqzblkinit(size_t size);


char sqz_fastqheadblk(sqzfastx_t *sqz, sqzblock_t *blk);


char sqz_fastqtailblk(sqzfastx_t *sqz, sqzblock_t *blk);


char sqz_fastaheadblk(sqzfastx_t *sqz, sqzblock_t *blk);


char sqz_fastatailblk(sqzfastx_t *sqz, sqzblock_t *blk);


void sqz_blkdestroy(sqzblock_t *blk);


size_t sqz_qualdecode(const uint8_t *codebuff,
                      uint8_t *qualstr,
                      size_t length);


unsigned char sqz_bit2decode(const uint64_t *mer,
                             uint8_t *decoded,
                             unsigned char len);


unsigned char sqz_writens(unsigned char numn, uint8_t *decoded);



size_t sqz_fastadecode(const uint8_t *buff, size_t size);


size_t sqz_seqdecode(const uint8_t *codebuff,
                     uint8_t *decodebuff,
                     size_t length,
                     char qflag,
                     uint64_t *wbytes);


static uint64_t sqz_blkdecode(const uint8_t *codebuff,
                              uint8_t       *decodebuff,
                              uint64_t      *wbytes,
                              uint64_t      blklen);


static uint64_t sqz_nblkcode(uint8_t *buff, uint64_t nnum);


static uint64_t sqz_blkcode(uint64_t *buff, const uint8_t *seq, uint64_t len);
