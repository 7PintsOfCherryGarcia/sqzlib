#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#include "sqz_data.h"


unsigned char sqz_getformat(const char *filename);


char sqz_kseqinit(sqzfastx_t *sqz);


/*
  When a buffer is filled (Can't hold the entire kseq sequence + the
  lenth of the next sequence), the length of the current sequence is
  stored as well as any bases the buffer can accomodate.
*/
static uint64_t sqz_fastqnblock(sqzfastx_t *sqz, kseq_t *seq);


static uint64_t sqz_fastqeblock(sqzfastx_t *sqz);


static size_t sqz_fastanblock(sqzfastx_t *sqz, kseq_t *seq);


static uint64_t sqz_fastaeblock(sqzfastx_t *sqz);


void sqz_fastxkill(sqzfastx_t *sqz);


unsigned char sqz_checksqz(const char *buf);


sqzfastx_t *sqz_fastxinit(const char *filename, uint64_t bsize);


sqzblock_t *sqz_sqzblkinit(size_t size);


void sqz_sqzblkkill(sqzblock_t *blk);

uint64_t sqz_filesize(FILE *fp);


char sqz_readblksize(sqzblock_t *blk, FILE *fp);


size_t sqz_inflate(sqzblock_t *blk);


void sqz_decode(sqzfastx_t *sqz, sqzblock_t *blk, uint64_t klibl);


size_t sqz_fastqdecode(sqzblock_t *blk);


size_t sqz_fastXdecode(sqzblock_t *blk,
                       uint8_t *klibbuff,
                       uint64_t klibl,
                       char fastqflag);


static void sqz_rememberseq(sqzfastx_t *sqz, kseq_t *seq, char fqflag);