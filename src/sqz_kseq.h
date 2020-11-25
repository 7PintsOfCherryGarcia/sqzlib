#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#include "sqz_data.h"


unsigned char sqz_getformat(const char *filename);


char sqz_kseqinit(sqzfastx_t *sqz);


uint64_t sqz_loadfastq(sqzfastx_t *sqz);


/*
  When a buffer is filled (Can't hold the entire kseq sequence + the
  lenth of the next sequence), the length of the current sequence is
  stored as well as any bases the buffer can accomodate.
*/
uint64_t sqz_fastqnblock(sqzfastx_t *sqz);


uint64_t sqz_fastqeblock(sqzfastx_t *sqz);


size_t sqz_loadfasta(sqzfastx_t *sqz);


size_t sqz_fastanblock(sqzfastx_t *sqz);


size_t sqz_fastaeblock(sqzfastx_t *sqz);


void sqz_kill(sqzfastx_t *sqz);


unsigned char sqz_checksqz(const char *buf);


uint64_t sqz_fastqwrap(sqzfastx_t *sqz, uint64_t offset);


uint64_t sqz_fastawrap(sqzfastx_t *sqz, uint64_t offset);

char sqz_loadname(sqzfastx_t *sqz, kstring_t name);
