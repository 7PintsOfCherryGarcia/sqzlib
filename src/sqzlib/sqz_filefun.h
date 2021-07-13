#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "sqz_data.h"

uint64_t chachasmooth(sqz_File *fp, uint8_t *buff, uint64_t start, uint64_t size);

char sqz_readblksize(sqzblock_t *blk, FILE *fp);


char sqz_getformat(const char *filename);


char sqz_kseqinit(sqzfastx_t *sqz);


uint64_t sqz_filesize(FILE *fp);


static void sqz_fastxreset(sqzfastx_t *sqz);


static void sqz_blkreset(sqzblock_t *blk);
