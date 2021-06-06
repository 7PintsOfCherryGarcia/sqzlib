#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#include "sqz_data.h"



char sqz_fastxinit(sqzfastx_t *sqz, unsigned char fmt, uint64_t bsize);


char sqz_getformat(const char *filename);


char sqz_kseqinit(sqzfastx_t *sqz);
