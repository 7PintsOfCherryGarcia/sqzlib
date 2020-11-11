#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#include "sqz_data.h"


unsigned char sqz_getformat(const char *filename);


char sqz_kseqinit(sqzfastx_t *sqz);


size_t sqz_loadfastq(sqzfastx_t *sqz);


size_t sqz_newblock(sqzfastx_t *sqz);


size_t sqz_endblock(sqzfastx_t *sqz);


size_t sqz_loadfasta(sqzfastx_t *sqz);


size_t sqz_fastanblock(sqzfastx_t *sqz);


size_t sqz_fastaeblock(sqzfastx_t *sqz);


void sqz_kill(sqzfastx_t *sqz);


unsigned char sqz_checksqz(const char *buf);
