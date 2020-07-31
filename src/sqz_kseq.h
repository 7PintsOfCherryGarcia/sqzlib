#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#include "sqz_data.h"


char sqz_getformat(const char *filename);


char sqz_kseqinit(sqzfastx_t *sqz);


size_t sqz_loadfastq(sqzfastx_t *sqz);


size_t sqz_newblock(sqzfastx_t *sqz);


size_t sqz_endblock(sqzfastx_t *sqz);


void sqz_kill(sqzfastx_t *sqz);
