#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#include "sqz_data.h"



sqzfastx_t *sqz_fastxinit(const char *filename, size_t buffersize);


char sqz_getformat(const char *filename);


char sqz_kseqinit(sqzfastx_t *sqz);
