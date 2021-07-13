#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include "sqz_data.h"

size_t sqz_deflate(sqzblock_t *blk, int level);


char sqz_zlibcmpdump(sqzblock_t *blk, size_t size, FILE *ofp);
