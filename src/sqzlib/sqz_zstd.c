#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zstd.h>
#include "sqz_data.h"

int64_t sqz_zstdcompress(sqzblock_t *blk, int level)
{
    int64_t ret = 0;

    ret = ZSTD_compress(blk->cmpbuff->data, blk->cmpbuff->size,
                        blk->blkbuff->data, blk->blkbuff->pos, level);
    if (ZSTD_isError(ret)) return -1;

    return ret;
}


uint64_t sqz_zstddecompress(sqzblock_t *blk)
{
    int64_t ret = 0;

    ret = ZSTD_decompress(blk->blkbuff->data, blk->blkbuff->size,
                          blk->cmpbuff->data, blk->cmpbuff->size);
    if (ZSTD_isError(ret)) return 0;
    return (uint64_t)ret;
}
