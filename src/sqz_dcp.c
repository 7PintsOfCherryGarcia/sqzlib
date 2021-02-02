#include <stdio.h>
#include <stdlib.h>

#define SQZLIB
#define KLIB
#include "sqz_data.h"
#include "sqz_dcp.h"


uint64_t sqz_filesize(FILE *fp)
{
    fseek(fp, 0, SEEK_END);
    long s = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    return s - 16;
}


sqzblock_t *sqz_sqzblkinit(size_t size)
{
    sqzblock_t *blk = malloc(sizeof(sqzblock_t));
    if (!blk) return NULL;
    blk->blkbuff = malloc(2*size);
    if (!blk->blkbuff) {
        fprintf(stderr, "[libsqz ERROR]: memory error.\n");
        free(blk);
        return NULL;
    }
    blk->blksize = 2*size;
    blk->blkpos  = 0;

    blk->namepos = 0;
    blk->newblk  = 1;
    //Compression buffer array
    blk->cmpbuff = malloc(2*size);
    if (!blk->cmpbuff) {
        fprintf(stderr, "[libsqz ERROR]: memory error.\n");
        free(blk->blkbuff);
        free(blk);
        return NULL;
    }
    blk->cmpsize = 2*size;
    blk->cmppos  = 0;
    return blk;
}



