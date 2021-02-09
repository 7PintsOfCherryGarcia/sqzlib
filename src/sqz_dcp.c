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


void sqzrewind(sqz_File *sqzfp)
{
    sqz_fastxreset(sqzfp->sqz);
    sqz_blkreset(sqzfp->blk);
    fseek(sqzfp->fp, HEADLEN, SEEK_SET);
    sqzfp->filepos = ftell(sqzfp->fp);
    sqzfp->ff = 0;
}


static void sqz_fastxreset(sqzfastx_t *sqz)
{
    sqz->endflag    = 0;
    sqz->cmpflag    = 0;
    sqz->offset     = 0;
    sqz->namepos    = 0;
    sqz->n          = 0;
    sqz->bases      = 0;
    sqz->rem        = 0;
    sqz->toread     = 0;
    sqz->prevlen    = 0;
}


static void sqz_blkreset(sqzblock_t *blk)
{
    blk->blkpos  = 0;
    blk->namepos = 0;
    blk->newblk  = 1;
    blk->cmppos  = 0;
}
