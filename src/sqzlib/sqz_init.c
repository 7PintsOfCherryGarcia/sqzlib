#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "sqz_data.h"


sqzfastx_t *sqz_fastxinit(uint8_t fmt, uint64_t size)
{
    uint8_t ret = 1;
    sqzfastx_t *sqz = calloc(1, sizeof(sqzfastx_t));
    if (!sqz) return NULL;
    sqz->namesize   = NAME_SIZE;
    sqz->endthread  = 128L;

    sqz->seq = malloc(size + 1);
    if (!sqz->seq)  goto exit;
    sqz->seq[size] = '\0';
    sqz->namebuffer = malloc(NAME_SIZE);
    if (!sqz->namebuffer) goto exit;
    sqz->pseq = malloc(16384);
    if (!sqz->pseq) goto exit;
    //Get file format if reading an sqz file (lower 3 bits of fmt)
    switch (fmt & 7) {
        case 0:
            goto exit;
        case 1:
            //FASTA do nothing, everything allocated
            break;
        case 2:
            //FASTQ
            sqz->qlt = malloc(size + 1);
            if (!sqz->qlt) goto exit;
            sqz->qlt[size] = 0;
            sqz->pqlt = malloc(16384);
            if (!sqz->pqlt) goto exit;
            break;
        case 5:
            sqz->readbuffer = malloc(size + 1);
            if (!sqz->readbuffer) goto exit;
            sqz->readbuffer[size] = 0;
            break;
        case 6:
            sqz->readbuffer = malloc(size + 1);
            if (!sqz->readbuffer) goto exit;
            sqz->readbuffer[size] = '\0';
            break;
    }
    ret = 0;
    exit:
        if (ret) {
            free(sqz->seq),free(sqz->qlt),free(sqz->namebuffer);
            free(sqz->readbuffer);
            return NULL;
        }
        return sqz;
}


void sqz_fastxkill(sqzfastx_t *sqz)
{
    if (sqz) {
        free(sqz->seq);
        free(sqz->namebuffer);
        free(sqz->readbuffer);
        free(sqz->qlt);
        free(sqz);
    }
}


sqzblock_t *sqz_sqzblkinit(uint64_t size)
{
    sqzblock_t *blk = malloc(sizeof(sqzblock_t));
    if (!blk) return NULL;
    //Encoding buffer
    blk->blkbuff = malloc(2*size);
    if (!blk->blkbuff) {
        free(blk);
        return NULL;
    }
    blk->blksize = 2*size;
    blk->mblksize = 2*size;
    blk->blkpos  = 0;

    blk->namepos = 0;
    blk->prevlen = 0;
    blk->newblk  = 1;

    //Compression buffer
    blk->cmpbuff = malloc(2*size);
    if (!blk->cmpbuff) {
        free(blk->blkbuff);
        free(blk);
        return NULL;
    }
    blk->cmpsize = 2*size;
    blk->mcmpsize = 2*size;
    blk->cmppos  = 0;
    return blk;
}


void sqz_blkkill(sqzblock_t *blk)
{
    if (blk) {
        free(blk->blkbuff);
        free(blk->cmpbuff);
        free(blk);
    }
}
