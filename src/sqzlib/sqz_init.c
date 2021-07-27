#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "sqz_data.h"


sqzfastx_t *sqz_fastxinit(uint8_t fmt, uint64_t size)
{
    uint8_t ret = 1;
    sqzfastx_t *sqz = calloc(1, sizeof(sqzfastx_t));
    if (!sqz) return NULL;
    sqz->seqbuffer  = NULL;
    sqz->qualbuffer = NULL;
    sqz->namebuffer = NULL;
    sqz->readbuffer = NULL;
    sqz->pseq       = NULL;
    sqz->pqlt       = NULL;
    sqz->namesize   = NAME_SIZE;
    sqz->plen = 102400L;
    sqz->endthread = 128L;
    //Get file format if reading an sqz file (lower 3 bits of fmt)
    switch (fmt & 7) {
        case 0:
            goto exit;
        case 1:
            //FASTA
            sqz->seqbuffer = malloc(size + 1);
            if (!sqz->seqbuffer)  goto exit;
            sqz->seqbuffer[size] = '\0';
            sqz->pseq = malloc(102400L);
            if (!sqz->pseq)       goto exit;
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) goto exit;
            break;
        case 2:
            //FASTQ
            sqz->qualbuffer = malloc(size + 1);
            if (!sqz->qualbuffer) goto exit;
            sqz->qualbuffer[size] = 0;
            sqz->seqbuffer = malloc(size + 1);
            if (!sqz->seqbuffer)  goto exit;
            sqz->seqbuffer[size] = 0;
            sqz->pseq = malloc(102400L);
            sqz->pqlt = malloc(102400L);
            if (!sqz->pseq || !sqz->pqlt) goto exit;
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) goto exit;
            break;
        case 5:
            sqz->readbuffer = malloc(size + 1);
            if (!sqz->readbuffer) goto exit;
            sqz->readbuffer[size] = 0;
            sqz->seqbuffer = malloc(size + 1);
            if (!sqz->seqbuffer)  goto exit;
            sqz->seqbuffer[size] = 0;
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) goto exit;
            break;
        case 6:
            sqz->readbuffer = malloc(size + 1);
            if (!sqz->readbuffer) goto exit;
            sqz->readbuffer[size] = '\0';
            sqz->qualbuffer = malloc(size + 1);
            if (!sqz->qualbuffer) goto exit;
            sqz->qualbuffer[size] = 0;
            sqz->seqbuffer = malloc(size + 1);
            if (!sqz->seqbuffer)  goto exit;
            sqz->seqbuffer[size] = 0;
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) goto exit;
            break;
    }
    ret = 0;
    exit:
        if (ret) {
            free(sqz->seqbuffer),free(sqz->qualbuffer),free(sqz->namebuffer);
            free(sqz->readbuffer),free(sqz->pseq),free(sqz->pqlt);
            return NULL;
        }
        return sqz;
}


void sqz_fastxkill(sqzfastx_t *sqz)
{
    if (sqz) {
        free(sqz->seqbuffer);
        free(sqz->namebuffer);
        free(sqz->readbuffer);
        free(sqz->qualbuffer);
        free(sqz->readbuffer);
        free(sqz->pseq);
        free(sqz->pqlt);
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
