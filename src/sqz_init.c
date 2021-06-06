#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#define SQZLIB
#define KLIB
#include "sqz_init.h"


char sqz_fastxinit(sqzfastx_t *sqz, unsigned char fmt, uint64_t bsize)
{
    char ret = 0;
    sqz->namesize   = NAME_SIZE;
    //Get file format if reading an sqz file
    switch (fmt & 7) {
        case 0:
            goto exit;
        case 1:
            //FASTA
            sqz->qualbuffer = NULL;
            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) goto exit;
            sqz->seqbuffer[bsize] = '\0';
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) {
                free(sqz->seqbuffer);
                sqz->seqbuffer = NULL;
            }
            break;
        case 2:
            //FASTQ
            sqz->qualbuffer = malloc(bsize + 1);
            if (!sqz->qualbuffer) goto exit;
            sqz->qualbuffer[bsize] = 0;
            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                free(sqz->qualbuffer);
                sqz->qualbuffer = NULL;
                goto exit;
            }
            sqz->seqbuffer[bsize] = 0;
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) {
                free(sqz->qualbuffer);
                free(sqz->seqbuffer);
                sqz->qualbuffer = NULL, sqz->seqbuffer  = NULL;
                goto exit;
            }
            break;
        case 5:
            sqz->readbuffer = malloc(bsize + 1);
            if (!sqz->readbuffer) goto exit;
            sqz->readbuffer[bsize] = 0;
            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                free(sqz->readbuffer);
                sqz->readbuffer = NULL;
                goto exit;
            }
            sqz->seqbuffer[bsize] = 0;
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) {
                free(sqz->readbuffer);
                free(sqz->seqbuffer);
                sqz->readbuffer = NULL, sqz->seqbuffer = NULL;
                goto exit;
            }
            break;
        case 6:
            sqz->readbuffer = malloc(bsize + 1);
            if (!sqz->readbuffer) goto exit;
            sqz->readbuffer[bsize] = '\0';
            sqz->qualbuffer = malloc(bsize + 1);
            if (!sqz->qualbuffer) {
                free(sqz->readbuffer);
                sqz->readbuffer = NULL;
                goto exit;
            }
            sqz->qualbuffer[bsize] = 0;
            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                free(sqz->readbuffer);
                free(sqz->qualbuffer);
                sqz->readbuffer = NULL, sqz->qualbuffer = NULL;
                goto exit;
            }
            sqz->seqbuffer[bsize] = 0;
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) {
                free(sqz->readbuffer);
                free(sqz->qualbuffer);
                free(sqz->seqbuffer);
                sqz->readbuffer=NULL,sqz->qualbuffer=NULL;sqz->seqbuffer=NULL;
                goto exit;
            }
            break;
    }
    ret = 1;
    exit:
        return ret;
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
    blk->cmppos  = 0;

    return blk;
}
