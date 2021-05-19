#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#define SQZLIB
#define KLIB
#include "sqz_init.h"


sqzfastx_t *sqz_fastxinit(const char *filename, uint64_t bsize)
{
    sqzfastx_t *sqz = malloc(sizeof(sqzfastx_t));
    if (!sqz) return NULL;
    {
    sqz->filename   = filename;
    sqz->fp         = NULL;
    sqz->fmt        = 0;
    sqz->endflag    = 0;
    sqz->cmpflag    = 0;
    sqz->offset     = 0;
    sqz->seq        = NULL;
    sqz->seqbuffer  = NULL;
    sqz->qualbuffer = NULL;
    sqz->namebuffer = NULL;
    sqz->readbuffer = NULL;
    sqz->namepos    = 0;
    sqz->namesize   = NAME_SIZE;
    sqz->n          = 0;
    sqz->bases      = 0;
    sqz->seqread    = 0;
    sqz->rem        = 0;
    sqz->toread     = 0;
    sqz->prevlen    = 0;
    }
    //Get file format if reading an sqz file
    unsigned char fmt = sqz_getformat(filename);

    //Initialize kseq objects
    if ( fmt <= 2 ) {
        if (!sqz_kseqinit(sqz)) {
            free(sqz);
            return NULL;
        }
    }
    switch (fmt & 7) {
        case 0:
            free(sqz);
            sqz = NULL;
            break;
        case 1:
            //FASTA
            sqz->fmt = fmt;
            sqz->qualbuffer = NULL;
            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                free(sqz);
                sqz = NULL;
                break;
            }
            sqz->seqbuffer[bsize] = '\0';
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) {
                free(sqz->seqbuffer);
                free(sqz);
                sqz = NULL;
            }
            break;
        case 2:
            //FASTQ
            sqz->fmt = fmt;
            sqz->qualbuffer = malloc(bsize + 1);
            if (!sqz->qualbuffer) {
                free(sqz);
                sqz = NULL;
                break;
            }
            sqz->qualbuffer[bsize] = 0;
            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                free(sqz->qualbuffer);
                free(sqz);
                sqz = NULL;
                break;
            }
            sqz->seqbuffer[bsize] = 0;
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) {
                free(sqz->qualbuffer);
                free(sqz->seqbuffer);
                free(sqz);
                sqz = NULL;
            }
            break;
        case 5:
            sqz->fmt = fmt;

            sqz->readbuffer = malloc(bsize + 1);
            if (!sqz->readbuffer) {
                free(sqz);
                sqz = NULL;
                break;
            }
            sqz->readbuffer[bsize] = 0;

            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                free(sqz->readbuffer);
                free(sqz);
                sqz = NULL;
                break;
            }
            sqz->seqbuffer[bsize] = 0;

            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) {
                free(sqz->readbuffer);
                free(sqz->seqbuffer);
                free(sqz);
                sqz = NULL;
            }
            break;
        case 6:
            sqz->fmt = fmt;

            sqz->readbuffer = malloc(bsize + 1);
            if (!sqz->readbuffer) {
                free(sqz);
                sqz = NULL;
                break;
            }
            sqz->readbuffer[bsize] = '\0';

            sqz->qualbuffer = malloc(bsize + 1);
            if (!sqz->qualbuffer) {
                free(sqz->readbuffer);
                free(sqz);
                sqz = NULL;
                break;
            }
            sqz->qualbuffer[bsize] = 0;

            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                free(sqz->readbuffer);
                free(sqz->qualbuffer);
                free(sqz);
                sqz = NULL;
                break;
            }
            sqz->seqbuffer[bsize] = 0;

            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) {
                free(sqz->readbuffer);
                free(sqz->qualbuffer);
                free(sqz->seqbuffer);
                free(sqz);
                sqz = NULL;
            }
            break;
    }
    return sqz;
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
