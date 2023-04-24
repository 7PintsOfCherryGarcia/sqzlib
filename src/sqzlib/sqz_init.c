#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "sqz_data.h"

static sqzbuff_t *sqz_buffinit(uint64_t size)
{
    sqzbuff_t *buff = calloc(1, sizeof(sqzbuff_t));
    if (!buff) return NULL;
    buff->data = calloc(size, sizeof(uint8_t));
    if (!buff->data) {
        free(buff);
        return NULL;
    }
    buff->size = size;
    buff->pos  = 0;
    return buff;
}

static void sqz_buffkill(sqzbuff_t *buff)
{
    if (buff) {
        if (buff->data)
            free(buff->data);
        free(buff);
    }
}

static sqzseq_t *sqz_seqinit(void)
{
    sqzseq_t *seq = calloc(1, sizeof(sqzseq_t));
    if (!seq) return NULL;
    seq->n = calloc(256, 1);
    if (!seq->n) {
        free(seq);
        return NULL;
    }
    seq->nlen = 0;
    return seq;
}

static void sqz_seqkill(sqzseq_t *seq)
{
    if (seq) {
        if (seq->n) free(seq->n);
        if (seq->s) free(seq->s);
        if (seq->q) free(seq->q);
        free(seq);
    }
}

static sqzblock_t *sqz_sqzblkinit(uint64_t size)
{
    sqzblock_t *blk = calloc(1, sizeof(sqzblock_t));
    if (!blk) return NULL;
    blk->blkbuff = sqz_buffinit(size);
    if (!blk->blkbuff) {
        free(blk);
        return NULL;
    }
    blk->cmpbuff = sqz_buffinit(size);
    if (!blk->cmpbuff) {
        free(blk->blkbuff);
        free(blk);
        return NULL;
    }
    return blk;
}

static void sqz_blkkill(sqzblock_t *blk)
{
    if (blk) {
        sqz_buffkill(blk->blkbuff);
        sqz_buffkill(blk->cmpbuff);
        free(blk);
    }
}

sqzseq_t *sqz_seqrealloc(sqzseq_t *seq, uint64_t newsize)
{
    seq->s = realloc(seq->s, newsize + 1);
    seq->q = realloc(seq->q, newsize + 1);
    if ( !seq->s || !seq->q ) {
        free(seq);
        return NULL;
    }
    seq->maxl = newsize;
    return seq;
}

sqzbuff_t *sqz_buffrealloc(sqzbuff_t *buff,  uint64_t size)
{
    buff->data = realloc(buff->data, size);
    if (!buff->data) return NULL;
    buff->size = size;
    if (buff->pos > size) buff->pos = size - 1;
    return buff;
}

sqzfastx_t *sqz_fastxinit(uint8_t fmt, uint64_t size)
{
    uint8_t ret = 1;
    sqzfastx_t *sqz = calloc(1, sizeof(sqzfastx_t));
    if (!sqz) return NULL;
    sqz->datflag    = 128U;
    sqz->size       = size;
    sqz->seq = calloc(size + 1, sizeof(uint8_t));
    sqz->namebuffer = sqz_buffinit(NAME_SIZE);
    sqz->lastseq    = sqz_seqinit();
    if (!sqz->namebuffer || !sqz->seq || !sqz->lastseq) goto exit;
    if ( !(sqz->blk = sqz_sqzblkinit(size)) )           goto exit;
    if (!(sqz->lseqbuff = sqz_buffinit(16384UL)))       goto exit;
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
            break;
        case 5:
            sqz->readbuffer = sqz_buffinit(size + 1);
            if (!sqz->readbuffer) goto exit;
            break;
        case 6:
            sqz->readbuffer = sqz_buffinit(size + 1);
            if (!sqz->readbuffer) goto exit;
            break;
    }
    ret = 0;
    exit:
        if (ret) {
            free(sqz->seq);
            free(sqz->qlt);
            free(sqz->namebuffer);
            free(sqz->readbuffer);
            free(sqz);
            return NULL;
        }
        return sqz;
}

void sqz_fastxkill(sqzfastx_t *sqz)
{
    if (sqz) {
        if (sqz->seq)
            free(sqz->seq);
        if (sqz->qlt)
            free(sqz->qlt);
        sqz_buffkill(sqz->namebuffer);
        sqz_buffkill(sqz->readbuffer);
        sqz_buffkill(sqz->lseqbuff);
        sqz_blkkill(sqz->blk);
        sqz_seqkill(sqz->lastseq);
        free(sqz);
    }
}

uint8_t sqz_blkrealloc(sqzblock_t *blk, uint64_t newsize)
{
    blk->blkbuff = sqz_buffrealloc(blk->blkbuff, newsize);
    if ( !(blk->blkbuff) ) return 1;
    blk->cmpbuff = sqz_buffrealloc(blk->cmpbuff, newsize);
    if ( !(blk->cmpbuff) ) return 1;
    return 0;
}



uint64_t sqz_seqsinblk(sqzblock_t *blk)
{
    return blk->blkbuff->pos;
    //uint64_t s = *(uint64_t *)( blk->blkbuff + ( blk->blksize - B64 ) );
    //uint64_t d = blk->blksize - B64 - s;
    //char *names = (char *)(blk->blkbuff + d);
    //uint64_t p = 0;
    //uint64_t n = 0;
    //while (p != s) {
    //    n++;
    //    p += strlen(names + p) + 1;
    //}
    //return n;
}
