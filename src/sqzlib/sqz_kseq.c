#include <stdio.h>
#include <stdint.h>

#include "sqz_data.h"
#include "klib/kseq.h"
KSEQ_INIT(sqzFile, sqz_gzread)


static uint8_t sqz_checksqz(sqzFile sqzfp)
{
    int32_t magic = 0;
    uint8_t fmt   = 0;
    uint8_t sqz   = 0;
    //Read magic
    sqz_gzread(sqzfp, (void *)&magic, 4);
    if (MAGIC ^ magic) {
        //Not sqz file
        sqz_gzrewind(sqzfp);
        return 0;
    }
    //Set sqz flag
    fmt |= 4;
    //Read sequence format
    sqz_gzread(sqzfp, &sqz, 1);
    fmt |= sqz;
    //Read compression library
    sqz_gzread(sqzfp, &sqz, 1);
    fmt |= (sqz << 3);
    sqz_gzrewind(sqzfp);
    sqzfp->fmt = fmt;
    return 1;
}

static uint8_t sqz_loadname(sqzbuff_t *buff, kseq_t *seq)
{
    uint8_t ret = 0;
    uint64_t pos = buff->pos;
    uint8_t *data = (uint8_t *)buff->data;
    //Check namebuffer is large enough
    if ( (pos + seq->comment.l + seq->name.l) >= buff->size) {
        buff = sqz_buffrealloc(buff, buff->size<<1);
        if (!(buff))
            goto exit;
        data = (uint8_t *)buff->data;
    }
    memcpy(data + pos, seq->name.s, seq->name.l + 1);
    pos += seq->name.l + 1;
    //If comment exists
    if (seq->comment.s) {
        //Substitute terminating null with space
        data[pos - 1] = ' ';
        //Append comment including terminating null
        memcpy(data + pos, seq->comment.s, seq->comment.l + 1);
        pos += seq->comment.l + 1;
    }
    ret = 1;
    exit:
        buff->pos = pos;
        return ret;
}

static uint8_t sqz_cpyseq(sqzseq_t *seq, kseq_t *kseq)
{
    //TODO think better the realloc
    if (!kseq) return 1;
    if (kseq->seq.l > seq->l) {
        seq->s = realloc(seq->s, kseq->seq.l + 1);
        if (!seq->s) return 1;
        if (kseq->qual.s) {
            seq->q = realloc(seq->q, kseq->seq.l + 1);
            if (!seq->q) return 1;
        }
    }
    memcpy(seq->s, kseq->seq.s, kseq->seq.l + 1);
    if (kseq->qual.s)
        memcpy(seq->q, kseq->qual.s, kseq->seq.l + 1);
    seq->l = kseq->seq.l;
    if ( (kseq->name.l + kseq->comment.l) > seq->nlen ) {
        seq->n = realloc(seq->n, kseq->name.l + kseq->comment.l + 1);
        if (!seq->n) return 1;
    }
    seq->nlen = kseq->name.l + kseq->comment.l;
    memcpy(seq->n, kseq->name.s, kseq->name.l + 1);
    if (kseq->qual.s) {
        seq->n[kseq->name.l] = ' ';
        memcpy(seq->q, kseq->qual.s, kseq->qual.l + 1);
    }
    return 0;
}

static uint32_t sqz_fastanblock(sqzfastx_t *sqz, kseq_t *kseq)
{
    uint64_t offset = 0;
    uint64_t bases  = 0;
    uint64_t maxlen = sqz->size - B64 - 1;
    uint64_t l;
    uint32_t n      = 0;
    uint8_t *seq = sqz->seq;
    while ( (kseq_read(kseq) >= 0) ) {
        n++;
        l = kseq->seq.l;
        bases += l;
        if (!sqz_loadname(sqz->namebuffer, kseq)) {
            offset = 0;
            goto exit;
        }
        if (l > maxlen) {
            if ( l/3 > sqz->lseqbuff->size)
                if ( !(sqz->lseqbuff = sqz_buffrealloc(sqz->lseqbuff, l/3)) ) {
                    offset = 0;
                    n = 0;
                    goto exit;
                }
            sqz->lseqflag = 1;
            if ( sqz_cpyseq(sqz->lastseq, kseq) ) {
                offset = 0;
                goto exit;
            }
            goto exit;
        }
        memcpy(seq + offset, &l, B64);
        offset += B64;
        memcpy(seq + offset, kseq->seq.s, l + 1);
        offset += l + 1;
        if ( maxlen <= (l + 1 + B64) ) break;
        maxlen -= (l + 1 + B64);
    }
    exit:
        sqz->n += n;
        sqz->bases  = bases;
        sqz->offset = offset;
        return n;
}

static uint32_t sqz_fastqnblock(sqzfastx_t *sqz, kseq_t *kseq)
{
    uint64_t offset = 0;
    uint64_t bases  = 0;
    uint64_t maxlen = sqz->size - B64 - 1;
    uint64_t l;
    uint32_t n      = 0;
    uint8_t  *seq = sqz->seq;
    uint8_t  *qlt = sqz->qlt;
    while ( kseq_read(kseq) >= 0 ) {
        n++;
        l = kseq->seq.l;
        bases += l;
        if (!sqz_loadname(sqz->namebuffer, kseq)) {
            offset = 0;
            goto exit;
        }
        if ( l > maxlen ) {
            if (l/3 > sqz->lseqbuff->size)
                if ( !(sqz->lseqbuff = sqz_buffrealloc(sqz->lseqbuff, l/3)) ) {
                    offset = 0;
                    goto exit;
                }
            sqz->lseqflag = 1;
            if ( sqz_cpyseq(sqz->lastseq, kseq) ) {
                offset = 0;
                goto exit;
            }
            goto exit;
        }
        memcpy(seq + offset, &l, B64);
        offset += B64;
        memcpy(seq + offset, kseq->seq.s, l + 1);
        memcpy(qlt + offset, kseq->qual.s, l + 1);
        offset += (l + 1);
        if ( maxlen <= (l + 1 + B64) ) break;
        maxlen -= (l + 1 + B64);
    }
    exit:
        sqz->n += n;
        sqz->bases  = bases;
        sqz->offset = offset;
        return offset;
}

void sqz_getformat(sqzFile sqzfp)
{
    if ( (sqz_checksqz(sqzfp)) ) return;

    kseq_t *seq = kseq_init(sqzfp);
    if (!seq) goto exit;

    int l = kseq_read(seq);
    //ERROR
    if (l < 0)
        goto exit;
    //FASTQ
    if (seq->qual.l > 0) {
        sqzfp->fmt = 2;
        goto exit;
    }
    //FASTA
    sqzfp->fmt = 1;
    exit:
        kseq_destroy(seq);
        sqz_gzrewind(sqzfp);
        return;
}

uint32_t sqz_loadfastX(sqzfastx_t *sqz, uint8_t fqflag, kseq_t *seq)
{
    if (fqflag) return sqz_fastqnblock(sqz, seq);
    return sqz_fastanblock(sqz, seq);
}


