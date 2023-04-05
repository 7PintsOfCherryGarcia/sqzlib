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

static uint8_t sqz_remeber(sqzfastx_t *sqz, kseq_t *seq, uint8_t fqflag)
{
    // TODO Pass sequence and length instead of kseq_t struct
    if (sqz->psize <= seq->seq.l) {
        if (fqflag) sqz->pqlt = realloc(sqz->pqlt, seq->seq.l*2);
        sqz->pseq = realloc(sqz->pseq, seq->seq.l*2);
        if ( ( fqflag && !(sqz->pqlt) ) || !(sqz->pseq)) return 1;
        sqz->psize = seq->seq.l*2;
    }
    memcpy(sqz->pseq, seq->seq.s, seq->seq.l + 1);
    if (fqflag) memcpy(sqz->pqlt, seq->qual.s, seq->seq.l + 1);
    return 0;
}

static uint8_t sqz_loadname(sqzfastx_t *sqz, kseq_t *seq)
{
    uint8_t ret = 0;
    uint8_t *namebuffer = sqz->namebuffer;
    uint64_t pos = sqz->namepos;
    //Check namebuffer is large enough
    if ( (pos + seq->comment.l + seq->name.l) >= sqz->namesize) {
        namebuffer = realloc(namebuffer, sqz->namesize*2);
        if (!(namebuffer))
            goto exit;
        sqz->namebuffer = namebuffer;
        sqz->namesize *= 2;
    }
    memcpy(namebuffer + pos, seq->name.s, seq->name.l + 1);
    pos += seq->name.l + 1;
    //If comment exists
    if (seq->comment.s) {
        //Substitute terminating null with space
        namebuffer[pos - 1] = ' ';
        //Append comment including terminating null
        memcpy(namebuffer + pos, seq->comment.s, seq->comment.l + 1);
        pos += seq->comment.l + 1;
    }
    ret = 1;
    exit:
        sqz->namepos = pos;
        return ret;
}

static uint32_t sqz_fastanblock(sqzfastx_t *sqz, kseq_t *kseq)
{
    uint64_t offset = 0;
    uint64_t l;
    uint64_t bases  = 0;
    uint64_t maxlen = sqz->size - B64 - 1;
    uint32_t n      = 0;
    uint8_t *seq = sqz->seq;
    uint64_t blksize = sqz->blk->blksize;
    while ( (kseq_read(kseq) >= 0) ) {
        l = kseq->seq.l;
        if ( l/4 > blksize) {
            sqz_sqzblkrealloc(sqz->blk, l/4);
            sqz->blk->blksize = l/4;
            blksize = l/4;
        }
        bases += l;
        n++;
        if (!sqz_loadname(sqz, kseq)) {
            offset = 0;
            goto exit;
        }
        if (l > maxlen) {
            fprintf(stderr, "EXIT from full buffer\n");
            sqz->endflag = 1;
            sqz->prevlen = l;
            sqz_remeber(sqz, kseq, 0);
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
        sqz->bases = bases;
        sqz->offset = offset;
        fprintf(stderr, "EXITING: %u\n", n);
        return n;
}

static uint64_t sqz_fastqnblock(sqzfastx_t *sqz, kseq_t *kseq)
{
    uint64_t offset = 0;
    uint64_t n      = 0;
    uint64_t l;
    uint64_t bases  = 0;
    uint64_t maxlen = sqz->size - B64 - 1;
    uint8_t  *seq = sqz->seq;
    uint8_t  *qlt = sqz->qlt;
    while ( kseq_read(kseq) >= 0 ) {
        l = kseq->seq.l;
        bases += l;
        n++;
        if (!sqz_loadname(sqz, kseq)) {
            offset = 0;
            goto exit;
        }
        if ( l > maxlen ) {
            sqz->endflag = 1;
            sqz->prevlen = l;
            sqz_remeber(sqz, kseq, 1);
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
        sqz->bases = bases;
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

uint64_t sqz_loadfastX(sqzfastx_t *sqz, uint8_t fqflag, kseq_t *seq)
{
    if (fqflag) return sqz_fastqnblock(sqz, seq);
    return sqz_fastanblock(sqz, seq);
}
