#include <stdio.h>
#include <stdint.h>

#include "sqz_data.h"
int64_t sqzread(sqzFile file, void *buff, uint64_t len);
void sqzclose(sqzFile file);
sqzFile sqzopen(const char *filename, const char *mode);
uint8_t sqz_checksqz(const char *filename);

#include "klib/kseq.h"
KSEQ_INIT(sqzFile, sqzread)


uint8_t sqz_getformat(const char *filename)
{
    uint8_t ret = 0;
    if ( (ret = sqz_checksqz(filename)) ) return ret;
    sqzFile fp = sqzopen(filename, "r");
    if (!fp) return ret;
    kseq_t *seq = kseq_init(fp);
    if (!seq) {
        sqzclose(fp);
        return ret;
    }
    int l = kseq_read(seq);
    //ERROR
    if (l < 0)
        goto exit;
    //FASTQ
    if (seq->qual.l > 0) {
        ret = 2;
        goto exit;
    }
    //FASTA
    ret = 1;
    exit:
        sqzclose(fp);
        kseq_destroy(seq);
        return ret;
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


static uint64_t sqz_fastqeblock(sqzfastx_t *sqz)
{
    kseq_t *pseq = sqz->pseq;
    uint64_t l = sqz->prevlen;
    uint64_t seqleft = l - sqz->seqread;
    char *seq = pseq->seq.s;
    char *qlt = pseq->qual.s;
    //buffer can be completely filled with current sequence
    if (seqleft >= LOAD_SIZE) {
        memcpy(sqz->seq, seq + sqz->seqread, LOAD_SIZE);
        memcpy(sqz->qlt, qlt + sqz->seqread, LOAD_SIZE);
        sqz->seqread += LOAD_SIZE;
        sqz->offset   = LOAD_SIZE;
        sqz->bases   += LOAD_SIZE;
        return LOAD_SIZE;
    }
    //Rest of sequence can go into buffer
    memcpy(sqz->seq, seq + sqz->seqread, seqleft);
    memcpy(sqz->qlt, qlt + sqz->seqread, seqleft);
    sqz->seq[seqleft] = 0;
    sqz->qlt[seqleft] = 0;
    seqleft++;
    sqz->endflag = 0;
    sqz->offset = seqleft;
    sqz->bases += seqleft;
    return seqleft;
}


static uint64_t sqz_fastaeblock(sqzfastx_t *sqz)
{
    kseq_t *pseq = sqz->pseq;
    uint64_t l = sqz->prevlen;
    uint64_t seqleft = l - sqz->seqread;
    char *seq = pseq->seq.s;
    //buffer can be completely filled with current sequence
    if (seqleft > LOAD_SIZE) {
        memcpy(sqz->seq, seq + sqz->seqread, LOAD_SIZE);
        sqz->seqread += LOAD_SIZE;
        sqz->offset   = LOAD_SIZE;
        sqz->bases   += LOAD_SIZE;
        return LOAD_SIZE;
    }
    //Rest of sequence can go into buffer
    memcpy(sqz->seq, seq + sqz->seqread, seqleft);
    sqz->seq[seqleft++] = '\0';
    sqz->endflag = 0;
    sqz->offset = seqleft;
    sqz->bases += seqleft;
    return seqleft;
}


static uint64_t sqz_fastqnblock(sqzfastx_t *sqz, kseq_t *kseq)
{
    uint64_t offset = 0;
    uint64_t n      = 0;
    uint64_t l;
    uint64_t bases  = 0;
    uint64_t maxlen = LOAD_SIZE - B64 - 1;
    uint8_t  *seq = sqz->seq;
    uint8_t  *qlt = sqz->qlt;
    while ( kseq_read(kseq) >= 0 ) {
        l = kseq->seq.l;
        n++;
        if (!sqz_loadname(sqz, kseq)) {
            offset = 0;
            goto exit;
        }
        if ( l > maxlen ) {
            memcpy(seq + offset, &l, B64);
            offset += B64;
            memcpy(seq + offset, kseq->seq.s, maxlen);
            memcpy(qlt + offset, kseq->qual.s, maxlen);
            offset += maxlen;
            seq[offset] = 0;
            qlt[offset] = 0;
            offset++;
            bases += maxlen;
            sqz->endflag = 1;
            sqz->seqread = maxlen;
            sqz->prevlen = l;
            sqz->pseq = kseq;
            goto exit;
        }
        bases += l;
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


static uint64_t sqz_fastanblock(sqzfastx_t *sqz, kseq_t *kseq)
{
    uint64_t offset = 0;
    uint64_t n      = 0;
    uint64_t l;
    uint64_t bases  = 0;
    uint64_t maxlen = LOAD_SIZE - B64 - 1;
    uint8_t *seq = sqz->seq;
    while ( (kseq_read(kseq) >= 0) ) {
        l = kseq->seq.l;
        n++;
        if (!sqz_loadname(sqz, kseq)) {
            offset = 0;
            goto exit;
        }
        if (l > maxlen) {
            memcpy(seq + offset, &l, B64);
            offset += B64;
            memcpy(seq + offset, kseq->seq.s, maxlen);
            offset += maxlen;
            seq[offset++] = 0;
            bases += maxlen;
            sqz->endflag = 1;
            sqz->seqread = maxlen;
            sqz->prevlen = l;
            sqz->pseq = kseq;
            goto exit;
        }
        bases += l;
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
        return offset;
}


uint64_t sqz_loadfastX(sqzfastx_t *sqz, uint8_t fqflag, kseq_t *seq)
{
    if (sqz->endflag) {
        if (fqflag) return sqz_fastqeblock(sqz);
        return sqz_fastaeblock(sqz);
    }
    if (fqflag) return sqz_fastqnblock(sqz, seq);
    return sqz_fastanblock(sqz, seq);
}
