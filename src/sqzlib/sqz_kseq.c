#include <stdio.h>
#include <stdint.h>

#include "sqz_data.h"
int64_t sqzread(sqzFile file, void *buff, uint64_t len);
void sqzclose(sqzFile file);
sqzFile sqzopen(const char *filename, const char *mode);

#include "klib/kseq.h"
KSEQ_INIT(sqzFile, sqzread)


static uint8_t sqz_checksqz(sqzFile sqzfp)
{
    uint64_t tmp   = 0;
    uint32_t magic = 0;
    uint8_t  fmt   = 0;
    uint8_t  sqz   = 0;
    //Read magic
    tmp += fread(&magic, 1, 4, sqzfp->fp);
    if (MAGIC ^ magic) {
        rewind(sqzfp->fp);
        return 0;
    }
    //Set sqz flag
    fmt |= 4;
    //Read sequence format
    tmp += fread(&sqz, 1, 1, sqzfp->fp);
    fmt |= sqz;
    //Read compression library
    tmp += fread(&sqz, 1, 1, sqzfp->fp);
    fmt |= (sqz << 3);
    rewind(sqzfp->fp);
    return fmt;
}


uint8_t sqz_getformat(sqzFile sqzfp)
{
    uint8_t ret = 0;
    if ( (ret = sqz_checksqz(sqzfp)) ) return ret;

    kseq_t *seq = kseq_init(sqzfp);
    if (!seq) return ret;
    int l = kseq_read(seq);
    //ERROR
    if (l < 0) goto exit;
    //FASTQ
    if (seq->qual.l > 0) {
        ret = 2;
        goto exit;
    }
    //FASTA
    ret = 1;
    exit:
        kseq_destroy(seq);
        if (ret) gzrewind(sqzfp->gzfp);
        return ret;
}


static uint8_t sqz_remeber(sqzfastx_t *sqz, kseq_t *seq, uint8_t fqflag)
{
    if (sqz->psize < seq->seq.l) {
        if (fqflag) sqz->pqlt = realloc(sqz->pqlt, seq->seq.l);
        sqz->pseq = realloc(sqz->pseq, seq->seq.l);
        if (!(sqz->pqlt) || !(sqz->pseq)) return 1;
        sqz->psize = seq->seq.l;
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
            sqz_remeber(sqz, kseq, 0);
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


static uint64_t sqz_fastaeblock(sqzfastx_t *sqz)
{
    uint64_t l = sqz->prevlen;
    uint64_t seqleft = l - sqz->seqread;
    uint8_t *seq = sqz->pseq;
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
            //TODO better errir
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
            sqz_remeber(sqz, kseq, 1);
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


static uint64_t sqz_fastqeblock(sqzfastx_t *sqz)
{
    uint64_t l = sqz->prevlen;
    uint64_t seqleft = l - sqz->seqread;
    uint8_t *seq = sqz->pseq;
    uint8_t *qlt = sqz->pqlt;
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
    sqz->qlt[seqleft++] = 0;
    sqz->endflag = 0;
    sqz->offset = seqleft;
    sqz->bases += seqleft;
    return seqleft;
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
