#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "sqz_data.h"


static uint8_t sqz_checksqz(const char *filename)
{
    uint64_t tmp   = 0;
    uint32_t magic = 0;
    uint8_t  fmt   = 0;
    uint8_t  sqz   = 0;
    FILE *fp = fopen(filename, "rb");
    if (!fp) return 0;
    //Read magic
    tmp += fread(&magic, 1, 4, fp);
    if (MAGIC ^ magic) {
        //Return 0 if not an sqz file
        fclose(fp);
        return 0;
    }
    //Set sqz flag (bit 6 of fmt)
    fmt |= 4;
    //Read format: 1 - fastA or 2 - fastQ (bits 7 and 8 of fmt)
    tmp += fread(&sqz, 1, 1, fp);
    fmt |= sqz;
    //Read compression library: 1 - zlib (bits 1,2,3,4,5 of fmt)
    tmp += fread(&sqz, 1, 1, fp);
    //fprintf(stderr, "libfmt: %u\n", sqz);
    fmt |= (sqz << 3);
    fclose(fp);
    return fmt;
}


uint8_t  sqz_getformat(const char *filename)
{
    uint8_t ret = 0;
    if ( (ret = sqz_checksqz(filename)) ) return ret;
    ret = 0;
    gzFile fp = gzopen(filename, "r");
    if (!fp) return ret;

    kseq_t *seq = kseq_init(fp);
    if (!seq) {
        gzclose(fp);
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
        gzclose(fp);
        kseq_destroy(seq);
        return ret;
}


static void sqz_rememberseq(sqzfastx_t *sqz, kseq_t *seq, char fqflag)
{
    if (sqz->plen < seq->seq.l) {
        sqz->plen = seq->seq.l + 1024;
        if (fqflag) sqz->pqlt = realloc(sqz->pqlt, sqz->plen);
        sqz->pseq = realloc(sqz->pseq, sqz->plen);
    }
    if (fqflag) strcpy(sqz->pqlt, seq->qual.s);
    strcpy(sqz->pseq, seq->seq.s);
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
    uint64_t l = sqz->prevlen;
    uint64_t seqleft = l - sqz->seqread;
    char *seq = sqz->pseq;
    char *qlt = sqz->pqlt;
    //buffer can be completely filled with current sequence
    if (seqleft >= LOAD_SIZE) {
        memcpy(sqz->seqbuffer, seq + sqz->seqread, LOAD_SIZE);
        memcpy(sqz->qualbuffer, qlt + sqz->seqread, LOAD_SIZE);
        sqz->seqread += LOAD_SIZE;
        sqz->offset   = LOAD_SIZE;
        sqz->bases   += LOAD_SIZE;
        return LOAD_SIZE;
    }
    //Rest of sequence can go into buffer
    memcpy(sqz->seqbuffer, seq + sqz->seqread, seqleft);
    memcpy(sqz->qualbuffer, qlt + sqz->seqread, seqleft);
    sqz->seqbuffer[seqleft] = 0;
    sqz->qualbuffer[seqleft] = 0;
    seqleft++;
    sqz->endflag = 0;
    sqz->offset = seqleft;
    sqz->bases += seqleft;
    return seqleft;
}


static uint64_t sqz_fastaeblock(sqzfastx_t *sqz)
{
    uint64_t l = sqz->prevlen;
    uint64_t seqleft = l - sqz->seqread;
    char *seq = sqz->pseq;
    //buffer can be completely filled with current sequence
    if (seqleft > LOAD_SIZE) {
        memcpy(sqz->seqbuffer, seq + sqz->seqread, LOAD_SIZE);
        sqz->seqread += LOAD_SIZE;
        sqz->offset   = LOAD_SIZE;
        sqz->bases   += LOAD_SIZE;
        return LOAD_SIZE;
    }
    //Rest of sequence can go into buffer
    memcpy(sqz->seqbuffer, seq + sqz->seqread, seqleft);
    sqz->seqbuffer[seqleft++] = '\0';
    sqz->endflag = 0;
    sqz->offset = seqleft;
    sqz->bases += seqleft;
    return seqleft;
}


static uint64_t sqz_fastqnblock(sqzfastx_t *sqz, kseq_t *seq)
{
    uint64_t offset = 0;
    uint64_t n      = 0;
    uint64_t l;
    uint64_t bases  = 0;
    uint64_t maxlen = LOAD_SIZE - B64;
    uint8_t  *seqbuffer = sqz->seqbuffer;
    uint8_t  *qltbuffer = sqz->qualbuffer;
    while ( kseq_read(seq) >= 0 ) {
        //fprintf(stderr, "%s\n", seq->name.s);
        l = seq->seq.l;
        n++;
        if (!sqz_loadname(sqz, seq)) {
            offset = 0;
            goto exit;
        }
        //Determine if current sequence can be loaded completely
        if ( l + 1  > maxlen) {
            memcpy(seqbuffer + offset, &l, B64);
            offset += B64;
            memcpy(seqbuffer + offset, seq->seq.s, maxlen);
            memcpy(qltbuffer + offset, seq->qual.s, maxlen);
            offset += maxlen;
            //TODO this is not strictly necessary, but it brings me peaceof mind
            seqbuffer[offset] = 0;
            qltbuffer[offset] = 0;
            //From previous TODO
            //offset should be LOAD_SIZE + 1 after this increment
            offset++;
            bases += maxlen;
            sqz->endflag = 1;
            sqz->seqread = maxlen;
            sqz->prevlen = l;
            /*
              The sqzfastx_t struct must remember the sequence (and qualities)
              of this partially decoded record. Thus the contents of seq->seq.s
              (and seq->qual.s) must be transfered to the corresponding members
              of sqzfastx_t. Making shure those buffers are big enough.
            */
            sqz_rememberseq(sqz, seq, 1);
            goto exit;
        }
        bases += l;
        memcpy(seqbuffer + offset, &l, B64);
        offset += B64;
        memcpy(seqbuffer + offset, seq->seq.s, l + 1);
        memcpy(qltbuffer + offset, seq->qual.s, l + 1);
        offset += l + 1;
        if ( maxlen <= l + 1 + B64 ) break;
        maxlen -= l + 1 + B64;
    }
    exit:
        sqz->n = n;
        sqz->bases = bases;
        sqz->offset = offset;
        return offset;
}


static uint64_t sqz_fastanblock(sqzfastx_t *sqz, kseq_t *seq)
{

    uint64_t offset = 0;
    uint64_t n      = 0;
    uint64_t l;
    uint64_t bases  = 0;
    uint64_t maxlen = LOAD_SIZE - B64;
    uint8_t *seqbuffer = sqz->seqbuffer;
    while ( (kseq_read(seq) >= 0) ) {
        l = seq->seq.l;
        n++;
        if (!sqz_loadname(sqz, seq)) {
            offset = 0;
            goto exit;
        }
        //Determine if current sequence can be loaded completely
        if (l + 1 > maxlen) {  //+1 for string end null
            memcpy(seqbuffer + offset, &l, B64);
            offset += B64;
            memcpy(seqbuffer + offset, seq->seq.s, maxlen);
            offset += maxlen;
            seqbuffer[offset++] = 0;
            bases += maxlen;
            sqz->endflag = 1;
            sqz->seqread = maxlen;
            sqz->prevlen = l;
            /*
              The sqzfastx_t struct must remember the sequence (and qualities)
              of this partially decoded record. Thus the contents of seq->seq.s
              (and seq->qual.s) must be transfered to the corresponding members
              of sqzfastx_t. Making shure those buffers are big enough.
            */
            sqz_rememberseq(sqz, seq, 0);
            goto exit;
        }
        bases += l;
        memcpy(seqbuffer + offset, &l, B64); //Copy sequence length
        offset += B64;
        memcpy(seqbuffer + offset, seq->seq.s, l + 1); //Copy sequence + NULL
        offset += l + 1;
        //Check that at least the length of the next sequence fits in buffer
        if ( maxlen <= l + 1 + B64 ) break;
        maxlen -= l + 1 + B64;
    }
    exit:
        sqz->n = n;
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
