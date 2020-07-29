#include <stdio.h>
#include <zlib.h>
#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread)
#include "sqz_kseq.h"


char sqz_getformat(const char *filename)
{
    char ret = 0;
    gzFile fp = gzopen(filename, "r");
    //ERROR
    if (!fp) {
        fprintf(stderr, "[libsqz: ERROR] Failed to open file: %s\n", filename);
        return ret;
    }
    kseq_t *seq = kseq_init(fp);
    //ERROR
    if (!seq) {
        gzclose(fp);
        return ret;
    }
    int l = kseq_read(seq);
    //ERROR
    if (l < 0) {
        fprintf(stderr, "[libsqz: ERROR] sequence file format not recognized\n");
        goto exit;
    }
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


char sqz_kseqinit(sqzfastx_t *sqz)
{
    char ret = 0;
    sqz->fp = gzopen(sqz->filename, "r");
    //ERROR
    if (!sqz->fp) {
        fprintf(stderr,
                "[libsqz ERROR: zlib] Failed to open file: %s\n", sqz->filename);
        goto exit;
    }
    sqz->seq = kseq_init(sqz->fp);
    //ERROR
    if (!sqz->seq) {
        gzclose(sqz->fp);
        fprintf(stderr,
                "[libsqz ERROR: kseq] Failed to init kseq: %s\n", sqz->filename);
        goto exit;
    }
    ret = 1;
    exit:
        return ret;
}


size_t sqz_loadfastq(sqzfastx_t *sqz)
{
    fprintf(stderr, "endflag: %d remaining: %lu\n", sqz->endflag, sqz->rem);
    if (!sqz->endflag) return sqz_newblock(sqz);
    return sqz_endblock(sqz);
}


size_t sqz_newblock(sqzfastx_t *sqz)
{
    size_t bleftover;
    size_t offset = 0;
    size_t lenbytes = sizeof(sqz->seq->seq.l);
    size_t n = 0;
    //int l;
    //Loop over sequences and load data into sqzfastx_t struct
    while ( kseq_read(sqz->seq) >= 0) {
        n++;
        /*
          When a buffer is filled (Can't hold the entire kseq sequence + the
          lenth of the next sequence), the length of the current sequence is
          stored as well as any bases the buffer can accomodate.
        */
        if (offset + sqz->seq->seq.l + 1 + lenbytes + lenbytes > LOAD_SIZE) {
            //Compute how much buffer is available
            bleftover = LOAD_SIZE - offset;
            //Copy sequence length data
            memcpy(sqz->seqbuffer + offset, &(sqz->seq->seq.l), lenbytes);
            offset += lenbytes;
            bleftover -= lenbytes;
            /*
            Compute how much sequence can be loaded to the buffer.
            There are two option: as much sequence as leftover buffer, or the
            entire sequence
            */
            sqz->rem = bleftover<sqz->seq->seq.l?bleftover:sqz->seq->seq.l + 1;
            //Copy as much seq data as we can fit in remaining buffer
            memcpy(sqz->seqbuffer + offset, sqz->seq->seq.s, sqz->rem);
            memcpy(sqz->qualbuffer + offset, sqz->seq->qual.s, sqz->rem);
            offset += sqz->rem;
            sqz->n = n;
            sqz->rem = sqz->seq->seq.l + 1 - sqz->rem;
            if (sqz->rem != 0) sqz->endflag = 1;
            return offset;

        }
        //copy sequence data into buffers
        else {
            //Copy sequence length data
            memcpy(sqz->seqbuffer + offset, &(sqz->seq->seq.l), lenbytes);
            offset += lenbytes;
            //Copy sequence string plus terminating null byte
            memcpy(sqz->seqbuffer+offset, sqz->seq->seq.s, sqz->seq->seq.l+1);
            //Copy quality string plus terminating null byte
            memcpy(sqz->qualbuffer+offset, sqz->seq->qual.s, sqz->seq->seq.l+1);
            offset += sqz->seq->seq.l + 1;
        }
    }
    if (n != 0) {
        sqz->n = n;
    //    sqz_encodencompress(sqz, offset);
    //    if(!sqz_cmpnflush(sqz)) return 0;
    }
    return offset;
}


size_t sqz_endblock(sqzfastx_t *sqz)
{
    size_t lsize;
    //buffer can be completely filled with current sequence
    if (sqz->rem >= LOAD_SIZE) {
        memcpy(sqz->seqbuffer,
               sqz->seq->seq.s + (sqz->seq->seq.l + 1 - sqz->rem),
               LOAD_SIZE);
        memcpy(sqz->qualbuffer,
               sqz->seq->qual.s + (sqz->seq->qual.l + 1 - sqz->rem),
               LOAD_SIZE);
        sqz->rem -= LOAD_SIZE;
        lsize = LOAD_SIZE;
        return lsize;
    }
    //Rest of sequence can go into buffer
    memcpy(sqz->seqbuffer,
           sqz->seq->seq.s + (sqz->seq->seq.l + 1 - sqz->rem),
           sqz->rem);
    memcpy(sqz->qualbuffer,
           sqz->seq->qual.s + (sqz->seq->seq.l + 1 - sqz->rem),
           sqz->rem);
    lsize = sqz->rem;
    sqz->endflag = 0;
    sqz->rem = 0;
    return lsize;
}


void sqz_kill(sqzfastx_t *sqz)
{
    if (sqz) {
        gzclose(sqz->fp);
        kseq_destroy(sqz->seq);
        free(sqz->seqbuffer);
        free(sqz->namebuffer);
        if (sqz->fmt == 2) free(sqz->qualbuffer);
        free(sqz);
    }
}
