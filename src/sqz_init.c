#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include "squeezmalib.h"
#include "sqz_init.h"

sqzfastx_t *sqz_fastxinit(const char *filename, size_t buffersize)
{
    sqzfastx_t *sqz = malloc(sizeof(sqzfastx_t));
    if (!sqz) {
        fprintf(stderr, "[squeezma ERROR] Memory error\n");
        return NULL;
    }
    sqz->n = 0;
    sqz->endflag = 0;
    sqz->toread = 0;
    sqz->prevlen = 0;
    //Get file format
    char fmt = sqz_getformat(filename);
    switch (fmt) {
        case 0:
            free(sqz);
            return NULL;
        //fasta
        case 1:
            fprintf(stderr, "[squeezma INFO] Detected FASTA format\n");
            sqz->fmt = fmt;
            //No need for quality buffer
            sqz->qualbuffer = NULL;
            sqz->seqbuffer = malloc(LOAD_SIZE + 1);
            if (!sqz->seqbuffer) {
                fprintf(stderr, "[squeezma ERROR] Failed to allocate sequence buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer[LOAD_SIZE] = 0;
            sqz->seqlen = 0;
            sqz->namebuffer = malloc(1*1024*1024);
            if (!sqz->namebuffer) {
                fprintf(stderr, "[squeezma ERROR] Failed to allocate name buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->namelen = 0;
            return sqz;
        case 2:
            fprintf(stderr, "[squeezma INFO] Detected FASTQ format\n");
            sqz->fmt = fmt;
            sqz->qualbuffer = malloc(LOAD_SIZE);
            if (!sqz->qualbuffer) {
                fprintf(stderr, "[squeezma ERROR] Failed to allocate quality buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer = malloc(LOAD_SIZE + 1);
            if (!sqz->seqbuffer) {
                fprintf(stderr, "[squeezma ERROR] Failed to allocate sequence buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer[LOAD_SIZE] = 0;
            sqz->seqlen = 0;
            sqz->namebuffer = malloc(1*1024*1024);
            if (!sqz->namebuffer) {
                fprintf(stderr, "[squeezma ERROR] Failed to allocate name buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->namelen = 0;
            return sqz;
    }
    return NULL;
}


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


void sqz_kill(sqzfastx_t *sqz)
{
    if (sqz) {
        free(sqz->seqbuffer);
        free(sqz->namebuffer);
        if (sqz->fmt == 2) free(sqz->qualbuffer);
        free(sqz);
    }
}
