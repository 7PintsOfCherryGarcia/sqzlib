#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "sqz_init.h"


sqzfastx_t *sqz_fastxinit(const char *filename, size_t buffersize)
{
    sqzfastx_t *sqz = malloc(sizeof(sqzfastx_t));
    if (!sqz) {
        fprintf(stderr, "[squeezma ERROR] Memory error\n");
        return NULL;
    }
    sqz->filename = filename;
    sqz->n = 0;
    sqz->offset  = 0;
    sqz->endflag = 0;
    sqz->rem =     0;
    sqz->toread =  0;
    sqz->prevlen = 0;
    //Get file format
    char fmt = sqz_getformat(filename);
    //Initialize kseq objects
    if (!sqz_kseqinit(sqz)) {
        free(sqz);
        return NULL;
    }
    switch (fmt) {
        case 0:
            free(sqz);
            sqz = NULL;
            break;
        //fasta
        case 1:
            fprintf(stderr, "[squeezma INFO] Detected FASTA format\n");
            sqz->fmt = fmt;
            //No need for quality buffer
            sqz->qualbuffer = NULL;
            sqz->seqbuffer = malloc(LOAD_SIZE + 1);
            if (!sqz->seqbuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate sequence buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer[LOAD_SIZE] = 0;
            sqz->namebuffer = malloc(1*1024*1024);
            if (!sqz->namebuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate name buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->namelen = 0;
            break;
        case 2:
            fprintf(stderr, "[squeezma INFO] Detected FASTQ format\n");
            sqz->fmt = fmt;
            sqz->qualbuffer = malloc(LOAD_SIZE);
            if (!sqz->qualbuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate quality buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer = malloc(LOAD_SIZE + 1);
            if (!sqz->seqbuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate sequence buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer[LOAD_SIZE] = 0;
            sqz->namebuffer = malloc(1*1024*1024);
            if (!sqz->namebuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate name buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->namelen = 0;
            break;
    }
    return sqz;
}


