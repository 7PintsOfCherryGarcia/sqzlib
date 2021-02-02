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
            {
            fprintf(stderr,
                    "[sqzlib ERROR]: File %s of unknown format\n",
                    filename);
            free(sqz);
            sqz = NULL;
            break;
            }
        case 1:
            {
            //fprintf(stderr, "[sqzlib INFO]: Detected fastA format %u\n", fmt);
            sqz->fmt = fmt;
            sqz->qualbuffer = NULL;
            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                fprintf(stderr,
                        "[sqzlib ERROR]: Failed to allocate sequence buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer[bsize] = '\0';
            sqz->namebuffer = malloc(NAME_SIZE);
            if (!sqz->namebuffer) {
                fprintf(stderr,
                        "[sqzlib ERROR] Failed to allocate name buffer\n");
                free(sqz);
                return NULL;
            }
            break;
            }
        case 2:
            {
            //fprintf(stderr, "[sqzlib INFO]: Detected fastQ format.\n");
            sqz->fmt = fmt;
            sqz->qualbuffer = malloc(bsize + 1);
            if (!sqz->qualbuffer) {
                fprintf(stderr,
                        "[sqzlib ERROR] Failed to allocate quality buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->qualbuffer[bsize] = 0;
            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate sequence buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer[bsize] = 0;
            sqz->namebuffer = malloc(1*1024*1024);
            if (!sqz->namebuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate name buffer\n");
                free(sqz);
                return NULL;
            }
            break;
            }
        case 5:
            {
            //fprintf(stderr,
            //        "[sqzlib INFO]: Detected sqz encoded FASTA format.\n");
            sqz->fmt = fmt;
            //Added when working on kseq compatibility
            {
            sqz->readbuffer = malloc(bsize + 1);
            if (!sqz->readbuffer) {
                fprintf(stderr,
                        "[sqzlib ERROR] Failed to allocate quality buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->readbuffer[bsize] = 0;
            }
            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate sequence buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer[bsize] = 0;
            sqz->namebuffer = malloc(1*1024*1024);
            if (!sqz->namebuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate name buffer\n");
                free(sqz);
                return NULL;
            }
            break;
            }
        case 6:
            {
            //fprintf(stderr,
            //        "[sqzlib INFO]: Detected sqz encoded FASTQ format.\n");
            sqz->fmt = fmt;
            //Added when working on kseq compatibility
            {
            sqz->readbuffer = malloc(bsize + 1);
            if (!sqz->readbuffer) {
                fprintf(stderr,
                        "[sqzlib ERROR] Failed to allocate quality buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->readbuffer[bsize] = '\0';
            }
            //
            sqz->qualbuffer = malloc(bsize + 1);
            if (!sqz->qualbuffer) {
                fprintf(stderr,
                        "[sqzlib ERROR] Failed to allocate quality buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->qualbuffer[bsize] = 0;
            sqz->seqbuffer = malloc(bsize + 1);
            if (!sqz->seqbuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate sequence buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer[bsize] = 0;
            sqz->namebuffer = malloc(1*1024*1024);
            if (!sqz->namebuffer) {
                fprintf(stderr,
                        "[squeezma ERROR] Failed to allocate name buffer\n");
                free(sqz);
                return NULL;
            }
            break;
            }
    }
    return sqz;
}


