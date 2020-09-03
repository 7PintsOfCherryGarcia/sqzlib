#include <stdio.h>
#include <stdlib.h>

#include "sqz_data.h"
#include "sqz_dcp.h"


long sqz_filesize(FILE *fp)
{
    fseek(fp, 0, SEEK_END);
    long s = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    return s;
}


sqzfastx_t *sqz_sqzinit(const char *filename, size_t bsize)
{
    sqzfastx_t *sqz = malloc(sizeof(sqzfastx_t));
    if (!sqz) {
        fprintf(stderr, "[squeezma ERROR] Memory error\n");
        return NULL;
    }
    sqz->filename = filename;
    sqz->n =       0;
    sqz->bases =   0;
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
    switch (fmt & 7) {
        case 0:
            fprintf(stderr,
                    "[sqzlib ERROR]: File %s of unknown format\n",
                    filename);
            free(sqz);
            sqz = NULL;
            break;
        //fasta
        case 1:
            fprintf(stderr,
                    "[sqzlib ERROR]: File %s already in decoded fasta format\n",
                    filename);
            free(sqz);
            sqz = NULL;
            break;
        case 2:
            fprintf(stderr,
                    "[sqzlib ERROR]: File %s already in decoded fastq format\n",
                    filename);
            free(sqz);
            sqz = NULL;
            break;
        case 5:
            fprintf(stderr, "[sqzlib INFO]: sqz encoded fasta file.\n");
            free(sqz);
            sqz = NULL;
            break;
        case 6:
            fprintf(stderr, "[sqzlib INFO]: sqz encoded fastq file.\n");
            sqz->fmt = fmt;
            sqz->qualbuffer = malloc(LOAD_SIZE + 1);
            if (!sqz->qualbuffer) {
                fprintf(stderr,
                        "[sqzlib ERROR] Failed to allocate quality buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->qualbuffer[LOAD_SIZE] = 0;
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



