#include <stdio.h>
#include <stdlib.h>

#include "sqz_dcp.h"
#include "sqz_data.h"


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
    switch (fmt & 7) {
        case 0:
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
            free(sqz);
            sqz = NULL;
    }
    return sqz;
}
