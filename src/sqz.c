#include <stdio.h>
#include <zlib.h>
#include "squeezmalib.h"
#include "sqz.h"


int main(int argc, char *argv[]) {
    int ret = -1;
    if (argc < 2) goto exit;

    if (!sqz_compress(argv[1])) goto exit;

    ret = 0;
    exit:
        return ret;
}


char sqz_compress(const char *filename)
{
    int ret = 0;
    //Initialize data main sqz data structure
    sqzfastx_t *sqz = sqz_fastxinit(filename, LOAD_SIZE);
    if (!sqz) {
        fprintf(stderr, "[sqz ERROR: INIT] Failed to start data structures.\n");
        goto exit;
    }
    //Check for format
    switch (sqz->fmt) {
    case 1:
        //if (!sqz_fasta(sqz, seq)) goto exit;
        break;
    case 2:
        if (!sqz_squeezefastq(sqz)) {
            fprintf(stderr, "[sqz ERROR: COMPRESS] Failed to compress data.\n");
            goto exit;
        }
        break;
    }
    ret = 1;
    exit:
        sqz_kill(sqz);
        return ret;
}


/*
  Loads fastq data to buffer and compresses it with deflate
*/
char sqz_squeezefastq(sqzfastx_t *sqz)
{
    size_t batchsize = 0;
    size_t numseqs = 0;
    int i = 0;
    while ( (batchsize = sqz_loadfastq(sqz)) > 0 ) {
        i++;
        fprintf(stderr, "^%lu loaded bytes batch size %lu.\n", batchsize, sqz->n);
        numseqs += sqz->n;
        if (!sqz->endflag)
            fprintf(stderr, "Data ready for compression and flushing\n");
        if ( i == 10 )break;
    }
    fprintf(stderr, "processed %lu seqs\n", numseqs);
    return 1;
}
