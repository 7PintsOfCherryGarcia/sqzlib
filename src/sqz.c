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
    char ret = 0;
    size_t batchsize = 0;
    size_t numseqs = 0;
    int i = 0;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) {
        goto exit;;
    }
    while ( (batchsize = sqz_loadfastq(sqz)) > 0 ) {
        i++;
        fprintf(stderr, "[sqz INFO]: %lu loaded bytes batch size %lu.\n",
                        2*batchsize,
                        sqz->n);
        if (!sqz_encode(sqz, blk)) {
            fprintf(stderr, "[sqz ERROR]: Encoding error.\n");
            goto exit;
        }
        numseqs += sqz->n;
        if (!sqz->endflag) {
            fprintf(stderr, "[sqz INFO]: Compress and flushing.\n");
            fprintf(stderr ,"%lu bases in block\n", sqz->bases);
            size_t cbytes = sqz_deflate(blk, 9);
            fprintf(stderr, "[sqz INFO]: Compressed to %lu bytes.\n", cbytes);
            blk->blksize = 0;
            sqz->bases = 0;
        }
        if ( i == 10 ) break;
    }
    fprintf(stderr, "processed %lu seqs\n", numseqs);
    ret = 1;
    exit:
        sqz_blkdestroy(blk);
        return ret;
}
