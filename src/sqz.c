#include <stdio.h>
#include <string.h>
#include "squeezmalib.h"
#include "sqz.h"


int main(int argc, char *argv[]) {
    int ret = -1;
    if (argc < 2) goto exit;
    char oname[256];
    strcpy(oname, argv[1]);
    strcat(oname, ".sqz");

    fprintf(stderr, "[sqz INFO]: output filename - %s\n", oname);

    if (!sqz_compress(argv[1], oname)) goto exit;

    ret = 0;
    exit:
        return ret;
}


char sqz_compress(const char *filename, const char *outname)
{
    int ret = 0;
    //Open output file handle
    //TODO move to function
    FILE *ofp = fopen(outname, "wb");
    if (!ofp) {
        fprintf(stderr, "[sqz ERROR]: Failed to open %s for writing.\n", outname);
        return 0;
    }

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
        if (!sqz_squeezefastq(sqz, ofp)) {
            fprintf(stderr, "[sqz ERROR: COMPRESS] Failed to compress data.\n");
            goto exit;
        }
        break;
    }
    ret = 1;
    exit:
        fclose(ofp);
        sqz_kill(sqz);
        return ret;
}


/*
  Loads fastq data to buffer and compresses it with deflate
*/
char sqz_squeezefastq(sqzfastx_t *sqz, FILE *ofp)
{
    char ret = 0;
    size_t batchsize = 0;
    size_t numseqs = 0;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) {
        goto exit;;
    }
    while ( (batchsize += sqz_loadfastq(sqz)) > 0 ) {
        numseqs += sqz->n;
        if (!sqz_encode(sqz, blk)) {
            fprintf(stderr, "[sqz ERROR]: Encoding error.\n");
            goto exit;
        }
        if (!sqz->endflag) {
            fprintf(stderr, "[sqz INFO]: Compress and flushing.\n");
            fprintf(stderr, "[sqz INFO]: Block info\n");
            fprintf(stderr, "\t\toffset - %lu bytes\n", batchsize);
            fprintf(stderr, "\t\tsequences - %lu\n", sqz->n);
            fprintf(stderr ,"\t\tbases - %lu\n", sqz->bases);
            size_t cbytes = sqz_deflate(blk, 9);
            if ( !sqz_zlibcmpdump(blk, cbytes, ofp) ) {
                fprintf(stderr, "[sqz ERROR]: IO error");
            }
            fprintf(stderr, "[sqz INFO]: Compressed to %lu bytes.\n", cbytes);
            blk->blksize = 0;
            sqz->bases = 0;
            batchsize = 0;
        }
    }
    fprintf(stderr, "processed %lu seqs\n", numseqs);
    ret = 1;
    exit:
        sqz_blkdestroy(blk);
        return ret;
}
