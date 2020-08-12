#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "squeezmalib.h"
#include "sqz.h"


int main(int argc, char *argv[]) {
    int ret = -1;
    if (argc < 2) goto exit;
    char oname[256];
    char *end;
    switch (sqz_ropts(argc, argv)) {
        case 0:
            goto exit;
        case 1:
            strcpy(oname, argv[1]);
            strcat(oname, ".sqz");
            fprintf(stderr, "[sqz INFO]: output filename - %s\n", oname);
            if (!sqz_compress(argv[1], oname)) goto exit;
            break;
        case 2:
            //TODO check for no file overwriting
            strcpy(oname, argv[argc - 1]);
            end = strrchr(oname, '.');
            if (end)
                *end = '\0';
            fprintf(stderr, "%s\n", oname);
            if (!sqz_decompress(argv[argc - 1], oname)) goto exit;
    }

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
                fprintf(stderr,
                        "[sqz ERROR: COMPRESS] Failed to compress data.\n");
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
    if ( !sqz_filehead(sqz, ofp) ) {
        fprintf(stderr, "[sqz ERROR]: IO error");
        goto exit;
    }

    while ( (batchsize += sqz_loadfastq(sqz)) > 0 ) {
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
            numseqs += sqz->n;

            size_t cbytes = sqz_deflate(blk, 9);
            if ( !sqz_zlibcmpdump(blk, cbytes, ofp) ) {
                fprintf(stderr, "[sqz ERROR]: IO error");
                goto exit;
            }
            fprintf(stderr, "[sqz INFO]: Compressed to %lu bytes.\n", cbytes);
            blk->blksize = 0;
            sqz->bases = 0;
            batchsize = 0;
        }
    }
    fprintf(stderr, "[sqz INFO]: processed %lu sequences\n", numseqs);
    if ( !sqz_filetail(numseqs, ofp) ) {
        fprintf(stderr, "[sqz ERROR]: IO error");
        goto exit;
    }
    ret = 1;
    exit:
        sqz_blkdestroy(blk);
        return ret;
}


char sqz_decompress(const char *filename, const char *outname)
{
    char ret = 0;
    fprintf(stderr, "%s %s\n", filename, outname);
    //Read head and tail
    sqzfastx_t *sqz = sqz_sqzinit(filename, LOAD_SIZE);
    if (!sqz) {
        fprintf(stderr, "[sqz ERROR: INIT] Failed to start data structures.\n");
        goto exit;
    }

    ret = 1;
    exit:
        return ret;
}

char sqz_ropts(int argc, char **argv)
{
    int elem;
    while (( elem = getopt(argc, argv, "d:h") ) >= 0) {
        switch(elem) {
        case 'd':
            return 2;
        case 'h':
            sqz_usage();
            return 0;
        case '?':
            sqz_usage();
            return 0;
        }

    }
    return 1;
}


void sqz_usage()
{
    fprintf(stderr, "[sqz INFO]: Usage:\n");
    fprintf(stderr, "\tsqz [options] <file>\n");
}
