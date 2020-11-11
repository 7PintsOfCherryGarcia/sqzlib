#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include <stdlib.h>
#include "squeezmalib.h"
#include "sqz.h"

long sqz_filesize(FILE *fp);


int main(int argc, char *argv[])
{
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
            fprintf(stderr, "[sqz INFO]: output filename - %s\n", oname);
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
    switch (sqz->fmt & 7) {
        case 1:
            if (!sqz_squeezefasta(sqz, ofp)) {
                fprintf(stderr,
                        "[sqz ERROR: COMPRESS] Failed to compress data.\n");
                goto exit;
            }
            break;
        case 2:
            if (!sqz_squeezefastq(sqz, ofp)) {
                fprintf(stderr,
                        "[sqz ERROR: COMPRESS] Failed to compress data.\n");
                goto exit;
            }
            break;
        case 5:
            fprintf(stderr, "File %s alredy sqz encoded.\n", filename);
            goto exit;
        case 6:
            fprintf(stderr, "File %s alredy sqz encoded.\n", filename);
            goto exit;
    }
    ret = 1;
    exit:
        fclose(ofp);
        sqz_kill(sqz);
        return ret;
}


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

    //main compression loop
    while ( (batchsize += sqz_loadfastq(sqz)) > 0 ) {
        if (!sqz_encode(sqz, blk)) {
            fprintf(stderr, "[sqz ERROR]: Encoding error.\n");
            goto exit;
        }
        if (!sqz->endflag) {
            fprintf(stderr, "[sqz INFO]: Compress and flushing.\n");
            fprintf(stderr, "[sqz INFO]: Block info\n");
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


char sqz_squeezefasta(sqzfastx_t *sqz, FILE *ofp)
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

    //main compression loop
    while ( (batchsize += sqz_loadfasta(sqz)) > 0 ) {
        if (!sqz_fastaencode(sqz, blk)) {
            fprintf(stderr, "[sqz ERROR]: Encoding error.\n");
            goto exit;
        }
        if (!sqz->endflag) {
            fprintf(stderr, "[sqz INFO]: Compress and flushing.\n");
            fprintf(stderr, "[sqz INFO]: Block info\n");
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
    FILE *ofp = fopen("/dev/stdout", "wb");
    if (!ofp) {
        fprintf(stderr, "[sqz ERROR]: Failed to open %s for writing.\n", outname);
        return 0;
    }
    //Read head and tail
    sqzfastx_t *sqz = sqz_fastxinit(filename, LOAD_SIZE);
    if (!sqz | !(sqz->fmt | 4) ) {
        //fprintf(stderr, "sqz format flag: %d\n", sqz->fmt);
        fprintf(stderr, "[sqz ERROR]: Input is not sqz encoded.\n");
        goto exit;
    }
    //Check for format
    switch (sqz->fmt & 7) {
        case 1:
            fprintf(stderr, "File %s already decoded.\n", filename);
            goto exit;
        case 2:
            fprintf(stderr, "File %s already decoded.\n", filename);
            goto exit;
        case 5:
            if (!sqz_spreadfasta(sqz, ofp)) {
                fprintf(stderr,
                        "[sqz ERROR]: Failed to decode data.\n");
                goto exit;
            }
            break;
        case 6:
            if (!sqz_spreadfastq(sqz, ofp)) {
                fprintf(stderr,
                        "[sqz ERROR]: Failed to decode data.\n");
                goto exit;
            }
            break;
        }
    ret = 1;
    exit:
        sqz_kill(sqz);
        fclose(ofp);
        return ret;
}


char sqz_spreadfastq(sqzfastx_t *sqz, FILE *ofp)
{
    char ret = 0;
    //read first blocks compressed and decompressed size
    size_t bytes = sizeof(size_t);
    size_t cmpsize;
    size_t dcpsize;
    size_t tmp = 0;

    FILE *fp = fopen(sqz->filename, "rb");
    if (!fp) return ret;

    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) {
        goto exit;;
    }

    long size = sqz_filesize(fp) - 16;
    fseek(fp, bytes, SEEK_SET);
    uint64_t nblock = 0;
    while ( ftell(fp) < size )
        {
            nblock++;
            fprintf(stderr, "[sqz INFO]: In block %lu\n", nblock);
            //TODO move reading to a separate routine that loads
            //all the necessary data to inflate a zlib block
            tmp = fread(&dcpsize, bytes, 1, fp);
            tmp += fread(&cmpsize, bytes, 1, fp);
            if ( (cmpsize != fread(blk->cmpbuff, 1, cmpsize, fp)) |
                 (tmp != 2) )
            {
                fprintf(stderr, "ERROR: tmp: %lu cmpsize: %lu\n", tmp, cmpsize);
                fprintf(stderr, "ERROR reading\n");
                break;
            }
            blk->cmpsize = cmpsize;
            blk->blksize = LOAD_SIZE*2;
            //Inflate
            size_t cbytes = sqz_inflate(blk);
            if (cbytes != dcpsize) {
                    fprintf(stderr, "[sqz ERROR]: Corrupt data encountered.\n");
                    goto exit;
            }
            fprintf(stderr, "\tDeflate size: %lu\n", cbytes);
            //Move all coded data to an decode buffer
            sqz_fastqdecode(blk->codebuff, dcpsize);
        }

    ret = 1;
    exit:
        sqz_blkdestroy(blk);
        fclose(fp);
    return ret;
}


char sqz_spreadfasta(sqzfastx_t *sqz, FILE *ofp)
{
    char ret = 0;
    //read first blocks compressed and decompressed size
    size_t bytes = sizeof(size_t);
    size_t cmpsize;
    size_t dcpsize;
    size_t tmp = 0;

    FILE *fp = fopen(sqz->filename, "rb");
    if (!fp) return ret;

    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) {
        goto exit;;
    }

    long size = sqz_filesize(fp) - 16;
    fseek(fp, bytes, SEEK_SET);
    uint64_t nblock = 0;
    while ( ftell(fp) < size )
        {
            nblock++;
            fprintf(stderr, "[sqz INFO]: In block %lu\n", nblock);
            //TODO move reading to a separate routine that loads
            //all the necessary data to inflate a zlib block
            tmp = fread(&dcpsize, bytes, 1, fp);
            tmp += fread(&cmpsize, bytes, 1, fp);
            if ( (cmpsize != fread(blk->cmpbuff, 1, cmpsize, fp)) |
                 (tmp != 2) )
            {
                fprintf(stderr, "ERROR: tmp: %lu cmpsize: %lu\n", tmp, cmpsize);
                fprintf(stderr, "ERROR reading\n");
                break;
            }
            blk->cmpsize = cmpsize;
            blk->blksize = LOAD_SIZE*2;
            //Inflate
            size_t cbytes = sqz_inflate(blk);
            if (cbytes != dcpsize) {
                    fprintf(stderr, "[sqz ERROR]: Corrupt data encountered.\n");
                    goto exit;
            }
            fprintf(stderr, "\tDeflate size: %lu\n", cbytes);
            //Move all coded data to an decode buffer
            sqz_fastadecode(blk->codebuff, dcpsize);
        }

    ret = 1;
    exit:
        sqz_blkdestroy(blk);
        fclose(fp);
    return ret;
}



char sqz_ropts(int argc, char **argv)
{
    int elem;
    char ret = 1;
    while (( elem = getopt(argc, argv, "d:h") ) >= 0) {
        switch(elem) {
        case 'd':
            ret = 2;
            continue;
        case 'h':
            ret = 0;
            sqz_usage();
            break;
        case '?':
            sqz_usage();
            ret = 0;
            break;
        }

    }
    return ret;
}


void sqz_usage()
{
    fprintf(stderr, "[sqz INFO]: Usage:\n");
    fprintf(stderr, "\tsqz [options] <file>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "\t-d <file>\tDecompress previously sqz compressed file.\n");
    fprintf(stderr, "\t-h  \t\tThis help message.\n");
}
