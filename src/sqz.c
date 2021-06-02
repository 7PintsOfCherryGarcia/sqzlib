#define SQZLIB
#define KLIB
#include "sqzlib.h"
typedef struct  {
    char ifile[256];
    char ofile[256];
} sqzopts_t;


char sqz_squeezefastX(sqzfastx_t *sqz, FILE *ofp, char fqflag)
{
    char ret = 0;
    uint64_t lbytes = 0;
    uint64_t cbytes = 0;
    uint64_t numseqs = 0;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;
    if ( !sqz_filehead(sqz, ofp) ) goto exit;
    while ( (lbytes += sqz_loadfastX(sqz, fqflag)) > 0 ) {
        if (!sqz_fastqencode(sqz, blk)) {
            fprintf(stderr, "[sqz ERROR]: Encoding error.\n");
            goto exit;
        }
        if (sqz->cmpflag) {
            numseqs += sqz->n;
            cbytes = sqz_deflate(blk, 9);
            if ( !sqz_zlibcmpdump(blk, cbytes, ofp) ) {
                fprintf(stderr, "[sqz ERROR]: Failed to write to output.\n");
                goto exit;
            }
            blk->blkpos  = 0;
            sqz->bases   = 0;
            sqz->namepos = 0;
            sqz->endflag = 0;
            sqz->cmpflag = 0;
            blk->newblk  = 1;
            sqz->n = 0;
            lbytes = 0;
        }
    }
    fprintf(stderr, "[sqz INFO]: processed %lu sequences\n", numseqs);
    if ( !sqz_filetail(numseqs, ofp) ) {
        fprintf(stderr, "[sqz ERROR]: Failed to finish sqz file.\n");
        goto exit;
    }
    ret = 1;
    exit:
        sqz_blkdestroy(blk);
        return ret;
}


char sqz_compress(const char *filename, const char *outname)
{
    char ret = 0;
    FILE *ofp = fopen(outname, "wb");
    if (!ofp) return ret;

    sqzfastx_t *sqz = sqz_fastxinit(filename, LOAD_SIZE);
    if (!sqz) goto exit;
    switch (sqz->fmt & 7) {
        case 1:
            if (!sqz_squeezefastX(sqz, ofp, 0)) goto exit;
            break;
        case 2:
            if (!sqz_squeezefastX(sqz, ofp, 1)) goto exit;
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


char sqz_inflatefastX(FILE *ifp, FILE *ofp, char fqflag)
{
    char ret = 0;
    uint8_t *outbuff = NULL;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;
    long size = sqz_filesize(ifp);
    outbuff = malloc(LOAD_SIZE);
    if (!outbuff) goto exit;
    uint64_t dsize = 0;
    fseek(ifp, HEADLEN, SEEK_SET);
    while ( ftell(ifp) < size )
        {
            if (!sqz_readblksize(blk, ifp)) goto exit;
            do {
                dsize = sqz_fastXdecode(blk, outbuff, LOAD_SIZE, fqflag);
                fwrite(outbuff, 1, dsize, ofp);
                fflush(ofp);
            } while (blk->newblk);
        }
    ret = 1;
    exit:
        sqz_blkdestroy(blk);
        free(outbuff);
        return ret;
}


char sqz_decompress(const char *filename, const char *outname)
{
    char ret = 0;
    FILE *ifp = NULL;
    FILE *ofp = NULL;
    ofp = fopen(outname, "wb");
    if (!ofp)
        goto exit;
    ifp = fopen(filename, "rb");
    if (!ifp)
        goto exit;
    uint8_t fmt = sqz_getformat(filename);
    //Check for format
    switch (fmt & 7) {
        case 1:
            fprintf(stderr, "File %s already decoded.\n", filename);
            goto exit;
        case 2:
            fprintf(stderr, "File %s already decoded.\n", filename);
            goto exit;
        case 5:
            if (!sqz_inflatefastX(ifp, ofp, 0)) {
                fprintf(stderr, "[sqz ERROR]: Failed to decode data.\n");
                goto exit;
            }
            break;
        case 6:
            if (!sqz_inflatefastX(ifp, ofp, 1)) {
                fprintf(stderr, "[sqz ERROR]: Failed to decode data.\n");
                goto exit;
            }
            break;
        }
    ret = 1;
    exit:
        if(ifp) fclose(ifp);
        if(ofp) fclose(ofp);
        return ret;
}


void sqz_usage()
{
    fprintf(stderr, "[sqz]: Usage:\n\n");
    fprintf(stderr, "\tsqz [options] <file>\n");
    fprintf(stderr, "\t\tOptions:\n");
    fprintf(stderr, "\t\t-d \t\tDecompress previously sqz compressed file.\n");
    fprintf(stderr, "\t\t-o <file>\tWrite to outputfile <file>.\n");
    fprintf(stderr, "\t\t-h \t\tThis help message.\n\n");
}


char *sqz_basename(char *namestr)
{
    char *tok = strtok(namestr, "/");
    char *ret = tok;
    while (NULL != (tok = strtok(NULL, "/")))
        ret = tok;
    return ret;
}


char sqz_ropts(int argc, char **argv, sqzopts_t *opts)
{
    int elem;
    char ret   = 1;      //Defaults to encode and compress
    char dflag = 0;      //Decompression flag
    char oflag = 1;      //Output flag
    while (( elem = getopt(argc, argv, "o:hd") ) >= 0) {
        switch(elem) {
            case 'd':
                dflag = 1;
                ret   = 2;
                continue;
            case 'o':
                oflag = 0;
                strcpy(opts->ofile, optarg);
                continue;
            case 'h':
                ret = 0;
                sqz_usage();
                goto exit;
            case '?':
                sqz_usage();
                ret = 0;
                break;
        }

    }
    if (dflag) {
        if (oflag) strcpy(opts->ofile,"/dev/stdout");
    }
    else {
        if (oflag) {
            char *bname;
            bname = sqz_basename( argv[argc - 1] );
            strcpy(opts->ofile, bname);
            strcat(opts->ofile, ".sqz");
        }
    }
    exit:
        return ret;
}


int main(int argc, char *argv[])
{
    int ret = 1;
    if (argc < 2) {
        sqz_usage();
        goto exit;
    }
    sqzopts_t opts;
    strcpy(opts.ifile, argv[argc - 1]);
    switch (sqz_ropts(argc, argv, &opts)) {
        case 0:
            goto exit;
        case 1:
            if (!sqz_compress(opts.ifile, opts.ofile)) goto exit;
            break;
        case 2:
            if (!sqz_decompress(opts.ifile, opts.ofile)) goto exit;
    }
    ret = 0;
    exit:
        if (ret)
            fprintf(stderr, "[sqz]: sqz failed.\n");
        return ret;
}
