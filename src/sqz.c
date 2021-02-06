#define SQZLIB
#define KLIB
#include "sqzlib.h"
#include "sqz.h"


int main(int argc, char *argv[])
{
    int ret = -1;
    if (argc < 2) {
        sqz_usage();
        goto exit;
    }
    char iname[256];
    strcpy(iname, argv[1]);
    char oname[256];
    char *bname;
    char *end;
    switch (sqz_ropts(argc, argv)) {
        case 0:
            goto exit;
        case 1:
            bname = sqz_basename( argv[argc - 1] );
            strcpy(oname, bname);
            strcat(oname, ".sqz");
            if (!sqz_compress(iname, oname)) goto exit;
            break;
        case 2:
            strcpy(oname, argv[argc - 1]);
            end = strrchr(oname, '.');
            if (end) *end = '\0';
            if (!sqz_decompress(argv[argc - 1], oname)) goto exit;
    }
    ret = 0;
    exit:
        return ret;
}


char sqz_compress(const char *filename, const char *outname)
{
    char ret = 0;
    FILE *ofp = fopen(outname, "wb");
    if (!ofp) return 0;
    sqzfastx_t *sqz = sqz_fastxinit(filename, LOAD_SIZE);
    if (!sqz) goto exit;
    switch (sqz->fmt & 7) {
        case 1:
            if (!sqz_squeezefasta(sqz, ofp)) goto exit;
            break;
        case 2:
            if (!sqz_squeezefastq(sqz, ofp)) goto exit;
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
    uint64_t lbytes = 0;
    uint64_t cbytes = 0;
    uint64_t numseqs = 0;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;

    if ( !sqz_filehead(sqz, ofp) ) goto exit;

    while ( (lbytes += sqz_loadfastq(sqz)) > 0 ) {
        if (!sqz_fastqencode(sqz, blk)) {
            fprintf(stderr, "[sqz ERROR]: Encoding error.\n");
            goto exit;
        }
        if (sqz->cmpflag) {
            numseqs += sqz->n;
            cbytes = sqz_deflate(blk, 9);
            if ( !sqz_zlibcmpdump(blk, cbytes, ofp) ) {
                fprintf(stderr, "[sqz ERROR]: IO error");
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
    uint64_t lbytes = 0;
    uint64_t cbytes  = 0;
    uint64_t numseqs = 0;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;;

    if ( !sqz_filehead(sqz, ofp) ) goto exit;
    while ( (lbytes += sqz_loadfasta(sqz)) > 0 ) {
        if (!sqz_fastaencode(sqz, blk)) {
            fprintf(stderr, "[sqz ERROR]: Encoding error.\n");
            goto exit;
        }
        if (sqz->cmpflag) {
            numseqs += sqz->n;
            cbytes = sqz_deflate(blk, 9);
            if ( !sqz_zlibcmpdump(blk, cbytes, ofp) ) {
                fprintf(stderr, "[sqz ERROR]: IO error");
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
    FILE *ifp = NULL;
    FILE *ofp = NULL;
    ofp = fopen("/dev/stdout", "wb");
    if (!ofp)
        goto exit;
    ifp = fopen(filename, "rb");
    if (!ifp)
        goto exit;

    uint8_t fmt = sqz_getformat(filename);
    //Check for format
    switch (fmt & 7) {
        case 1:
            {
            fprintf(stderr, "File %s already decoded.\n", filename);
            goto exit;
            }
        case 2:
            {
            fprintf(stderr, "File %s already decoded.\n", filename);
            goto exit;
            }
        case 5:
            {
            if (!sqz_spreadfasta(ifp, ofp)) {
                fprintf(stderr,
                        "[sqz ERROR]: Failed to decode data.\n");
                goto exit;
            }
            break;
            }
        case 6:
            {
            if (!sqz_spreadfastq(ifp, ofp)) {
                fprintf(stderr,
                        "[sqz ERROR]: Failed to decode data.\n");
                goto exit;
            }
            break;
            }
        }
    ret = 1;
    exit:
        if(ifp) fclose(ifp);
        if(ofp) fclose(ofp);
        return ret;
}


char sqz_spreadfastq(FILE *ifp, FILE *ofp)
{
    char ret = 0;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;
    long size = sqz_filesize(ifp);
    uint8_t *outbuff = malloc(LOAD_SIZE);
    uint64_t dsize = 0;
    if (!outbuff) goto exit;
    fseek(ifp, HEADLEN, SEEK_SET);
    while ( ftell(ifp) < size )
    {
        if (!sqz_readblksize(blk, ifp)) goto exit;
        //Decode block
        do {
            dsize = sqz_fastXdecode(blk, outbuff, LOAD_SIZE, 1);
            fwrite(outbuff, 1, dsize, ofp);
        } while (blk->blkpos);
    }
    ret = 1;
    exit:
        sqz_blkdestroy(blk);
    return ret;
}


char sqz_spreadfasta(FILE *ifp, FILE *ofp)
{
    char ret = 0;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;;
    long size = sqz_filesize(ifp) - 16;
    uint8_t *outbuff = malloc(LOAD_SIZE);
    uint64_t dsize = 0;
    if (!outbuff) goto exit;
    fseek(ifp, B64, SEEK_SET);
    while ( ftell(ifp) < size ) {
        //Decompress blk
        if (!sqz_readblksize(blk, ifp)) goto exit;
        //Decode block
        do {
            dsize = sqz_fastXdecode(blk, outbuff, LOAD_SIZE, 0);
            fwrite(outbuff, 1, dsize, ofp);
        } while (blk->blkpos);
    }
    ret = 1;
    exit:
        sqz_blkdestroy(blk);
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


char *sqz_basename(char *namestr)
{
    char *tok = strtok(namestr, "/");
    char *ret = tok;
    while (NULL != (tok = strtok(NULL, "/")))
        ret = tok;
    return ret;
}


