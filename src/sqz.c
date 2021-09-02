#include "sqzlib.h"


typedef struct  {
    char ifile[256];
    char ofile[256];
    int  nthread;
    char  libfmt;
    char fqflag;
} sqzopts_t;


char sqz_threadlauncher(FILE *ofp,
                        const char *filename,
                        char fqflag,
                        int nthread,
                        uint8_t libfmt,
                        unsigned char fmt);

char sqz_deflatefastX(const char *ifile,
                      const char *ofile,
                      char fqflag,
                      uint8_t fmt,
                      uint8_t libfmt,
                      int nthread)
{
    char ret = 0;
    FILE *ofp = fopen(ofile, "wb");
    if ( !(ofp) ) return ret;

    if ( !sqz_filehead(fmt, libfmt, ofp) )
        goto exit;
    fflush(ofp);
    sqz_threadlauncher(ofp,
                       ifile,
                       fqflag,
                       nthread,
                       libfmt,
                       fmt);
    //TODO Replace the 0 (# sequences) sequence writting
    if ( !sqz_filetail(0, ofp) )
        goto exit;

    ret = 1;
    exit:
        fclose(ofp);
        return ret;
}


char sqz_compress(sqzopts_t opts)
{
    char ret = 0;
    uint8_t fmt = sqz_getformat(opts.ifile);
    switch (fmt & 7) {
        case 1:
            if (!sqz_deflatefastX(opts.ifile,
                                  opts.ofile,
                                  0,
                                  fmt,
                                  opts.libfmt,
                                  opts.nthread)) goto exit;
            break;
        case 2:
            if (!sqz_deflatefastX(opts.ifile,
                                  opts.ofile,
                                  1,
                                  fmt,
                                  opts.libfmt,
                                  opts.nthread)) goto exit;
            break;
        case 5:
            fprintf(stderr, "File %s alredy sqz encoded.\n", opts.ifile);
            goto exit;
        case 6:
            fprintf(stderr, "File %s alredy sqz encoded.\n", opts.ifile);
            goto exit;
    }
    ret = 1;
    exit:
        return ret;
}


uint8_t sqz_inflatefastX(FILE *ifp, FILE *ofp, char fqflag, uint8_t libfmt)
{
    uint8_t ret      = 0;
    uint8_t *outbuff = NULL;
    uint64_t dsize   = 0;
    int64_t size    = 0;
    sqzblock_t *blk  = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;
    size = sqz_filesize(ifp);
    outbuff = malloc(LOAD_SIZE);
    if (!outbuff) goto exit;
    fseek(ifp, HEADLEN, SEEK_SET);
    while ( ftell(ifp) < size )
        {
            if (!sqz_readblksize(blk, ifp, libfmt)) goto exit;
            do {
                dsize = sqz_fastXdecode(blk, outbuff, LOAD_SIZE, fqflag);
                fwrite(outbuff, 1, dsize, ofp);
                fflush(ofp);
            } while (blk->newblk);
        }
    ret = 1;
    exit:
        sqz_blkkill(blk);
        free(outbuff);
        return ret;
}


char sqz_decompress(sqzopts_t opts)
{
    char ret = 0;
    FILE *ifp = NULL;
    FILE *ofp = NULL;
    ofp = fopen(opts.ofile, "wb");
    if (!ofp)
        goto exit;
    ifp = fopen(opts.ifile, "rb");
    if (!ifp)
        goto exit;
    uint8_t fmt = sqz_getformat(opts.ifile);
    uint8_t libfmt = fmt >> 3;
    //Check for format
    switch (fmt & 7) {
        case 1:
            fprintf(stderr, "File %s already decoded.\n", opts.ifile);
            goto exit;
        case 2:
            fprintf(stderr, "File %s already decoded.\n", opts.ifile);
            goto exit;
        case 5:
            if (!sqz_inflatefastX(ifp, ofp, 0, libfmt)) {
                fprintf(stderr, "[sqz ERROR]: Failed to decode data.\n");
                goto exit;
            }
            break;
        case 6:
            if (!sqz_inflatefastX(ifp, ofp, 1, libfmt)) {
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
    fprintf(stderr, "\tsqz [options] <file>\n\n");
    fprintf(stderr, "\t\tOptions:\n");
    fprintf(stderr, "\t\t-d \t\tDecompress previously sqz compressed file.\n");
    fprintf(stderr, "\t\t-o <file>\tWrite to outputfile <file>.\n");
    fprintf(stderr,
            "\t\t-t <n>\t\tUse <n> compression/decompression threads.\n");
    fprintf(stderr,
            "\t\t-l <lib>\tUse <lib> compression library."\
            "\n\t\t\t\t\tDefault: -l zlib\n"\
            "\t\t\t\t\tOptions: zlib, zstd\n");
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
    char ret    = 1;      //Defaults to encode and compress
    char dflag  = 0;      //Decompression flag
    char oflag  = 1;      //Output flag
    char *lib   = NULL;   //Compression library to use
    int nthread = 1;
    int minargs = 2;
    while (( elem = getopt(argc, argv, "o:t:l:hd") ) >= 0) {
        switch(elem) {
            case 'd':
                dflag = 1;
                ret   = 2;
                minargs += 1;
                continue;
            case 'o':
                oflag = 0;
                strcpy(opts->ofile, optarg);
                minargs += 2;
                continue;
            case 't':
                nthread = atoi(optarg);
                minargs += 2;
                continue;
            case 'l':
                lib = optarg;
                minargs += 2;
                continue;
            case 'h':
                ret = 0;
                sqz_usage();
                goto exit;
            case '?':
                sqz_usage();
                ret = 0;
                goto exit;
        }

    }
    if (argc < minargs) {
        fprintf(stderr, "[sqz]: Please provide input file.\n\n");
        sqz_usage();
        ret = 0;
        goto exit;
    }
    if ( access(argv[argc - 1], F_OK) ) {
        fprintf(stderr, "[sqz]: %s, no such file.\n\n", argv[argc - 1]);
        sqz_usage();
        ret = 0;
        goto exit;
    }
    strcpy(opts->ifile, argv[argc - 1]);

    if (dflag) {
        if (oflag) strcpy(opts->ofile, "/dev/stdout");
    }
    else {
        if (oflag) {
            char *bname = sqz_basename( argv[argc - 1] );
            strcpy(opts->ofile, bname);
            strcat(opts->ofile, ".sqz");
        }
    }
    if (nthread <= 0) opts->nthread = 1;
    else opts->nthread = nthread;

    if (!lib) opts->libfmt = 1;
    else {
        if (!strcmp(lib, "zlib")) opts->libfmt = 1;
        else if (!strcmp(lib, "zstd")) opts->libfmt = 2;
        else {
            fprintf(stderr, "[sqz]: %s compression lib not supported\n\n", lib);
            sqz_usage();
            ret = 0;
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
    switch (sqz_ropts(argc, argv, &opts)) {
        case 0:
            goto exit;
        case 1:
            if ( !sqz_compress(opts) ) goto exit;
            break;
        case 2:
            if ( !sqz_decompress(opts)) goto exit;
            break;
    }
    ret = 0;
    exit:
        return ret;
}
