#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>
#include <unistd.h>

typedef struct  {
    char ifile[256];
    char ofile[256];
    int  nthread;
    char libfmt;
    int  idx;
} sqzopts_t;


uint8_t sqz_compress(const char *i,
                     const char *o,
                     uint8_t lib,
                     uint8_t nthread,
                     uint8_t idx);

uint8_t sqz_decompress(const char *i, const char *o, uint8_t nthread);


char sqz_cmpr(sqzopts_t opts)
{
    return sqz_compress(opts.ifile,
                        opts.ofile,
                        opts.libfmt,
                        opts.nthread,
                        opts.idx);
}


char sqz_dcmp(sqzopts_t opts)
{
    return sqz_decompress(opts.ifile, opts.ofile, opts.nthread);
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
    fprintf(stderr, "\t\t-i \t\tAdditionally index file.\n");
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
    int  idx    = 0;      //Index file flag
    int  nthread = 1;
    int  minargs = 2;
    while (( elem = getopt(argc, argv, "o:t:l:hdi") ) >= 0) {
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
            case 'i':
                idx = 1;
                continue;
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
        fprintf(stderr, "[sqz]: %s, file not found.\n\n", argv[argc - 1]);
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
    opts->idx = idx;
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
            if ( sqz_cmpr(opts) ) goto exit;
            break;
        case 2:
            if ( sqz_dcmp(opts)) goto exit;
            break;
    }
    ret = 0;
    exit:
        return ret;
}
