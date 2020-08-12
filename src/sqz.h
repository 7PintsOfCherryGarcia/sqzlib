#include <stdio.h>


char sqz_compress(const char *filename, const char *outname);


char sqz_squeezefastq(sqzfastx_t *sqz, FILE *ofp);


char sqz_decompress(const char *filename, const char *outname);


char sqz_ropts(int argc, char **argv);


void sqz_usage();
