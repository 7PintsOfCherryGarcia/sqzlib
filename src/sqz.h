

char sqz_compress(const char *filename, const char *outname);


char sqz_squeezefastq(sqzfastx_t *sqz, FILE *ofp);


char sqz_spreadfasta(FILE *ifp, FILE *ofp);


char sqz_squeezefasta(sqzfastx_t *sqz, FILE *ofp);


char sqz_decompress(const char *filename, const char *outname);


char sqz_spreadfastq(FILE *ifp, FILE *ofp);


char sqz_ropts(int argc, char **argv);


void sqz_usage();


char *sqz_basename(char *namestr);


char sqz_readblksize(sqzblock_t *blk, FILE *fp);
