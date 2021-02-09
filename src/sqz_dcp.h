char sqz_getformat(const char *filename);


char sqz_kseqinit(sqzfastx_t *sqz);


uint64_t sqz_filesize(FILE *fp);


static void sqz_fastxreset(sqzfastx_t *sqz);


static void sqz_blkreset(sqzblock_t *blk);
