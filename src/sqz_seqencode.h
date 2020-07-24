#define TWO_BIT_MASK (3)
#define NEND (128)


unsigned char bit2decode(const uint64_t *mer, char *decoded, uint32_t len);


uint64_t bit2encode(const unsigned char *str, uint32_t strlen);


size_t sqz_seqencode(const unsigned char *str, uint32_t strlen, sqzcodeblock_t *codeblk);


const unsigned char *sqz_findn(const unsigned char *strseq);


size_t sqz_seqdecode(const uint8_t *buff);


unsigned char sqz_writens(unsigned char numn, char *decoded);
