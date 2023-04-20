#define LOAD_SIZE 4UL*1024UL*1024UL
#define NAME_SIZE 1UL*1024UL*1024UL
#define B64       8U
#define B32       4U
#define HEADLEN   8U
#define NEND      128U
#define MAGIC     151324677U

#define NBLK      255U
#define QBLK       63U
#define EBLK        0U

#define FQH       '@'
#define FAH       '>'
#define FQS       '+'
#define NL        '\n'

#define sqz_fileformat(f) ((f)->fmt >> 3)
#define sqz_isfastq(f) (((f)->fmt & 3) == 2 ? 1 : 0)

typedef struct {
    char    *n;
    uint32_t nlen;
    uint8_t *s;
    uint8_t *q;
    uint64_t l;
    uint64_t maxl;
} sqzseq_t;

typedef struct {
    void     *data;
    uint64_t  size;
    uint64_t  pos;
} sqzbuff_t;

typedef struct {
    //Encoding
    sqzbuff_t  *blkbuff;
    //Compression
    sqzbuff_t  *cmpbuff;
    //Generic
    uint32_t   n;        //Generic counter
    uint64_t   datasize; //How much of current sequence had been decoded
    uint8_t    newblk;   //flag
} sqzblock_t;

typedef struct {
    //flags
    uint8_t     datflag;
    uint8_t     cmpflag;
    //data members
    uint64_t    size;
    uint64_t    offset;
    uint8_t     *seq;
    uint8_t     *qlt;
    sqzbuff_t   *namebuffer;
    sqzbuff_t   *readbuffer;
    sqzblock_t  *blk;
    //return members
    uint32_t    n;
    uint64_t    bases;
    uint32_t    blks;
    //Last sequence loaded
    sqzseq_t    *lastseq;
    sqzbuff_t   *lseqbuff;
    uint8_t     lseqflag; //Sequece has not completely been read into a buffer flag
    //Partially decoded sequences
    uint64_t    rem;     //Length of sequence remaining to be read
    uint64_t    prevlen; //Size of sequence currently being read
} sqzfastx_t;

typedef struct sqzFile_s {
    const char  *name;
    void        *gzfp;
    sqzfastx_t  *sqz;
    uint64_t    size;
    uint8_t     ff;
    uint8_t     fmt;
    uint64_t    filepos;
} *sqzFile;


uint8_t    sqz_gzopen(const char *filename, sqzFile sqzfp, const char *mode);
void       sqz_gzclose(sqzFile sqzfp);
int32_t    sqz_gzread(sqzFile file, void *buff, uint32_t len);
void       sqz_gzrewind(sqzFile sqzfp);
uint64_t   sqz_gztell(sqzFile sqzfp);
void       sqz_gzseek(sqzFile sqzfp, uint64_t OFS, uint8_t WH);
void       sqz_getformat(sqzFile sqzfp);
sqzfastx_t *sqz_fastxinit(uint8_t fmt, uint64_t size);
void       sqz_fastxkill(sqzfastx_t *sqz);
void       sqz_blkkill(sqzblock_t *blk);
uint64_t   sqz_inflate(sqzblock_t *blk);
size_t     sqz_deflate(sqzblock_t *blk, int level);
int64_t    sqz_zstdcompress(sqzblock_t *blk, int level);
uint64_t   sqz_zstddecompress(sqzblock_t *blk);
uint8_t    sqz_blkrealloc(sqzblock_t *blk, uint64_t newsize);
sqzbuff_t  *sqz_buffrealloc(sqzbuff_t *buff,  uint64_t size);
sqzseq_t   *sqz_seqrealloc(sqzseq_t *seq, uint64_t newsize);
