
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

typedef struct {
    char    *n;
    uint32_t nlen;
    uint8_t *s;
    uint8_t *q;
    uint64_t l;
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
    uint64_t   namepos; //Position within namebuffer TODO Rethink this memeber
    uint64_t   prevlen; //How much of current sequence had been decoded
    uint8_t    newblk;  //flag
} sqzblock_t;

/*
  "sqzfastx_t"
  libsqueezma main data loading structure. Defines the buffers and flags for
  reading sequencing data into.
*/
typedef struct {
    //flags
    uint8_t     datflag;
    char        endflag; //Sequece has not completely been read into a buffer flag
    char        cmpflag;
    //data members
    uint64_t    size; //Amount of sequence read.
    uint64_t    offset;
    uint8_t     *seq;
    uint8_t     *qlt;
    sqzbuff_t   *namebuffer;
    uint8_t     *readbuffer;
    uint64_t    namesize;
    uint64_t    namepos;
    sqzblock_t  *blk;
    //return members
    uint64_t    n;
    uint64_t    bases;
    uint64_t    blks;
    //Last sequence loaded
    sqzseq_t    *lastseq;
    sqzbuff_t   *lseqbuff;
    //Partially decoded sequences
    uint64_t    rem;     //Length of sequence remaining to be read
    uint64_t    prevlen; //Size of sequence currently being read
} sqzfastx_t;

typedef struct sqzFile_s {
    const char  *namestr;
    void        *gzfp;
    sqzfastx_t  *sqz;
    uint64_t    size;
    uint8_t     ff;
    uint8_t     fmt;
    uint8_t     libfmt;
    uint64_t    filepos;
} *sqzFile;


uint8_t  sqz_gzopen(const char *filename, sqzFile sqzfp, const char *mode);
void sqz_gzclose(sqzFile sqzfp);
int32_t  sqz_gzread(sqzFile file, void *buff, uint32_t len);
void     sqz_gzrewind(sqzFile sqzfp);
uint64_t sqz_gztell(sqzFile sqzfp);
void     sqz_gzseek(sqzFile sqzfp, uint64_t OFS, uint8_t WH);

void sqz_getformat(sqzFile sqzfp);


sqzfastx_t *sqz_fastxinit(uint8_t fmt, uint64_t size);
sqzblock_t *sqz_sqzblkinit(uint64_t size);
void sqz_fastxkill(sqzfastx_t *sqz);
void sqz_blkkill(sqzblock_t *blk);
uint64_t sqz_inflate(sqzblock_t *blk);
uint64_t sqz_fastXdecode(sqzblock_t *blk,
                         uint8_t *buff,
                         uint64_t buffsize,
                         uint8_t fqflag);
size_t sqz_deflate(sqzblock_t *blk, int level);
int64_t sqz_zstdcompress(sqzblock_t *blk, int level);
uint64_t sqz_zstddecompress(sqzblock_t *blk);
uint8_t sqz_sqzblkrealloc(sqzblock_t *blk, uint64_t newsize);
sqzbuff_t *sqz_buffrealloc(sqzbuff_t *buff,  uint64_t size);

