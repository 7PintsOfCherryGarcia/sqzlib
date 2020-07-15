#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread)

#define LOAD_SIZE 8*1024*1024
#define CHUNK 4*16384
#define TWO_BIT_MASK (3)
#define NEND (128)

//Table to change "ACGT" to 0123 else to 4
unsigned char seq_nt4_table[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


//Table to change 01234 to ACGTN
unsigned char seq_dec_table[128] = {
    'A', 'C','G', 'T',  'N', 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


//squeezma main struct
typedef struct {
    uint8_t *seqbuffer;
    size_t  seqlen;
    uint8_t *namebuffer;
    size_t  namelen;
    uint8_t *qualbuffer;
    size_t  quallen;
    char fmt;
    size_t n;
    char endflag;        //Indicate a sequece has not completely been read into a buffer
    size_t toread;       //Size of sequence still needed to be read
    size_t prevlen;      //Size of sequence currently being read
} sqzfastx_t;


typedef struct {
    void *codebuff;
    size_t codesize;
} sqzcodeblock_t;


typedef struct {
    void *cmpbuff;
    size_t cmpsize;
    FILE *outfile;
} sqzcmpblock_t;


//djb2 http://www.cse.yorku.ca/~oz/hash.html
uint64_t djb2(char *str);


uint64_t bit2encode(const unsigned char *str, uint32_t strlen);



/*
Returns number of bases decoded
*/
unsigned char bit2decode(const uint64_t *mer, char *decoded, uint32_t len);


//Comprsion functions shamelessly copypasata from https://zlib.net/zlib_how.html
/*
Compresses an array of data with zlibs deflate into the provided dest array
returns the number of compressed bytes
*/
size_t zlibsqueezze(void *nuts, size_t nutlength, uint8_t *dest, size_t destlen,int level);


/*
Decompress a zlib data stream stored in cmpbuff into dest
*/
size_t zlibpunch(void *cmpbuff, size_t cmplen, uint8_t *dest, size_t destlen);


/*
Returns pointer to first non ACGTacgt base in *strseq, If no nonACGTacgt bases are present
in the stirng, findn returns a pointer to the teminating null byte of strseq like in
strchrnul() in GNU
*/
const unsigned char *findn(const unsigned char *strseq);

/*
    Sequence compression loop
    This function takes a sequence and encodes it via 2 bit enconding for bases
    AaCcGgTt and runlength encoding for N. The following format is used to store the
    encoded sequence:
        [slen|32bit|][blen|32bit]*[bblock|n*64bit](n*[nblock|n*8bit])
    slen   - sequence length
    blen   - block length: Sequence length up to an N bases
    bblock - 64 bit integers encoding up to 32 base kmers each for a total of blen bases
    nblock - 7 bit integers encoding up to 127 consecutive Ns. First bit indicates if
             current nblock is the last one before next blen,bblock or end of sequence

    Encoding algorithm:
        1.Read DNA sequence
        2.Write sequence length to [slen]
        While there is sequence to compress:
            3.Find first occurance of non ACGTacgt base
            If a non ACGTacgt base is found
                4.Compute sequence length up to found non ACGTacgt base
                5.Compute number of consecutive Ns
                6.Write sequence length up to found N to [blen]
                Loop over sequence (up to found N in chunks of 32 bases)
                    7.Encode 32mer or whatever is left of sequence into [bblock]
                Loop over number of Ns in chunks of 127
                    If no more Ns left
                        8.Encode 127 Ns or number of Ns left in [nblock]
                        9.Set bit 7 of currnt [nblock]
                    If more Ns need to be coded
                        8.Encode 127Ns as byte=127 in [nblock]
            10.Compute how much sequence there is left
            If there is sequence to be encoded (trailing bases after last occurance of N)
                11.Encode sequence as in 6
*/
size_t loopencode(const unsigned char *str, uint32_t strlen, uint8_t *cmpbuff);


unsigned char writens(unsigned char numn, char *decoded);


size_t loopdecode(const uint8_t *buff);


int sqz_encodencompress(sqzfastx_t *sqz, size_t size);


char sqz_getformat(const char *filename);


sqzfastx_t *sqz_fastxinit(const char *filename, size_t buffersize);


void sqz_kill(sqzfastx_t *sqz);


char sqz_loadfasta(sqzfastx_t *sqz, kseq_t *seq);


char sqz_loadfastq(sqzfastx_t *sqz, kseq_t *seq);
