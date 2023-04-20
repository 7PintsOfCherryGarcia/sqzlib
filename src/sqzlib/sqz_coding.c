#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "sqz_data.h"
#include "klib/kseq.h"
KSEQ_INIT(sqzFile, sqz_gzread)

#define TWO_BIT_MASK 3

const uint8_t nblk  = 255U;
const uint8_t qblk  =  63U;
const uint8_t eblk  =   0U;


//Table to change "ACGT" to 0123 else to 4
const unsigned char seq_nt4_tableSQZ[128] = {
    128, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
      4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
      4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
      4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
      4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
      4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
      4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
      4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

//Table to change 01234 to ACGTN
const unsigned char seq_dec_tableSQZ[128] = {
    'A', 'C','G', 'T',  'N', 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

//Table to change 0,1,2,3,4,5,6,7 to quality values
const unsigned char qual_val_table[8] = {33,39,48,55,60,66,70,73};

static inline const uint8_t *sqz_findn(const uint8_t *seq)
{
    while (seq_nt4_tableSQZ[*seq] < 4)
        seq++;
    return seq;
}

static inline uint64_t sqz_bit2encode(const uint8_t *seq, uint8_t seqlen)
{
    uint64_t result = 0;
    for (uint64_t i = 0; i < seqlen; i++)
        result = (result << 2) | seq_nt4_tableSQZ[seq[i]];
    return result;
}

static void sqz_go2namen(sqzbuff_t *names, sqzseq_t  *seq, uint32_t n)
{
    uint32_t idx = 0;
    char *namelist = (char *)names->data;
    while (idx < n) {
        namelist += strlen(namelist) + 1;
        idx++;
    }
    seq->nlen = strlen(namelist);
    memcpy(seq->n, namelist, seq->nlen + 1);
}

static uint64_t sqz_getlength(sqzbuff_t *buff)
{
    uint64_t l =  *(uint64_t *)( (uint8_t *)buff->data + buff->pos );
    buff->pos += B64;
    return l;
}

static uint64_t sqz_getblklen(sqzbuff_t *buff)
{
    uint8_t *data = (uint8_t *)buff->data;
    uint64_t blklen = *(uint64_t *)(data + buff->pos);
    buff->pos += B64;
    return blklen;
}

static uint8_t sqz_fasta2buff(sqzseq_t *seq, sqzbuff_t *buff)
{
    uint64_t pos  = buff->pos;
    uint64_t size = buff->size;
    uint8_t *data = buff->data;
    if ( (seq->l + 3 + seq->nlen) > (size - pos) ) {
        buff = sqz_buffrealloc(buff, pos + seq->l + 3 + seq->nlen);
        data = buff->data;
        if (!data) return 1;
    }
    data[pos++] = FAH;
    memcpy(data + pos, seq->n, seq->nlen);
    pos += seq->nlen;
    data[pos++] = NL;
    memcpy(data + pos, seq->s, seq->l);
    pos += seq->l;
    data[pos++] = NL;
    buff->pos = pos;
    return 0;
}

static uint64_t sqz_blkcode(uint64_t *buff, const uint8_t *seq, uint64_t len)
{
    if (!len) return 0;
    uint64_t nblks = ( (len / 32) + ( (len % 32) > 0 ) ) - 1;
    uint64_t code = 0;
    uint64_t buffpos;
    uint64_t nbytes = 0;
    uint8_t lenleft = 0;
    //1 less than the total number of iterations needed.
    for (buffpos = 0; buffpos < nblks; buffpos++) {
        //Pack 32 mer into 64 bit int
        code = sqz_bit2encode(seq, 32);
        memcpy(buff + buffpos, &code, B64);
        nbytes += B64;
        seq += 32;
    }
    lenleft = ( len - ( nblks * 32) );
    code = sqz_bit2encode(seq, lenleft);
    memcpy(buff + buffpos, &code, B64);
    nbytes += B64;
    return nbytes;
}

static uint64_t sqz_nblkcode(uint8_t *buff, uint64_t len)
{
    if (!len) return 0;
    //1 less than the total number of iterations needed.
    uint64_t nblks = ( (len / 127) + ( (len % 127) > 0 ) ) - 1;
    uint64_t buffpos = 0;
    for (buffpos = 0; buffpos < nblks; buffpos++)
        buff[buffpos] = 127; // 127 consecutive non ACGT bases
    //Last itteration needs special treatment
    //Set bit 7 if last N block before next base
    buff[buffpos++] = (len - (nblks * 127)) | NEND;
    return buffpos;
}

static uint8_t sqz_8binqual(uint8_t q)
{
    if (q >= 73) return 7<<5;
    if ( (q < 73) && (q >= 68) ) return 6<<5;
    if ( (q < 68) && (q >= 63) ) return 5<<5;
    if ( (q < 63) && (q >= 58) ) return 4<<5;
    if ( (q < 58) && (q >= 53) ) return 3<<5;
    if ( (q < 53) && (q >= 43) ) return 2<<5;
    if ( (q < 43) && (q >= 35) ) return 1<<5;
    if ( (q < 35) && (q > 0)   ) return 1;
    return 0;
}

static uint8_t sqz_1writebuff(sqzbuff_t *buff, uint8_t uint)
{
    memcpy((uint8_t *)buff->data + buff->pos, &uint, 1);
    buff->pos += 1;
    return 1;
}

static uint64_t sqz_qualencode(const uint8_t *qlt, uint64_t l, sqzbuff_t *buff)
{
    uint64_t pos = 0;
    uint8_t q       = sqz_8binqual(*qlt);
    uint8_t c       = 0;
    uint8_t code    = 0;
    while(l) {
        //If next qual val is the same
        if ( sqz_8binqual(*(qlt + 1)) == q) {
            //Increase counter
            c++;
            //If counter reached
            if (c == 31) {
                //Encode
                qlt++;
                l--;
                code = code | (q & 224);
                code = code | c;
                pos += sqz_1writebuff(buff, code);
                c = 0;
                code = 0;
                q = sqz_8binqual(*(qlt + 1));
            }
            qlt++;
            l--;
            continue;
        }
        //Encode
        code = code | (q & 224);
        code = code | c;
        pos += sqz_1writebuff(buff, code);
        qlt++;
        l--;
        q = sqz_8binqual(*qlt);
        c = 0;
        code = 0;
    }
    buff->pos += pos;
    return pos;
}

static uint64_t sqz_seqencode(const uint8_t *seq, uint64_t l, sqzbuff_t *buff)
{
    uint64_t pos = 0;
    const uint8_t *lstop  = seq;
    const uint8_t *nptr   = seq + l;
	  const uint8_t *npos   = NULL;
    uint64_t nn     = 0;
	  uint64_t blen   = 0;
    uint8_t *inbuff = (uint8_t *)buff->data + buff->pos;
    do {
        //Get position of first non ACGTacgt base
        npos = sqz_findn(lstop);
        if (*npos) {
            //Determine block length up to first non ACGT base
            blen = npos - lstop;
            //Determine number of consecutive Ns until next base
            nn = 0;
	          while ( seq_nt4_tableSQZ[*npos] == 4) {
	              nn++;
		            npos++;
	          }
            //Write block length
            memcpy(inbuff + pos, &blen, B64);
            pos += B64;
            //Encode sequence
            pos += sqz_blkcode((uint64_t *)(inbuff + pos), lstop, blen);
            //Trace next base after sequence of Ns
            lstop = npos;
            //Indicate that a non-ACGT block follows
            memcpy(inbuff + pos, &nblk, 1);
            pos++;
            //Encode non-ACGT block
            pos += sqz_nblkcode(inbuff + pos, nn);

        }
    } while (*npos);
    //Detect and encode trailing bases
    blen = nptr - lstop;
    if (blen) {
        memcpy(inbuff + pos, &blen, B64);
        pos += B64;
        pos += sqz_blkcode((uint64_t *)(inbuff + pos), lstop, blen);
    }
    buff->pos += pos;
    return pos;
}

static uint64_t sqz_appendbuff(sqzbuff_t *dest, sqzbuff_t *src)
{
    memcpy((uint8_t *)dest->data + dest->pos, src->data, src->pos);
    dest->pos += src->pos;
    return src->pos;
}

static uint64_t sqz_64writebuff(sqzbuff_t *buff, uint64_t uint)
{
    memcpy((uint8_t *)buff->data + buff->pos, &uint, B64);
    buff->pos += B64;
    return B64;
}

static uint8_t sqz_fastaheadblk(sqzfastx_t *sqz)
{
    sqzblock_t *blk = sqz->blk;
    sqzbuff_t  *buff = blk->blkbuff;
    uint8_t *seqb    = sqz->seq;
    uint64_t sqzsize   = sqz->offset;
    uint8_t *seq = NULL;
    uint64_t seqlen = 0;
    uint64_t k = 0;
    while ( k < sqzsize ) {
        seqlen = *(uint64_t *)( seqb + k );
        k += B64;
        sqz_64writebuff(buff, seqlen);
        seq = seqb + k;
        sqz_seqencode(seq, seqlen, buff);
        sqz_1writebuff(buff, eblk);
        k += seqlen + 1;
    }
    if (sqz->lseqflag) {
        sqzseq_t *sqzseq = sqz->lastseq;
        sqz_64writebuff(buff, sqzseq->l);
        sqz_seqencode((uint8_t *)sqzseq->s, sqzseq->l, sqz->lseqbuff);
        if ( (sqz->lseqbuff->pos + buff->pos) > buff->size)
            if ( sqz_blkrealloc(blk, sqz->lseqbuff->pos + buff->pos + 1) )
                return 1;
        sqz_appendbuff(buff, sqz->lseqbuff);
        sqz_1writebuff(buff, eblk);
    }
    if ( (buff->pos + sqz->namebuffer->pos) >= buff->size )
        if ( sqz_blkrealloc(blk, buff->pos + sqz->namebuffer->pos + B64) )
            return 1;
    sqz_appendbuff(buff, sqz->namebuffer);
    sqz_64writebuff(buff, sqz->namebuffer->pos);
    sqz->cmpflag = 1;
    return 1;
}

static uint8_t sqz_fastqheadblk(sqzfastx_t *sqz)
{
    sqzblock_t *blk = sqz->blk;
    sqzbuff_t  *buff  = blk->blkbuff;
    uint64_t  blkpos     = 0;
    uint8_t  *seqb       = sqz->seq;
    uint8_t  *qltb       = sqz->qlt;
    uint64_t  sqzsize    = sqz->offset;
    uint8_t  *seq     = NULL;
    uint8_t  *qlt     = NULL;
    uint64_t  seqlen  = 0;
    uint64_t  k       = 0;
    while ( k < sqzsize ) {
        seqlen = *(uint64_t *)( seqb + k );
        k += B64;
        blkpos += sqz_64writebuff(buff, seqlen);
        seq  = seqb + k;
        qlt  = qltb + k;
        blkpos += sqz_seqencode(seq, seqlen, buff);
        blkpos += sqz_1writebuff(buff, qblk);
        blkpos += sqz_qualencode(qlt, seqlen, buff);
        k += seqlen + 1;
    }
    if (sqz->lseqflag) {
        sqzseq_t *sqzseq = sqz->lastseq;
        blkpos += sqz_64writebuff(buff, sqzseq->l);
        sqz_seqencode((uint8_t *)sqzseq->s, sqzseq->l, sqz->lseqbuff);
        sqz_1writebuff(sqz->lseqbuff, qblk);
        sqz_qualencode((uint8_t *)sqzseq->q, sqzseq->l, sqz->lseqbuff);
        if ( (sqz->lseqbuff->pos + blkpos) > buff->size) {
            buff = sqz_buffrealloc(buff, sqz->lseqbuff->pos + blkpos + 1);
            blk->blkbuff = buff;
        }
        blkpos += sqz_appendbuff(buff, sqz->lseqbuff);

    }
    if ( (buff->pos + sqz->namebuffer->pos) >= buff->size ) {
        buff = sqz_buffrealloc(buff, buff->pos + sqz->namebuffer->pos + B64);
        blk->blkbuff = buff;
    }
    blkpos += sqz_appendbuff(buff, sqz->namebuffer);
    blkpos += sqz_64writebuff(buff, sqz->namebuffer->pos);
    sqz->cmpflag = 1;
    return 1;
}

char sqz_fastXencode(sqzfastx_t *sqz, uint8_t fqflag)
{
    if (fqflag) return sqz_fastqheadblk(sqz);
    return sqz_fastaheadblk(sqz);
}

/*
################################################################################
  Decoding functions
################################################################################
*/
static inline uint64_t sqz_blksize(uint64_t blklen)
{
    return ( (blklen / 32) + ((blklen % 32) > 0) ) * B64;
}

static uint64_t sqz_qdecode(const uint8_t  *codebuff,
                            uint8_t  *outbuff,
                            uint64_t rquals,
                            uint64_t offset)
{
    uint64_t pos = 0;
    uint64_t nbyte = 0;
    uint8_t  q;
    uint8_t  prevqn;
    uint8_t  i = 0;
    uint8_t  j;
    //Go to byte of current offset
    while (pos < offset) {
        pos += 1 + (*(codebuff + nbyte) & 31);
        nbyte++;
    }
    prevqn = pos - offset;
    if (prevqn) {
        q = (*(codebuff - 1) & 224) >> 5;
        for (i = 0; i < prevqn; i++)
            *outbuff++ = qual_val_table[q];
        rquals -= prevqn;
    }
    while (rquals) {
        i = 1 + (*(codebuff + nbyte) & 31);
        q = (*(codebuff + nbyte) & 224) >> 5;
        nbyte++;
        i = i > rquals ? rquals : i;
        rquals -= i;
        for (j = 0; j < i; j++)
            *outbuff++ = qual_val_table[q];
    }
    return nbyte;
}

/*
  Counts number of Ns in an N block
  Returns number of Ns
  Modifies pos with number of bytes decoded
*/
static uint64_t sqz_countnblk(sqzbuff_t *buff)
{
    uint8_t *codebuff = (uint8_t *)buff->data + buff->pos;
    uint64_t count = 0;
    uint64_t i = 0;
    while ( !(codebuff[i++] & 128) ) {
        count += 127;
    }
    count += codebuff[i - 1] & 127;
    buff->pos += i;
    return count;
}

static void sqz_writens(uint64_t  numn, uint8_t *decoded)
{
    while (numn-- > 0)
        decoded[numn] = 'N';
}

static void sqz_merdecode(uint64_t mer,
                          uint8_t *outbuff,
                          uint8_t  discard,
                          uint8_t  nbases)
{
    mer >>= (discard * 2);
    uint8_t byte;
    while(nbases--) {
        byte = mer & TWO_BIT_MASK;
        outbuff[nbases] = seq_dec_tableSQZ[byte];
        mer >>= 2;
    }
}

/*
Decodes blklen encoded bases
Returns number of bytes decoded
Updates wbytes with number of bytes written
*/
static uint64_t sqz_blkdecode(sqzbuff_t *buff, uint8_t *outbuff, uint64_t blklen)
{
    if (!blklen) return 0;
    uint8_t *codebuff = (uint8_t *)buff->data + buff->pos;
    uint64_t codepos = 0;
    //1 less than the total number of iterations needed
    uint64_t blknum  = ( (blklen / 32) + ( (blklen % 32) > 0) ) - 1;
    uint64_t outpos  = 0;
    uint64_t blkidx  = 0;
    uint64_t blkleft  = (blklen % 32) == 0 ? 32 : blklen % 32;
    uint64_t *mer;
    mer = (uint64_t *)codebuff;
    for (blkidx = 0; blkidx < blknum; blkidx++) {
        sqz_merdecode( *(mer + blkidx), outbuff + outpos, 0, 32);
        codepos += B64;
        outpos += 32;
    }
    sqz_merdecode( *(mer + blkidx), outbuff + outpos, 0, blkleft);
    codepos += B64;
    outpos += blkleft;
    buff->pos += codepos;
    return outpos;
}

static uint8_t sqz_seqdecode(sqzbuff_t *buff, sqzseq_t  *sqzseq, char qflag)
{
    uint64_t l = sqz_getlength(buff);
    if (l > sqzseq->maxl) {
        sqzseq = sqz_seqrealloc(sqzseq, l);
        if (!sqzseq) return 1;
    }
    uint64_t blklen = 0;
    uint64_t bases  = 0;
    uint64_t nnum = 0;
    while (bases < l) {
        blklen = sqz_getblklen(buff);
        bases += sqz_blkdecode(buff,
                               sqzseq->s + bases,
                               blklen);
        if (bases < l) {
            buff->pos++; //nflag byte
            nnum = sqz_countnblk(buff);
            sqz_writens(nnum, sqzseq->s + bases);
            bases += nnum;
            if (bases == l) buff->pos++;  //Sequence ends in non-ACGT block
            continue;
        }
        //Flag byte (q or end os seq)
        buff->pos++;
    }
    //Decode quality strings for fastq data
    //TODO deal with this later
    if (qflag)
        sqz_qdecode((uint8_t *)buff->data, sqzseq->q, l, 0);
    if (bases != l) return 1;
    sqzseq->l = bases;
    return 0;
}

static uint64_t sqz_fastXdecode(sqzfastx_t *sqz, sqzbuff_t *outbuff, uint8_t fqflag)
{
    sqzbuff_t *blkbuff = sqz->blk->blkbuff;
    uint64_t  datasize  = sqz->blk->datasize;
    sqzbuff_t *namebuff = sqz->namebuffer;
    sqzseq_t  *sqzseq   = sqz->lastseq;
    uint32_t  n         = 0;
    while ( blkbuff->pos < datasize ) {
        sqz_go2namen(namebuff, sqzseq, n);
        if ( sqz_seqdecode(blkbuff, sqzseq, fqflag) )
            goto exit;
        if ( sqz_fasta2buff(sqzseq, outbuff) ) {
            outbuff->pos = 0;
            n = 0;
            goto exit;
        }
        n++;
    }
    exit:
        sqz->n = n;
        return outbuff->pos;
}


uint64_t sqz_decode(sqzFile sqzfp)
{
    sqzbuff_t *buff = sqzfp->sqz->readbuffer;
    return sqz_fastXdecode(sqzfp->sqz, buff, sqz_isfastq(sqzfp));
}

//WARNINGS AVIDING
void dummy_fun(kseq_t *seq)
{
    kseq_read(seq);
    kseq_destroy(seq);
    kseq_init(NULL);
}
