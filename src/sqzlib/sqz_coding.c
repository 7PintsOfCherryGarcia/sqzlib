#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>

#include "sqz_data.h"

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


static uint64_t sqz_qualencode(const uint8_t *qlt, uint64_t l, uint8_t *blkbuff)
{
    uint64_t blkpos = 0;
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
                memcpy(blkbuff + blkpos, &code, 1);
                blkpos++;
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
        memcpy(blkbuff + blkpos, &code, 1);
        blkpos++;
        qlt++;
        l--;
        q = sqz_8binqual(*qlt);
        c = 0;
        code = 0;
    }
    return blkpos;
}


static uint64_t sqz_seqencode(const uint8_t *seq, uint64_t l, uint8_t *blkbuff)
{
    uint64_t blkpos = 0;
    const uint8_t *lstop  = seq;
    const uint8_t *nptr   = seq + l;
	  const uint8_t *npos   = NULL;
    uint64_t nn     = 0;
	  uint64_t blen   = 0;
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
            memcpy(blkbuff + blkpos, &blen, B64);
            blkpos += B64;
            //Encode sequence
            blkpos += sqz_blkcode((uint64_t *)(blkbuff + blkpos), lstop, blen);
            //Trace next base after sequence of Ns
            lstop = npos;
            //Indicate that a non-ACGT block follows
            memcpy(blkbuff + blkpos, &nblk, 1);
            blkpos++;
            //Encode non-ACGT block
            blkpos += sqz_nblkcode(blkbuff + blkpos, nn);

        }
    } while (*npos);
    //Detect and encode trailing bases
    blen = nptr - lstop;
    if (blen) {
        memcpy(blkbuff + blkpos, &blen, B64);
        blkpos += B64;
        blkpos += sqz_blkcode((uint64_t *)(blkbuff + blkpos), lstop, blen);
    }
    return blkpos;
}


static uint8_t sqz_fastqheadblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    uint8_t  *blkbuff    = blk->blkbuff;
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
        memcpy(blkbuff + blkpos, &seqlen, B64);
        blkpos += B64;
        seq  = seqb + k;
        qlt  = qltb + k;
        blkpos += sqz_seqencode(seq, seqlen, blkbuff + blkpos);
        memcpy(blkbuff + blkpos, &qblk, 1);
        blkpos++;
        blkpos += sqz_qualencode(qlt, seqlen, blkbuff + blkpos);
        k += seqlen + 1;
    }
    //Last sequence, loaded into kseq, but not copied into buffer
    if (sqz->endflag) {
        memcpy(blkbuff + blkpos, &sqz->prevlen, B64);
        blkpos += B64;
        blkpos += sqz_seqencode(sqz->pseq, sqz->prevlen, blkbuff + blkpos);
        memcpy(blkbuff + blkpos, &qblk, 1);
        blkpos++;
        blkpos += sqz_qualencode(sqz->pqlt, sqz->prevlen, blkbuff + blkpos);
    }
    //Copy name data making sure it fits
    uint64_t blksize = blk->blksize;
    if ( (blkpos + sqz->namepos) >= blksize ) {
        blksize += sqz->namepos + B64;
        blkbuff = realloc(blk->blkbuff, blksize);
        blk->blkbuff = blkbuff;
        blk->blksize = blksize;
        blk->cmpsize = blksize;
        blk->cmpbuff = realloc(blk->cmpbuff, blk->cmpsize);
    }
    memcpy(blkbuff + blkpos, sqz->namebuffer, sqz->namepos);
    blkpos += sqz->namepos;
    memcpy(blkbuff + blkpos, &(sqz->namepos), B64);
    blkpos += B64;
    sqz->cmpflag = 1;
    blk->blkpos = blkpos;
    return 1;
}


static uint8_t sqz_fastaheadblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    uint8_t *blkbuff   = blk->blkbuff;
    uint64_t blkpos    = blk->blkpos;
    uint8_t *seqb      = sqz->seq;
    uint64_t sqzsize   = sqz->offset;
    uint8_t *seq = NULL;
    uint64_t seqlen = 0;
    uint64_t k = 0;
    while ( k < sqzsize ) {
        seqlen = *(uint64_t *)( seqb + k );
        k += B64;
        memcpy(blkbuff + blkpos, &seqlen, B64);
        blkpos += B64;
        seq = seqb + k;
        blkpos += sqz_seqencode(seq, seqlen, blkbuff + blkpos);
        memcpy(blkbuff + blkpos, &eblk, 1);
        blkpos++;
        k += seqlen + 1;
    }
    //Last sequence, loaded into kseq, but not copied into buffer
    if (sqz->endflag) {
        memcpy(blkbuff + blkpos, &sqz->prevlen, B64);
        blkpos += B64;
        blkpos += sqz_seqencode(sqz->pseq, sqz->prevlen, blkbuff + blkpos);
        memcpy(blkbuff + blkpos, &eblk, 1);
        blkpos++;
    }
    //Copy name data making sure it fits
    uint64_t blksize = blk->blksize;
    if ( (blkpos + sqz->namepos) >= blksize ) {
        blksize += sqz->namepos + B64;
        blkbuff = realloc(blk->blkbuff, blksize);
        blk->blkbuff = blkbuff;
        blk->blksize = blksize;
        blk->cmpsize = blksize;
        blk->cmpbuff = realloc(blk->cmpbuff, blk->cmpsize);
    }
    memcpy(blkbuff + blkpos, sqz->namebuffer, sqz->namepos);
    blkpos += sqz->namepos;
    memcpy(blkbuff + blkpos, &(sqz->namepos), B64);
    blkpos += B64;
    sqz->cmpflag = 1;
    blk->blkpos = blkpos;
    return 1;
}


char sqz_fastXencode(sqzfastx_t *sqz, sqzblock_t *blk, uint8_t fqflag)
{
    if (fqflag) return sqz_fastqheadblk(sqz, blk);
    return sqz_fastaheadblk(sqz, blk);
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
static uint64_t sqz_countnblk(const uint8_t *buff, uint64_t *pos)
{
    uint64_t count = 0;
    uint64_t i = 0;
    while ( !(buff[i++] & 128) ) {
        count += 127;
    }
    count += buff[i - 1] & 127;
    *pos = i;
    return count;
}


static uint64_t sqz_qbytes(uint8_t *codebuff, uint64_t bases)
{
    uint64_t rbases = 0;
    uint64_t nbytes = 0;
    while (rbases < bases) {
        rbases += 1 + (*(codebuff++) & 31);
        nbytes++;
    }
    return nbytes;
}


static uint64_t sqz_gotoqblk(uint8_t *codebuff)
{
    uint64_t codepos = 0;
    uint64_t blklen;
    uint8_t  nflag;
    uint64_t bases = 0;
    uint64_t pos;
    while (1) {
        blklen = *(uint64_t *)(codebuff + codepos);
        codepos += B64;
        if (blklen)
            codepos += sqz_blksize(blklen);
        bases += blklen;
        nflag = *(codebuff + codepos);
        codepos++;
        switch (nflag) {
        case NBLK:
            bases += sqz_countnblk(codebuff + codepos, &pos);
            codepos += pos;
            continue;
        case QBLK:
            goto exit;
        case EBLK:
            goto exit;
        }
    }
    exit:
        return codepos;
}


static uint64_t sqz_gotoblk(uint8_t  *codebuff,
                            uint64_t *retpos,
                            uint64_t  offset)
{
    uint64_t bases    = 0;
    uint64_t rbases = 0;
    uint64_t codepos  = 0;
    uint64_t rpos     = 0;
    uint64_t npos     = 0;
    uint8_t  nflag    = 0;
    uint64_t blklen;
    while (offset >= bases ) {
        //Store how many bases have been counted
        rbases = bases;
        rpos = codepos;
        //Count block
        blklen = *(uint64_t *)(codebuff + codepos);
        codepos += B64;
        bases += blklen;
        if (blklen)
            codepos += sqz_blksize(blklen);
        //Count N block
        nflag = *(codebuff + codepos);
        codepos++;
        switch (nflag) {
        case NBLK:
            bases += sqz_countnblk(codebuff + codepos, &npos);
            codepos += npos;
            continue;
        case QBLK:
            break;
        case EBLK:
            break;
        }
    }
    *retpos += rpos;
    return rbases;
}


static void sqz_writens(uint64_t  numn, uint8_t *decoded)
{
    while (numn-- > 0)
        decoded[numn] = 'N';
}


static uint64_t sqz_ndecode(uint8_t *codebuff,
                            uint8_t *outbuff,
                            uint64_t offset,
                            uint64_t len,
                            uint64_t *retpos,
                            uint64_t *blkoffset)
{
    uint64_t codepos = 0;
    uint64_t nlength = 0;
    uint64_t nbases  = 0;
    //Compute N blk length and number of bytes it needed for encoding
    nlength = sqz_countnblk(codebuff, &codepos);
    //Determine how many Ns will fit in buffer
    nbases = len < ( nlength - offset ) ? len : ( nlength - offset );
    sqz_writens(nbases, outbuff);
    //Move to next block if no more Ns are left to decode
    if ( !(nlength - offset - nbases)  ) {
        *retpos += codepos;
        *blkoffset = 0;
    }
    return nbases;
}


static void sqz_merdecode(uint64_t mer,
                          uint8_t *outbuff,
                          uint8_t  startpos,
                          uint8_t  nbases)
{
    uint8_t discard = 32U - startpos - nbases;
    mer >>= (discard * 2);
    uint8_t byte;
    while(nbases--) {
        byte = mer & TWO_BIT_MASK;
        outbuff[nbases] = seq_dec_tableSQZ[byte];
        mer >>= 2;
    }
}


static void sqz_mblk(const uint64_t *codebuff,
                     uint8_t *outbuff,
                     uint64_t offset,
                     uint64_t len)
{
    uint64_t codepos = 0;
    uint64_t outpos  = 0;
    uint64_t blknum  = 0;
    //Compute start position within 32mer
    uint8_t  startpos = offset % 32;
    //Number of bases to decode from this first 32mer
    uint8_t  nbases   = len < (uint64_t)(32 - startpos) ?
                        (uint8_t)len : (32 - startpos);
    if (nbases) blknum++;
    //Compute how many blks encode the number of bases requested
    blknum += ( (len - nbases) / 32) + ( ( (len - nbases)  % 32) > 0 ) - 1;
    //Move to corresponding 32mer
    const uint64_t *mer = codebuff + (offset / 32);
    uint64_t blkidx;
    for (blkidx = 0; blkidx < blknum; blkidx++) {
        sqz_merdecode(*(mer + blkidx), outbuff + outpos, startpos, nbases);
        codepos += B64;
        outpos += nbases;
        startpos = 0;
        nbases = 32;
    }
    len -= outpos;
    //How large is this last mer?
    sqz_merdecode(*(mer + blkidx), outbuff + outpos, startpos, len);
    codepos += B64;
    outpos += len;
}


/*
  Returns size of a sequence code block
*/
static uint64_t sqz_codeblksize(uint8_t *codebuff)
{
    uint64_t bases    = 0;
    uint64_t codepos  = 0;
    uint64_t blkpos2  = 0;
    uint64_t blklen   = 0;
    uint8_t  nflag    = 0;
    uint64_t seqlen = *(uint64_t *)codebuff;
    codepos += B64;
    while (bases != seqlen) {
        blklen = *(uint64_t *)(codebuff + codepos);
        codepos += B64;
        if (blklen) {
            bases += blklen;
            codepos += sqz_blksize(blklen);
        }
        //Read N flag
        nflag = *(codebuff + codepos);
        codepos++;
        //Check and account for Ns
        switch (nflag) {
        case NBLK:
            bases += sqz_countnblk(codebuff + codepos, &blkpos2);
            codepos += blkpos2;
            if (bases == seqlen) codepos++;  //Sequence ends in non-ACGT block
            continue;
        case QBLK:
            codepos += sqz_qbytes(codebuff + codepos, bases);
            continue;
        case EBLK:
            break;
        }
    }
    return codepos;
}


/*
Decodes blklen encoded bases
Returns number of bytes decoded
Updates wbytes with number of bytes written
*/
static uint64_t sqz_blkdecode(const uint8_t *codebuff,
                              uint8_t       *outbuff,
                              uint64_t      *wbytes,
                              uint64_t      blklen)
{
    if (!blklen) return 0;
    //1 less than the total number of iterations needed
    uint64_t blknum  = ( (blklen / 32) + ( (blklen % 32) > 0) ) - 1;
    uint64_t codepos = 0;
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
    *wbytes += outpos;
    return codepos;
}


static void sqz_pblk(uint8_t *codebuff,
                     uint64_t offset,
                     uint64_t len,
                     uint8_t *outbuff)
{
    uint64_t outpos = 0;
    uint64_t nnum = 0;
    uint64_t availbases  = len;
    uint64_t codepos   = 0;
    uint64_t blklen   = 0;
    uint64_t noffset  = 0;
    uint64_t blkbases = 0;
    //Go to corresponding block (block corresponding to current offset)
    blkbases = sqz_gotoblk(codebuff, &codepos, offset);
    offset -= blkbases;
    while (outpos < len) {
        //Determine seq blk, n blk, or quality block
        blklen = *(uint64_t *)(codebuff + codepos);
        //Offset within N block?
        if (offset >= blklen) {
            codepos += (B64 + sqz_blksize(blklen) + 1);
            noffset = offset - blklen;
            nnum = sqz_ndecode(codebuff + codepos,
                               outbuff + outpos,
                               noffset,
                               availbases,
                               &codepos,
                               &offset);
            outpos += nnum;
            availbases -= nnum;
            continue;
        }
        //How much of the block is available
        uint64_t blkavail = blklen - offset;
        uint64_t bases;
        bases =  blkavail < len - outpos ? blkavail : len - outpos;
        sqz_mblk( (uint64_t*)(codebuff + codepos + B64),
                 outbuff + outpos,
                 offset,
                 bases);
        outpos += bases;
        offset += bases;
        availbases -= bases;
    }
}


/*
  Given a partially decoded sequence, continue decoding from where it was left
  Retrns number of decoded bytes
*/
static uint64_t sqz_rdecode(uint8_t  *codebuff,
                            uint8_t  *outbuff,
                            uint64_t outsize,
                            uint64_t seqlen,
                            uint64_t decoded,
                            uint8_t  fqflag)
{
    uint64_t bleft     = 0;
    uint64_t seq2decode  = 0;
    uint64_t qlt2decode = 0;
    uint64_t qdecoded  = 0;
    seq2decode = seqlen < decoded ? 0 : seqlen - decoded;
    seq2decode = seq2decode < outsize ? seq2decode : outsize;
    qdecoded = seqlen < decoded ? decoded - seqlen - 3 : 0;
    //Compute how much buffer will be available
    bleft = outsize - seq2decode;
    if (seq2decode) {
        sqz_pblk(codebuff, decoded, seq2decode, outbuff);
        if (bleft < 3) goto exit;
    }
    //In case of fastq, we have to deal with qualities
    if (fqflag) {
        if (!qdecoded) {
            outbuff[seq2decode++] = NL;
            outbuff[seq2decode++] = '+';
            outbuff[seq2decode++] = NL;
            bleft -= 3;
        }
        //How many qualities can be decoded
        qlt2decode = seqlen - qdecoded < bleft ? (seqlen - qdecoded) : bleft;
        sqz_qdecode(codebuff + sqz_gotoqblk(codebuff),
                    outbuff + seq2decode,
                    qlt2decode,
                    qdecoded);
    }
    exit:
        return seq2decode + qlt2decode;
}


/*
  Decodes outsize bytes from blkbuff
  Returns number of bytes written
*/
static uint64_t sqz_pdecode(uint8_t  *codebuff,
                            uint8_t  *outbuff,
                            uint64_t outsize,
                            uint64_t seqlen,
                            uint8_t  fqflag)
{
    uint64_t outpos = 0;
    //Determine number of bases that fit in buffer
    uint64_t todecode = (seqlen < outsize) ? seqlen : outsize;
    sqz_pblk(codebuff, 0, todecode, outbuff + outpos);
    outpos += todecode;
    outsize -= todecode;
    //Decode quality if enough space
    if ( (outsize > 3) && fqflag ) {
        outbuff[outpos++] = NL;
        outbuff[outpos++] = '+';
        outbuff[outpos++] = NL;
        outsize -= 3;
        sqz_qdecode(codebuff + sqz_gotoqblk(codebuff),
                    outbuff + outpos,
                    outsize,
                    0);
        outpos += outsize;
    }
    return outpos;
}


/*
Loads length bases into decodebuffer
Returns number of bytes decoded
Updates wbytes with number of bytes written
*/
static uint64_t sqz_seqdecode(const uint8_t *codebuff,
                              uint8_t       *outbuff,
                              uint64_t       length,
                              char           qflag,
                              uint64_t      *wbytes)
{
    uint64_t seqlen     = length;
    uint64_t codepos    = 0;
    uint64_t blklen     = 0;
    uint64_t qltnum     = 0;
    uint64_t outpos     = 0;
    uint64_t nnum = 0;
    uint64_t nbytes;
    while (length > 0) {
        blklen = *(uint64_t *)(codebuff + codepos);
        qltnum += blklen;
        codepos += B64;
        codepos += sqz_blkdecode(codebuff + codepos,
                                 outbuff + outpos,
                                 &outpos,
                                 blklen);
        length -= blklen;
        if (length) {
            codepos++; //nflag byte
            nnum = sqz_countnblk(codebuff + codepos, &nbytes);
            qltnum += nnum;
            codepos += nbytes;
            sqz_writens(nnum, outbuff + outpos);
            outpos += nnum;
            length -= nnum;
            if (!length) codepos++;  //Sequence ends in non-ACGT block
            continue;
        }
        //Flag byte (q or end os seq)
        codepos++;
    }
    //Decode quality strings for fastq data
    if (qflag) {
        outbuff[outpos++] = NL;
        outbuff[outpos++] = '+';
        outbuff[outpos++] = NL;
        codepos += sqz_qdecode(codebuff + codepos, outbuff + outpos, seqlen, 0);
        outpos += seqlen;
    }
    outbuff[outpos++] = NL;
    *wbytes += outpos;
    return codepos;
}


uint64_t sqz_fastXdecode(sqzblock_t *blk,   //Data block
                         uint8_t *outbuff,  //Array to place decoded data
                         uint64_t outsize,  //Size of outbuff
                         char fqflag)       //0 - fasta, 1 - fastq
{
    uint8_t  *codebuff   = blk->blkbuff;
    uint64_t  codepos    = blk->blkpos;
    uint64_t  outpos     = 0;
    //Last 8 bytes store size of name buffer
    uint64_t  namesize  = *(uint64_t *)( codebuff + ( blk->blksize - B64 ) );
    //Size of sequence data in block
    uint64_t  datasize  = blk->blksize - B64 - namesize;
    char     *namebuff  = (char *)( codebuff + datasize );
    uint64_t  namepos   = blk->namepos;
    uint64_t  prevbytes = blk->prevlen;
    //Read length of first sequence
    uint64_t seqlen     = *(uint64_t *)( codebuff + codepos );;
    //Length of name of first sequence
    uint64_t namelen    = strlen(namebuff + namepos);
    //Extra bytes per seq
    uint8_t E = fqflag ? 6 : 3;
    uint64_t buffneed = ( seqlen * (1 + fqflag) )  + (E + namelen);
    //Check if there is sequence that was not completely decoded
    if (prevbytes) {
        uint64_t decoded = prevbytes - namelen - 2;
        outpos = sqz_rdecode(codebuff + codepos + B64,
                             outbuff,
                             outsize,
                             seqlen,
                             decoded,
                             fqflag);
        prevbytes += outpos;
        //Check that no more sequence needs decoding. (Move to next sequence)
        if (prevbytes < buffneed - 1)
            goto exit;
        outbuff[outpos++] = NL;
        outsize -= outpos;
        prevbytes = 0;
        //Update code buffer position to next sequence
        codepos += sqz_codeblksize(codebuff + codepos);
        //Get new sequence info
        //TODO: Overflow read when end of block has been reached
        namepos += namelen + 1;
        namelen  = strlen(namebuff + namepos);
        seqlen   = *(uint64_t *)( codebuff + codepos );
        buffneed = ( seqlen * (1 + fqflag) )  + (E + namelen);
    }
    //While there is data to decode
    while ( codepos < datasize ) {
        //Test if there is enough space for current sequence
        if ( outsize < buffneed ) {
            //Test if at least: ">"|"@" + name length + '\n' + 1 base fits
            if (outsize < namelen + 3) {
                //Can't decode this sequence, just return buffer
                prevbytes    = 0;
                blk->blkpos  = codepos;
                blk->namepos = namepos;
                goto exit;
            }
            outbuff[outpos++] = fqflag ? FQH : FAH;
            memcpy(outbuff + outpos, namebuff + namepos, namelen);
            outpos += namelen;
            outbuff[outpos++] = NL;
            outsize -= (2 + namelen);
            prevbytes = sqz_pdecode(codebuff + codepos + B64,
                                    outbuff + outpos,
                                    outsize,
                                    seqlen,
                                    fqflag);
            outpos += prevbytes;
            prevbytes += (2 + namelen);
            blk->blkpos  = codepos;
            blk->namepos = namepos;
            goto exit;
        }
        //Sequence fits and there is still space in buffer
        outbuff[outpos++] = fqflag ? FQH : FAH;
        memcpy(outbuff + outpos, namebuff + namepos, namelen);
        outpos += namelen;
        outbuff[outpos++] = NL;
        codepos += B64;
        //fprintf(stderr, "%s\nlen: %lu\n", namebuff + namepos, seqlen);
        codepos += sqz_seqdecode(codebuff + codepos,
                                 outbuff + outpos,
                                 seqlen,
                                 fqflag,
                                 &outpos);
        outsize -= buffneed;
        namepos += namelen + 1;
        //New sequence
        seqlen   = *(uint64_t *)( codebuff + codepos );
        namelen  = strlen(namebuff + namepos);
        buffneed = ( seqlen * (1 + fqflag) )  + (E + namelen);
    }
    blk->newblk  = 0;
    blk->blkpos  = 0;
    blk->namepos = 0;
    exit:
        blk->prevlen = prevbytes;
        return outpos;
}
