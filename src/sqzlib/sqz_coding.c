#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>

#include "sqz_data.h"

#define TWO_BIT_MASK 3

//Table to change "ACGT" to 0123 else to 4
unsigned char seq_nt4_tableSQZ[128] = {
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
unsigned char seq_dec_tableSQZ[128] = {
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
unsigned char qual_val_table[8] = {33,39,48,55,60,66,70,73};


//TODO Make inline
static  const uint8_t *sqz_findn(const uint8_t *seq)
{
    while (seq_nt4_tableSQZ[*seq] < 4)
        seq++;
    return seq;
}


static uint64_t sqz_bit2encode(const uint8_t *seq, uint8_t seqlen)
{
    uint64_t result = 0;
    for (uint64_t i = 0; i < seqlen; i++)
        result = (result << 2) | seq_nt4_tableSQZ[seq[i]];
    return result;
}


static uint8_t
sqz_bit2decode(const uint64_t *mer, uint8_t *decoded, unsigned char len)
{
    unsigned char byte;
    unsigned char nbase = 0;
    uint64_t code = *mer;
    --len;
    do {
        byte = code & TWO_BIT_MASK;
        decoded[len] = seq_dec_tableSQZ[byte];
        nbase++;
        code >>= 2;
    } while (len-- != 0);
    return nbase;
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


static uint64_t
sqz_seqencode(const uint8_t *seq, uint64_t seqlen, uint8_t *blkbuff)
{
    uint64_t       blkpos = 0;             //position in blkbuff
    const uint8_t *lstop  = seq;           //Track position in sequence
    const uint8_t *nptr   = seq + seqlen;  //End of string
	  const uint8_t *npos   = NULL;          //Track positions where N occurs
    uint64_t nn     = 0;                   //Number of Ns
	  uint64_t blen   = 0;                   //Length of segment before first N
    uint8_t  flag1 = 255;
    uint8_t  flag2 = 0;
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
            //Indicate that a non-ACGT block follows
            memcpy(blkbuff + blkpos, &flag1, 1);
            blkpos++;
            //Encode non-ACGT block
            blkpos += sqz_nblkcode(blkbuff + blkpos, nn);
            //Trace next base after sequence of Ns
            lstop = npos;
        }
    } while (*npos);
    //Detect and encode trailing bases
    blen = nptr - lstop;
    if (blen) {
        memcpy(blkbuff + blkpos, &blen, B64);
        blkpos += B64;
        blkpos += sqz_blkcode((uint64_t *)(blkbuff + blkpos), lstop, blen);
        /*
          When a sequence has only been partialy loaded when encoding loop
          finishes, a flag with a value of zero is written to the buffer.
          This is so that when decoding, a quality block can be distinguished
          from more sequence blocks. This flag is needed only when the loaded
          sequence is truncated at ACGT bases. If the sequence is truncated at
          non ACGT bases there will be no trailing bases to decode (blen == 0)
          and the corresponding flag would have already been written.
        */
        memcpy(blkbuff + blkpos, &flag2, 1);
        blkpos++;
    }
    return blkpos;
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


static uint64_t sqz_qualencode(const uint8_t *qual, uint8_t *blkbuff)
{
    uint64_t blkpos = 0;
    uint8_t q       = sqz_8binqual(*qual);
    uint8_t c       = 0;
    uint8_t code    = 0;
    while(*qual) {
        if ( sqz_8binqual(*(qual + 1)) == q) {
            c++;
            if (c == 31) {
                //Encode
                qual++;
                code = code | (q & 224);
                code = code | c;
                memcpy(blkbuff + blkpos, &code, 1);
                blkpos += 1;
                c = 0;
                code = 0;
                q = sqz_8binqual(*(qual + 1));
            }
            qual++;
            continue;
        }
        //Encode
        code = code | (q & 224);
        code = code | c;
        memcpy(blkbuff + blkpos, &code, 1);
        blkpos += 1;
        qual++;
        q = sqz_8binqual(*qual);
        c = 0;
        code = 0;
    }
    return blkpos;
}


//TODO change to #define
static uint64_t getblkbytesize(uint64_t blklen)
{
    return ( (blklen / 32) + ((blklen % 32) > 0) ) * B64;
}


static uint8_t sqz_fastqheadblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    uint8_t  *blkbuff    = blk->blkbuff;
    uint64_t  blkpos     = blk->blkpos;
    uint8_t  *seqbuffer  = sqz->seqbuffer;
    uint8_t  *qltbuffer  = sqz->qualbuffer;
    uint64_t  sqzsize    = sqz->offset;
    uint8_t  *seq     = NULL;
    uint8_t  *qlt     = NULL;
    uint64_t  seqlen  = 0;
    uint64_t  seqread = 0;
    uint64_t  k       = 0;
    while ( k < sqzsize ) {
        seqlen = *(uint64_t *)( seqbuffer + k );
        k += B64;
        memcpy(blkbuff + blkpos, &seqlen, B64);
        blkpos += B64;
        seq  = seqbuffer + k;
        qlt  = qltbuffer + k;
        seqread = seqlen < (sqzsize - k) ? seqlen : sqzsize - k - 1;
        k += seqread + 1;
        blkpos += sqz_seqencode(seq, seqread, blkbuff + blkpos);
        blkpos += sqz_qualencode(qlt, blkbuff + blkpos);
    }
    //More sequence to encode?
    if (sqz->endflag) {
        //Unset new block flag.
        blk->newblk = 0;
        goto exit;
    }
    memcpy(blkbuff + blkpos, sqz->namebuffer, sqz->namepos);
    blkpos += sqz->namepos;
    memcpy(blkbuff + blkpos, &(sqz->namepos), B64);
    blkpos += B64;
    sqz->cmpflag = 1;
    exit:
        blk->blkpos = blkpos;
        return 1;
}


static uint8_t sqz_fastqtailblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    uint64_t seqlen  = sqz->prevlen;
    uint64_t seqleft = seqlen - sqz->seqread;
    uint8_t *seqbuff = sqz->seqbuffer;
    uint8_t *qltbuff = sqz->qualbuffer;
    uint8_t *blkbuff = blk->blkbuff;
    uint64_t blksize = blk->blksize;
    uint64_t blkpos  = blk->blkpos;
    if (sqz->endflag) seqleft = LOAD_SIZE;
    if ( (blkpos + getblkbytesize(seqleft)) > blksize ) {
        blksize *= 2;
        blkbuff  = realloc(blk->blkbuff, blksize);
        blk->blkbuff  = blkbuff;
        blk->blksize  = blksize;
        blk->cmpsize *= 2;
        blk->cmpbuff  = realloc(blk->cmpbuff, blk->cmpsize);
    }
    blkpos += sqz_seqencode(seqbuff, seqleft, blkbuff + blkpos);
    blkpos += sqz_qualencode(qltbuff, blkbuff + blkpos);
    if (sqz->endflag) {
        blk->newblk = 0;
        goto exit;
    }
    blk->newblk = 1;
    //Copy name data making sure it fits
    if ( (blkpos + sqz->namepos) >= blksize ) {
        blksize += sqz->namepos + B64;
        blkbuff = realloc(blk->blkbuff, blksize);
        blk->blkbuff = blkbuff;
        blk->blksize = blksize;
    }
    memcpy(blkbuff + blkpos, sqz->namebuffer, sqz->namepos);
    blkpos += sqz->namepos;
    memcpy(blkbuff + blkpos, &(sqz->namepos), B64);
    blkpos += B64;
    sqz->cmpflag = 1;
    exit:
        blk->blkpos = blkpos;
        return 1;
}


static uint8_t sqz_fastaheadblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    uint8_t *blkbuff   = blk->blkbuff;
    uint64_t blkpos    = blk->blkpos;
    uint8_t *seqbuffer = sqz->seqbuffer;
    uint64_t sqzsize   = sqz->offset;
    uint8_t *seq = NULL;
    uint64_t seqlen = 0;
    uint64_t seqread = 0;
    uint64_t k = 0;
    while ( k < sqzsize ) {
        seqlen = *(uint64_t *)( seqbuffer + k );
        k += B64;
        memcpy(blkbuff + blkpos, &seqlen, B64);
        blkpos += B64;
        seq = seqbuffer + k;
        seqread = seqlen < (sqzsize - k) ? seqlen : sqzsize - k - 1;
        k += seqread + 1;
        blkpos += sqz_seqencode(seq, seqread, blkbuff + blkpos);
    }
    //Indicate if there is more sequence to read
    if (sqz->endflag) {
        //Unset new block flag.
        blk->newblk = 0;
        goto exit;
    }
    memcpy(blkbuff + blkpos, sqz->namebuffer, sqz->namepos);
    blkpos += sqz->namepos;
    memcpy(blkbuff + blkpos, &(sqz->namepos), B64);
    blkpos += B64;
    sqz->cmpflag = 1;
    exit:
        blk->blkpos = blkpos;
        return 1;
}


static uint8_t sqz_fastatailblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    uint64_t seqlen  = sqz->prevlen;
    uint64_t seqleft = seqlen - sqz->seqread;
    uint8_t *seqbuff = sqz->seqbuffer;
    uint8_t *blkbuff = blk->blkbuff;
    uint64_t blksize = blk->blksize;
    uint64_t blkpos  = blk->blkpos;

    if (sqz->endflag) seqleft = LOAD_SIZE;
    if ( (blkpos + getblkbytesize(seqleft)) > blksize ) {
        blksize *= 2;
        blkbuff = realloc(blk->blkbuff, blksize);
        blk->blkbuff = blkbuff;
        blk->blksize = blksize;
        blk->cmpsize *= 2;
        blk->cmpbuff = realloc(blk->cmpbuff, blk->cmpsize);
    }
    //encode sequence
    blkpos += sqz_seqencode(seqbuff, seqleft, blkbuff + blkpos);
    //Indicate if there is more sequence to read
    if (sqz->endflag) {
        blk->newblk = 0;
        goto exit;
    }
    blk->newblk = 1;
    //Copy name data making sure it fits
    if ( (blkpos + sqz->namepos) >= blksize ) {
        blksize += sqz->namepos + B64;
        blkbuff = realloc(blk->blkbuff, blksize);
        blk->blkbuff = blkbuff;
        blk->blksize = blksize;
    }
    memcpy(blkbuff + blkpos, sqz->namebuffer, sqz->namepos);
    blkpos += sqz->namepos;
    memcpy(blkbuff + blkpos, &(sqz->namepos), B64);
    blkpos += B64;
    sqz->cmpflag = 1;
    exit:
        blk->blkpos = blkpos;
        return 1;
}


char sqz_fastXencode(sqzfastx_t *sqz, sqzblock_t *blk, uint8_t fqflag)
{
    if (blk->newblk) {
        if (fqflag) return sqz_fastqheadblk(sqz, blk);
        return sqz_fastaheadblk(sqz, blk);
    }
    if (fqflag) return sqz_fastqtailblk(sqz, blk);
    return sqz_fastatailblk(sqz, blk);
}


/*
################################################################################
  Decoding functions
################################################################################
*/
static void sqz_qdecode(uint8_t  *blkbuff,
                        uint8_t  *buff,
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
        pos += 1 + (*(blkbuff++) & 31);
        nbyte++;
    }
    prevqn = pos - offset;
    rquals -= prevqn;
    q = (*(blkbuff - 1) & 224) >> 5;
    for (i = 0; i < prevqn; i++)
        *buff++ = qual_val_table[q];
    while (rquals) {
        i = 1 + (*blkbuff & 31);       //Get count
        q = (*(blkbuff++) & 224) >> 5; //Get symbol
        //TODO change to ternary operator ?:
        if (i > rquals)
            i = rquals;
        rquals -= i;
        //Write symbols
        for (j = 0; j < i; j++)
            *buff++ = qual_val_table[q];
    }
}

/*
  Counts number of Ns in an N block
  Returns number of Ns
  Modifies pos with number of bytes decoded
*/
static uint64_t countnblk(const uint8_t *buff, uint64_t *pos)
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


static uint64_t countqbytes(uint8_t *codebuff, uint64_t bases)
{
    uint64_t rbases = 0;
    uint64_t nbytes = 0;
    while (rbases < bases) {
        rbases += 1 + (*(codebuff++) & 31);
        nbytes++;
    }
    return nbytes;
}


static uint64_t sqz_gotoqblk(uint8_t *blkbuff,
                             uint64_t qoffset,
                             uint64_t *rblklen,
                             uint64_t *blkbases)
{
    uint64_t blklen;
    uint64_t blkpos = 0;
    uint8_t  nflag;
    uint64_t bases = 0;
    uint64_t dbases = *blkbases;
    uint64_t pos;
    while (1) {
        blklen = *(uint64_t *)(blkbuff + blkpos);
        blkpos += B64;
        if (blklen)
            blkpos += getblkbytesize(blklen);
        bases += blklen;
        nflag = *(blkbuff + blkpos);
        blkpos++;
        if (nflag) {
            bases += countnblk(blkbuff + blkpos, &pos);
            blkpos += pos;
        }
        else {
            if (bases + dbases >= qoffset)
                break;
            blkpos += countqbytes(blkbuff + blkpos, bases);
            dbases += bases;
            bases = 0;
        }
    }
    *rblklen  = bases;
    *blkbases = dbases;
    return blkpos;
}


static void sqz_blkqdecode(uint8_t *blkbuff,
                           uint8_t *buff,
                           uint64_t nquals,
                           uint64_t qoffset)
{
    uint64_t ndecoded  = 0;
    uint64_t blkpos    = 0;
    uint64_t prevq     = 0;
    uint64_t qinblk    = 0;
    uint64_t blkoffset = 0;
    uint64_t rquals    = 0;
    while (nquals) {
        //Given offset, go to corresponding quality block. Storing the number
        //of qualities available to extract in current block and the number
        //of qualities from all previous blocks.
        blkpos = sqz_gotoqblk(blkbuff, qoffset, &qinblk, &prevq);
        //Compute how many bases can be decoded from current block
        blkoffset = qoffset - prevq;
        rquals = nquals < qinblk - blkoffset ? nquals : qinblk - blkoffset;
        //Load qualities from block
        sqz_qdecode(blkbuff + blkpos, buff, rquals, blkoffset);
        //keep track of how many qualities are loaded from current blk
        buff += rquals;
        ndecoded += rquals;
        qoffset += rquals;
        nquals  -= rquals;
        blkpos += countqbytes(blkbuff + blkpos, qinblk);
        //TODO possible bug if the entire block has not been decoded
        //Why am I doing this?
        prevq += qinblk;
    }
}


static uint64_t sqz_gotoblk(uint8_t  *blkbuff,
                            uint64_t *blkpos,
                            uint8_t   fqflag,
                            uint64_t  offset)
{
    uint64_t bases  = 0;
    uint64_t pbases = 0;
    uint64_t quals  = 0;
    uint64_t pos    = 0;
    uint64_t pos2   = 0;
    uint64_t ppos   = 0;
    uint8_t  nflag  = 0;
    uint64_t blklen;
    while (offset >= bases ) {
        //Store how many bases have been counted
        pbases = bases;
        ppos = pos;
        blklen = *(uint64_t *)(blkbuff + pos);
        pos += B64;
        bases += blklen;
        if (blklen)
            pos += getblkbytesize(blklen);
        nflag = *(blkbuff + pos);
        pos++;
        if (nflag) {
            bases += countnblk(blkbuff + pos, &pos2);
            pos += pos2;
        }
        //If flag is set, there are quality bytes to decode
        else if (fqflag) {
            //Should only pass number of qualities present in block
            quals = bases - quals;
            pos += countqbytes(blkbuff + pos, quals);
        }
    }
    *blkpos += ppos;
    return pbases;
}


/*
  Writes numn Ns into *decoded
*/
static void sqz_writens(uint64_t  numn, uint8_t *decoded)
{
    while (numn-- > 0)
        decoded[numn] = 'N';
}


uint64_t nblkdecode(uint8_t *blkbuff,
                    uint8_t *buff,
                    uint64_t offset,
                    uint64_t buffsize,
                    uint64_t *blkpos,
                    uint64_t *blkoffset,
                    uint64_t  blklen)
{
    uint64_t blkpos2 = 0;
    uint64_t nlength = 0;
    uint64_t nbases  = 0;
    //Jump to begining of N block
    uint64_t pos =  B64 + getblkbytesize(blklen) + 1;
    //Compute N blk length and number of bytes it needed for encoding
    nlength = countnblk(blkbuff + pos, &blkpos2);
    //Determine how many Ns will fit in buffer
    nbases = buffsize < ( nlength - offset ) ? buffsize : ( nlength - offset );
    sqz_writens(nbases, buff);
    //Move to next block if no more Ns are left to decode taking care of
    //quality data if fastq.
    if ( !(nlength - offset - nbases)  ) {
        *blkpos += pos + blkpos2;
        *blkoffset = 0;
    }
    return nbases;
}


static void bit2decode(uint64_t mer,
                       uint8_t *outbuff,
                       uint8_t  startpos,
                       uint8_t  nbases,
                       uint8_t  len)
{
    uint8_t discard = len - startpos - nbases;
    mer >>= (discard * 2);
    uint8_t byte;
    while(nbases--) {
        byte = mer & TWO_BIT_MASK;                  //Extract lowest 2 bits
        outbuff[nbases] = seq_dec_tableSQZ[byte];   //Write corresponding base
        mer >>= 2;                                  //Discard lowest 2 bits
    }
}


void blkdecode(const uint8_t *codebuff,
               uint8_t       *decodebuff,
               uint64_t       blkoffset,
               uint64_t       len,
               uint64_t       blklen)
{
    uint64_t codepos   = 0;
    uint64_t decodepos = 0;
    uint64_t blknum    = 0;
    //Compute start position within mer
    uint8_t  startpos = blkoffset % 32;
    //Number of bases to decode from these first 64 bits
    //TODO check these casts
    uint8_t  nbases   = len < (uint64_t)(32 - startpos) ?
                        (uint8_t)len : (32 - startpos);
    if (nbases) blknum++;
    //Compute how many blks encode the number of bases requested
    blknum += ( (len - nbases) / 32) + ( ( (len - nbases)  % 32) > 0 ) - 1;
    //Place mer on corresponding byte
    uint64_t *mer = (uint64_t *)codebuff + blkoffset / 32;
    //fprintf(stderr, "\tstarting mer: %lu\n", *mer);
    uint64_t blkidx;
    for (blkidx = 0; blkidx < blknum; blkidx++) {
        bit2decode(*(mer + blkidx),
                   decodebuff + decodepos,
                   startpos,
                   nbases,
                   32);
        //Keep traked of amount of decoded data: 64bit per mer
        codepos += B64;
        //Move sequence pointer by amount of bases decoded
        decodepos += nbases;
        //Next iterations whould always start at base 0
        startpos = 0;
        //Next iterations will always be 32 bases
        nbases = 32;
    }
    len -=decodepos;
    //How large is this last mer?
    uint64_t lefy = blklen - blkoffset - decodepos + startpos;
    uint8_t lastmer = lefy < 32 ? lefy : 32;
    bit2decode(*(mer + blkidx),
               decodebuff + decodepos,
               startpos,
               len,
               lastmer);
    codepos += B64;
    decodepos += len;
}


static uint64_t pdecode(uint8_t *codebuff,
                        uint64_t seqoffset,
                        uint64_t buffsize,
                        uint8_t *outbuff,
                        uint8_t fqflag)
{
    uint64_t outpos = 0; //decoded bases
    uint64_t nnum = 0;   //decoded Ns
    uint64_t availbases  = buffsize; //available bases to decode
    uint64_t blkpos   = 0;
    uint64_t blklen   = 0;
    uint64_t noffset  = 0;
    uint8_t  nflag    = 0;
    uint64_t bases    = 0;
    uint64_t tmpbases = 0;
    //Go to corresponding block (block corresponding to current offset)
    tmpbases = sqz_gotoblk(codebuff, &blkpos, fqflag, seqoffset);
    //Compute number of bases current offset is from block
    seqoffset -= tmpbases;
    //Main decoding loop. Starting from the beginning of the blk, decode data
    //until number of requested bases is satisfied
    while (outpos < buffsize) {
        //Determine seq blk, n blk, or quality block
        blklen = *(uint64_t *)(codebuff + blkpos);
        if (seqoffset >= blklen) {
            nflag = *( codebuff + blkpos + B64 + getblkbytesize(blklen) );
            if (nflag) {
                //Compute how many Ns have been decoded
                noffset = seqoffset - blklen;
                bases += noffset;
                nnum = nblkdecode(codebuff + blkpos,
                                  outbuff + outpos,
                                  noffset,
                                  availbases,
                                  &blkpos,
                                  &seqoffset,
                                  blklen);
                bases += nnum;
                outpos += nnum;
                availbases -= nnum;
            }
            else if (fqflag) {
                blkpos += B64 + getblkbytesize(blklen) + 1;
                uint64_t qbytes = countqbytes(codebuff + blkpos, bases);
                blkpos += qbytes;
                seqoffset = 0;
            }
            else {
                blkpos += B64 + getblkbytesize(blklen) + 1;
                seqoffset = 0;
                blklen = *(uint64_t *)(codebuff + blkpos);
            }
            continue;
        }
        bases += blklen;
        //ACGT blk
        //How much of the block is available
        uint64_t blkavail = blklen - seqoffset;
        uint64_t acgt;
        acgt =  blkavail < buffsize - outpos ? blkavail : buffsize - outpos;
        blkdecode(codebuff + blkpos + B64,
                  outbuff + outpos,
                  seqoffset,
                  acgt,
                  blklen);
        outpos += acgt;
        seqoffset += acgt;
        availbases -= acgt;
    }
    return outpos;
}


/*
  Decodes outsize bytes from blkbuff
  Returns number of bytes written
*/
static uint64_t sqz_pdecode(uint8_t *codebuff,
                            uint8_t *outbuff,
                            uint64_t outsize,
                            uint64_t seqlen,
                            char *name,
                            uint8_t  fqflag)
{
    uint64_t outpos = 0;
    outbuff[outpos++] = fqflag ? FQH : FAH;        //> or @
    memcpy(outbuff + outpos, name, strlen(name));  //seq name
    outpos += strlen(name);
    outbuff[outpos++] = NL;                        //newline
    outsize -= outpos;                             //( 2 + strlen(name) );
    //Determine number of bases that fit in buffer
    uint64_t todecode = (seqlen <= outsize) ? seqlen : outsize;
    pdecode(codebuff,
            0,
            todecode,
            outbuff + outpos,
            fqflag);
    outpos += todecode;
    outsize -= todecode;
    //Decode quality if enough space
    if ( (outsize > 3) && fqflag ) {
        outbuff[outpos++] = NL;
        outbuff[outpos++] = '+';
        outbuff[outpos++] = NL;
        outsize -= 3;
        sqz_blkqdecode(codebuff,         //data buffer
                       outbuff + outpos, //output buffer
                       outsize,          //requested qualities/size of outbuff
                       0);               //already decoded qualities
        outpos += outsize;
    }
    return outpos;
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
    //Compute how much sequence can be decoded
    seq2decode = seqlen < decoded ? 0 : seqlen - decoded;
    //Check how much of leftover sequence fits in output buffer
    seq2decode = seq2decode < outsize ? seq2decode : outsize;
    //Compute how many qualities have been decoded
    qdecoded = seqlen < decoded ? decoded - seqlen - 3 : 0;
    //Check how much of leftover qualities fit in output buffer
    //account for decoded sequences filling outbuff
    //qdecoded = qdecoded < (outsize - seq2decode) ?
    //           qdecoded : (outsize - seq2decode);
    //Compute how much buffer will be available
    bleft = outsize - seq2decode;
    if (seq2decode) {
        pdecode(codebuff,
                decoded,
                seq2decode,
                outbuff,
                fqflag);
        if (bleft < 3) goto exit; //Exit function if new lines don't fit
    }
    //In case of fastq, we have to deal with qualities
    if (fqflag) {
        //Add corresponfing new lines if begining if quality string
        if (!qdecoded) {
            outbuff[seq2decode++] = NL;
            outbuff[seq2decode++] = '+';
            outbuff[seq2decode++] = NL;
            bleft -= 3;
            if (!bleft) goto exit; //Exit function if no more space in outbuffer
        }
        //How many qualities can be decoded
        qlt2decode = seqlen - qdecoded < bleft ? (seqlen - qdecoded) : bleft;
        sqz_blkqdecode(codebuff,                     //data buffer
                       outbuff + seq2decode,         //output buffer
                       qlt2decode,                   //requested qualities
                       qdecoded);                    //already decoded qualities

    }
    exit:
        return seq2decode + qlt2decode;
}


/*
  Returns size of a sequence code block
*/
static uint64_t sqz_codeblksize(uint8_t *codebuff, uint8_t fqflag)
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
            //Advance to after the block
            codepos += getblkbytesize(blklen);
        }
        //Read N flag
        nflag = *(codebuff + codepos);
        codepos++;
        //Check and account for Ns
        if (nflag) {
            bases += countnblk(codebuff + codepos, &blkpos2);
            codepos += blkpos2;
        }
        //Account for quality data
        else if (fqflag)
            codepos += countqbytes(codebuff + codepos, bases);
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
    //uint64_t blknum = blklen / 32;
    //1 less than the total number of iterations needed
    uint64_t blknum  = ( (blklen / 32) + ( (blklen % 32) > 0) ) - 1;
    //fprintf(stderrm "\t\t$$blknum: %lu\n");
    uint64_t codepos = 0;
    uint64_t outpos  = 0;
    uint64_t blkidx  = 0;
    uint64_t blkleft  = (blklen % 32) == 0 ? 32 : blklen % 32;
    uint64_t *mer;
    mer = (uint64_t *)codebuff;
    for (blkidx = 0; blkidx < blknum; blkidx++) {
        sqz_bit2decode(mer + blkidx, outbuff + outpos, 32);
        //Keep traked of amount of decoded data: 64bit per mer
        codepos += B64;
        //Move sequence pointer by amount of bases decoded
        outpos += 32;
    }
    //if (blkleft) {
    sqz_bit2decode(mer + blkidx, outbuff + outpos, blkleft);
    codepos += B64;
    outpos += blkleft;
    //}
    *wbytes += outpos;
    return codepos;
}


/*
  Decodes length quality values from codebuff
  Returns: number of bytes decoded
 */
static uint64_t sqz_qualdecode(const uint8_t *codebuff,
                               uint8_t *qualstr,
                               uint64_t length)
{
    uint64_t byte    = 0;
    uint64_t offset  = 0;
    uint64_t decoded = 0;
    uint8_t  code    = 0;
    uint8_t  count   = 0;
    uint8_t  q       = 0;
    while (decoded != length) {
        code = *(codebuff + byte); //get byte value
        count = code & 31;         //get current symbol count
        q = code >> 5;             //get symbol index value
        for (int i = 0; i <= count; i++)
            qualstr[offset++] = qual_val_table[q];
        decoded += (count + 1); //Count[0-31], real length [1-32]
        byte++;
    }
    //fprintf(stderr, "\t\treturning %lu bytes decoded\n", byte);
    return byte;
}

/*
Loads length bases into decodebuffer
Returns number of bytes decoded
Updates wbytes with number of bytes written
*/
uint64_t sqz_seqdecode(const uint8_t *codebuff,
                       uint8_t       *outbuff,
                       uint64_t       length,
                       char           qflag,
                       uint64_t      *wbytes)
{
    uint64_t seqlen     = length;
    uint64_t codepos    = 0;
    uint64_t blklen     = 0;
    uint64_t qltdecoded = 0;
    uint64_t qltnum     = 0;
    uint64_t outpos     = 0;
    uint8_t nflag;
    uint64_t nnum = 0;
    uint64_t nbytes;
    //fprintf(stderr, "**Decoding sequence of length: %lu\n", length);
    while (length > 0) {
        blklen = *(uint64_t *)(codebuff + codepos);
        qltnum += blklen;
        //fprintf(stderr, "\t**blklen: %lu\n", blklen);
        //fprintf(stderr, "\t\t**prev: %lu\n",
        //        *(uint64_t *)(codebuff + codepos - 1));
        //fprintf(stderr, "\t\t**next: %lu\n",
        //        *(uint64_t *)(codebuff + codepos + 1));
        codepos += B64;
        codepos += sqz_blkdecode(codebuff + codepos,
                                 outbuff + outpos,
                                 &outpos,
                                 blklen);
        length -= blklen;
        //fprintf(stderr, "\t**length left: %lu\n", length);
        if (length) {     //Partial sequence decoding
            nflag = *(codebuff + codepos); //Read N flag
            codepos++;
            if (nflag) {
                //Deal with n block
                nnum = countnblk(codebuff + codepos, &nbytes);
                qltnum += nnum;
                //fprintf(stderr, "\t**nblock: %lu\n", nnum);
                codepos += nbytes;
                sqz_writens(nnum, outbuff + outpos);
                outpos += nnum;
                length -= nnum;
            }
            else {
                //fprintf(stderr, "\t**not nflag\n");
                //Deal with something else which can be
                if (qflag) {
                    //fprintf(stderr, "\t**qualflag\n");
                    //Deal with quality string
                    /*
                    Even though entire sequence has not been decoded completely,
                    there is quality data to be decoded. This happens because in
                    this particular sequence, during the encoding process a
                    buffer was filled before entire sequence could be loaded.
                    In this case, what was loaded of the sequence is encoded,
                    and then the rest of the sequence finishes loading.
                    Resulting in the following patter in the code block:
                        blklen:seqcode:qualcode:blklen:seqcode:qualcode
                    Instead of:
                        blklen:seqcode:qualcode
                    */
                    //TODO Completely change this hacky way of decoding
                    //Decode nucleotides first
                    //Put qualities at end of nucleotides after extra characters
                    //codepos += sqz_qualdecode(codebuff + codepos,
                    //                          outbuff + seqlen + 3, //!!!!!TODO
                    //                          blklen);
                    //fprintf(stderr,
                    //      "\t**decoding %lu qualities nnum: %lu qltnum: %lu\n",
                    //        blklen, nnum, qltnum);
                    codepos += sqz_qualdecode(codebuff + codepos,
                                              outbuff + seqlen + 3 + qltdecoded,
                                              qltnum);
                    //Change sequence length to reflect part of the quality
                    //string has been decoded
                    //seqlen -= blklen;
                    //qltdecoded += blklen;
                    qltdecoded += qltnum;
                    qltnum = 0;
                }
            }
            continue;
        }
        //nflag = *(codebuff + codepos);
        //Skip the nflag
        codepos++;
    }
    //Decode quality strings for fastq data
    if (qflag) {
        outbuff[outpos++] = NL;
        outbuff[outpos++] = '+';
        outbuff[outpos++] = NL;
        codepos += sqz_qualdecode(codebuff + codepos,
                                  outbuff + outpos + qltdecoded,
                                  seqlen - qltdecoded);
        //outpos += (seqlen - qltdecoded);
        outpos += seqlen;
    }
    outbuff[outpos++] = NL;
    *wbytes += outpos;
    return codepos;
}


uint64_t sqz_fastXdecode(sqzblock_t *blk,   //Data block
                         uint8_t *outbuff,  //Array to place decoded data
                         uint64_t buffsize, //Size of outbuff
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
    //New lines per sequence
    uint8_t E = fqflag ? 6 : 3;
    /*
      To fit current sequence, length of sequence bytes is needed (2x if fastq
      for quality values). Plus 3 header/newline bytes (6 for fastq).
    */
    uint64_t buffneed = ( seqlen * (1 + fqflag) )  + (E + namelen);
    //Check if there is sequence that was not completely decoded
    if (prevbytes) {
        prevbytes -= (namelen + 2);
        outpos = sqz_rdecode(codebuff + codepos + B64,
                             outbuff,
                             buffsize,
                             seqlen,
                             prevbytes,
                             fqflag);
        prevbytes += outpos;
        //Check that no more sequence needs decoding. (Move to next sequence)
        if (prevbytes < buffneed - E - namelen)
            goto exit;
        outbuff[outpos++] = NL;
        buffsize -= outpos;
        prevbytes = 0;
        //Update code buffer position to next sequence
        codepos += sqz_codeblksize(codebuff + codepos, fqflag);
        //Get new sequence info
        //TODO: Overflow read when end of block has been reached
        namepos += namelen + 1;
        namelen  = strlen(namebuff + namepos);
        seqlen   = *(uint64_t *)( codebuff + codepos );
        buffneed = ( seqlen * (1 + fqflag) )  + (E + namelen);
    }
    //While there is data to decode
    //fprintf(stderr, "codepos: %lu datasize: %lu\n", codepos, datasize);
    while ( codepos < datasize ) {
        //Test if there is enough space for current sequence
        if ( buffsize < buffneed ) {
            //Test if at least: ">"|"@" + name length + '\n' + 1 base fits
            if (buffsize < namelen + 3) {
                //Can't decode this sequence, just return buffer
                prevbytes    = 0;
                blk->blkpos  = codepos;
                blk->namepos = namepos;
                goto exit;
            }
            //Decode as much as possible
            prevbytes = sqz_pdecode(codebuff + codepos + B64,
                                    outbuff + outpos,
                                    buffsize,
                                    seqlen,
                                    namebuff + namepos,
                                    fqflag);
            outpos += prevbytes;
            blk->blkpos  = codepos;
            blk->namepos = namepos;
            goto exit;
        }
        //Sequence fits and there is still space in buffer
        //Add header symbol: ">" or "@"
        outbuff[outpos++] = fqflag ? FQH : FAH;
        //Copy sequence name
        memcpy(outbuff + outpos, namebuff + namepos, namelen);
        outpos += namelen;
        outbuff[outpos++] = NL;
        codepos += B64;
        codepos += sqz_seqdecode(codebuff + codepos,
                                 outbuff + outpos,
                                 seqlen,
                                 fqflag,
                                 &outpos);
        buffsize -= buffneed;
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
