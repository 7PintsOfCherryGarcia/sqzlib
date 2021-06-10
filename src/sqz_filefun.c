#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SQZLIB
#define KLIB
#include "sqz_filefun.h"

#define TWO_BIT_MASK (3)

//Table to change "ACGT" to 0123 else to 4
unsigned char seq_nt4_tableSQZ2[128] = {
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
unsigned char seq_dec_tableSQZ2[128] = {
    'A', 'C','G', 'T',  'N', 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
     4,   4,  4,   4,    4,  4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};




uint64_t pdecode(uint8_t *blkbuff,
                 uint64_t start,
                 uint64_t bases,
                 uint8_t *buff,
                 uint8_t fqflag);


void bit2decode(uint64_t mer,
                uint8_t *decoded,
                uint8_t  startpos,
                uint8_t  nbases,
                uint8_t  len);


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
    uint8_t  nbases   = len < (uint64_t)(32 - startpos) ? (uint8_t)len : (32 - startpos);
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


void writens(uint64_t  numn,
             uint8_t *decoded)
{
    while (numn-- > 0) {
        decoded[numn] = 'N';
    }
}


uint64_t bit2encode(const uint8_t *seq,
                    uint8_t seqlen)
{
    uint64_t result = 0;
    for (uint8_t i = 0; i < seqlen; i++) {
        result = (result << 2) | seq_nt4_tableSQZ2[seq[i]];
    }
    return result;
}


void bit2decode(uint64_t mer,
                uint8_t *decoded,
                uint8_t  startpos,
                uint8_t  nbases,
                uint8_t  len)
{
    uint8_t discard = len - startpos - nbases;
    mer >>= (discard * 2);
    uint8_t byte;
    while(nbases--) {
        byte = mer & TWO_BIT_MASK;                  //Extract lowest 2 bits
        decoded[nbases] = seq_dec_tableSQZ2[byte]; //Write corresponding base
        mer >>= 2;                                 //Discard lowest 2 bits
    }
}


//TODO change to #define
uint64_t getblkbytesize(uint64_t blklen)
{
    return ( (blklen / 32) + ((blklen % 32) > 0) ) * B64;
}


uint64_t countnblk(const uint8_t *buff, uint64_t *pos)
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


uint64_t countqbytes(uint8_t *blkbuff, uint64_t bases)
{
    uint64_t rbases = 0;
    uint64_t nbytes = 0;
    while (rbases < bases) {
        rbases += 1 + (*(blkbuff++) & 31);
        nbytes++;
    }
    return nbytes;
}


uint64_t sqz_gotoblk(uint8_t  *blkbuff,
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


uint64_t sqz_getblkbytes(uint8_t *blkbuff, uint64_t seqlength)
{
    uint64_t bases = 0;
    uint64_t blkpos = 0;
    while (bases != seqlength) {
        bases += sqz_gotoblk(blkbuff + blkpos, &blkpos, 0, 0);
    }
    return blkpos;
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
    //Check nflag
    //uint8_t nflag = *(blkbuff + pos);
    //if (!nflag) {
    //    fprintf(stderr, "\tMust be q blcok if flag is present\n");

    //    goto exit; //No Nflag so more blk follows
    //}
    //fprintf(stderr, "\tNflag from within: %u\n", nflag);
    //Compute N blk length and number of bytes it needed for encoding
    nlength = countnblk(blkbuff + pos, &blkpos2);
    //Determine how many Ns will fit in buffer
    nbases = buffsize < ( nlength - offset ) ? buffsize : ( nlength - offset );
    //fprintf(stderr, "\t%lu Ns\n", nbases);
    writens(nbases, buff);
    //exit:
    //Move to next block if no more Ns are left to decode taking care of
    //quality data if fastq.
    if ( !(nlength - offset - nbases)  ) {
        *blkpos += pos + blkpos2;
        *blkoffset = 0;
    }
    else {
        fprintf(stderr, "NO Ns maybe Qs?");
    }
    return nbases;
}


uint64_t pdecode(uint8_t *blkbuff,
                 uint64_t seqoffset,
                 uint64_t buffsize,
                 uint8_t *buff,
                 uint8_t fqflag)
{
    uint64_t dbytes; //decoded bases
    uint64_t nnum;   //decoded Ns
    uint64_t availbases  = buffsize; //available bases to decode
    uint64_t blkpos   = 0;
    uint64_t blklen   = 0;
    uint64_t noffset  = 0;
    uint8_t  nflag    = 0;
    uint64_t bases    = 0;
    uint64_t tmpbases = 0;
    //Go to corresponding block (block corresponding to current offset)
    tmpbases = sqz_gotoblk(blkbuff, &blkpos, fqflag, seqoffset);
    //Compute number of bases current offset is from block
    seqoffset -= tmpbases;
    //Main decoding loop. Starting from the beginning of the blk, decode data
    //until number of requested bases is satisfied
    dbytes  = 0;
    nnum    = 0;
    while (dbytes < buffsize) {
        //fprintf(stderr, "[pdecode] DECODING\n");
        //Determine seq blk, n blk, or quality block
        blklen = *(uint64_t *)(blkbuff + blkpos);
        //fprintf(stderr, "[pdecode]\tblock size: %lu\t\toffset: %lu\n",
        //        blklen, seqoffset);
        if (seqoffset >= blklen) {
            //fprintf(stderr, "[pdecode]\tNblk\n");
            nflag = *( blkbuff + blkpos + B64 + getblkbytesize(blklen) );
            //fprintf(stderr, "[pdecode]\tnflag: %u\n", nflag);
            //Compute how many Ns have been decoded
            noffset = seqoffset - blklen;
            if (nflag) {
                nnum = nblkdecode(blkbuff + blkpos,
                                  buff + dbytes,
                                  noffset,
                                  availbases,
                                  &blkpos,
                                  &seqoffset,
                                  blklen);
                dbytes += nnum;
                availbases -= nnum;
            }
            else if (fqflag) {
                blkpos += B64 + getblkbytesize(blklen) + 1;
                uint64_t qbytes =  countqbytes(blkbuff + blkpos, bases + dbytes);
                blkpos += qbytes;
                seqoffset = 0;
            }
            else {
                blkpos += B64 + getblkbytesize(blklen) + 1;
                seqoffset = 0;
                blklen = *(uint64_t *)(blkbuff + blkpos);
                //fprintf(stderr, "[pdecode]\tnew blklen: %lu\n", blklen);
            }
            continue;
        }
        //ACGT blk
        //How much of the block is available
        uint64_t blkavail = blklen - seqoffset;
        uint64_t acgt;
        acgt =  blkavail < buffsize - dbytes ? blkavail : buffsize - dbytes;
        //fprintf(stderr, "[pdecode] decoding %lu bases\n", acgt);
        blkdecode(blkbuff + blkpos + B64,
                  buff + dbytes,
                  seqoffset,
                  acgt,
                  blklen);
        dbytes += acgt;
        seqoffset += acgt;
        availbases -= acgt;
    }
    return dbytes;
}


uint64_t sqz_codeblksize(uint8_t *blkbuff, uint8_t fqflag)
{
    uint64_t bases   = 0;
    uint64_t blkpos  = 0;
    uint64_t blkpos2 = 0;
    uint64_t quals   = 0;
    uint64_t blklen;
    uint8_t  nflag;
    uint64_t seqlen = *(uint64_t *)blkbuff;
    //fprintf(stderr, "$$\t\tseqlen: %lu\n", seqlen);
    blkpos += B64;
    while (bases != seqlen) {
        blklen = *(uint64_t *)(blkbuff + blkpos);
        //fprintf(stderr, "$$\t\tblklen: %lu\n", seqlen);
        blkpos += B64;
        if (blklen) {
            bases += blklen;
            blkpos += getblkbytesize(blklen);
        }
        nflag = *(blkbuff + blkpos);
        //fprintf(stderr, "$$\t\tnflag: %u\n", nflag);
        blkpos++;
        if (nflag) {
            bases += countnblk(blkbuff + blkpos, &blkpos2);
            blkpos += blkpos2;
        }
        else if (fqflag) {
            quals = bases - quals;
            blkpos += countqbytes(blkbuff + blkpos, quals);
        }
    }
    return blkpos;
}
