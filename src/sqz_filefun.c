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
                 uint64_t seqlen,
                 uint64_t start,
                 uint64_t bases,
                 uint8_t *buff);


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
    uint8_t  nbases   = len < (32 - startpos) ? len : (32 - startpos);
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


uint64_t sqz_getblkbases(uint8_t  *blkbuff,
                         uint64_t *blkpos)
{
    //fprintf(stderr, "\t\t\tIn decoding\n");
    uint64_t bases    = 0;
    uint64_t blkbytes = 0;
    uint64_t pos      = 0;
    uint64_t pos2     = 0;
    uint8_t  nflag    = 0;
    uint64_t blklen;
    blklen = *(uint64_t *)(blkbuff + pos);
    //fprintf(stderr, "\t\t\tblklen: %lu\n", blklen);
    pos += B64;
    bases += blklen;
    if (blklen)
        blkbytes = getblkbytesize(blklen);
    pos += blkbytes;
    nflag = *(blkbuff + pos);
    pos++;
    //fprintf(stderr, "\t\t\tnflag: %u\n\t\t\tpos: %lu\n", nflag, pos);
    if (nflag)
        bases += countnblk(blkbuff + pos, &pos2);
    pos += pos2;
    *blkpos += pos;
    //blklen = *(uint64_t *)(blkbuff + pos);
    //fprintf(stderr, "\t\t\tNEXT blklen: %lu\n\t\t\tpos %lu\n", blklen, pos);
    return bases;
}

uint64_t sqz_getblkbytes(uint8_t *blkbuff, uint64_t seqlength) {
    uint64_t bases = 0;
    uint64_t blkpos = 0;
    while (bases != seqlength) {
        bases += sqz_getblkbases(blkbuff + blkpos, &blkpos);
    }
    return blkpos;
}


uint64_t chachasmooth(sqz_File *fp, uint8_t *buff, uint64_t start, uint64_t size)
{
    sqz_readblksize(fp->blk, fp->fp);
    uint8_t *blkbuff  = fp->blk->blkbuff;
    uint64_t blksize  = fp->blk->blksize;
    uint64_t blkpos   = fp->blk->blkpos;
    uint64_t namepos  = fp->blk->namepos;
    uint64_t namesize = *(uint64_t *)( blkbuff + ( blksize - B64 ) );
    uint64_t datasize = blksize - B64 - namesize;
    uint64_t seqlen   = 0;
    char *namebuff    = (char *)( blkbuff + datasize );

    uint64_t decodesize;

    //Read sequence length
    seqlen = *(uint64_t *)( blkbuff + blkpos );
    if (start >= seqlen) return 0;
    blkpos += B64;
    decodesize = size;

    //Check if how much sequence fits in buffer
    if ( (seqlen - start) < size) decodesize = (seqlen - start);
    fprintf(stdout, ">%s\n", namebuff + namepos);

    fprintf(stderr, "Decoding %lu bases\n", decodesize);
    uint64_t numbases = 0;
    numbases = pdecode(blkbuff + blkpos, seqlen, start, decodesize, buff);
    if (numbases != decodesize)
        fprintf(stderr, "ERROR\n");

    fprintf(stderr, "blkbytesize: %lu\n",
            sqz_getblkbytes(blkbuff + blkpos, seqlen));
    blkpos += sqz_getblkbytes(blkbuff + blkpos, seqlen);

    namepos += strlen(namebuff + namepos) + 1;
    seqlen = *(uint64_t *)( blkbuff + blkpos);
    fprintf(stderr, "Next seq2:\n\t%s\n", namebuff + namepos);
    fprintf(stderr, "Len:\n\t%lu\n", seqlen);

    return numbases;
}


uint64_t nblkdecode(uint8_t *blkbuff,
                    uint8_t *buff,
                    uint64_t noffset,
                    uint64_t bases,
                    uint64_t *blkpos,
                    uint64_t *blkoffset,
                    uint64_t  blklen)
{
    uint64_t blkpos2 = 0;
    uint64_t nlength = 0;
    uint64_t nbases  = 0;
    //Jump to begining of N block
    uint64_t pos =  B64 + getblkbytesize(blklen);
    //Check nflag
    uint8_t nflag = *(blkbuff + pos);
    if (!nflag) {
        goto exit; //No Nflag so more blk follows
    }
    //Compute N blk length and number of bytes it needed fro encoding
    nlength = countnblk(blkbuff + pos + 1, &blkpos2);
    //Determine how many Ns will fit in buffer
    nbases = bases < ( nlength - noffset ) ? bases : ( nlength - noffset );
    writens(nbases, buff);
    exit:
        //Move to next block if no more Ns are left to decode
        if ( !(nlength - noffset - nbases)  ) {
            *blkpos += pos + blkpos2 + 1;
            *blkoffset = 0;
        }
        return nbases;
}


uint64_t pdecode(uint8_t *blkbuff,
                 uint64_t seqlen,
                 uint64_t start,
                 uint64_t bases,
                 uint8_t *buff)
{
    uint64_t dbases  = 0; //already decoded bases
    uint64_t ebases  = 0; //bases from previous blocks
    uint64_t blkpos  = 0;
    uint64_t blkpos2 = 0;
    uint64_t blklen  = 0;
    uint64_t blkoffset = 0;
    uint64_t noffset   = 0;
    //Go to corresponding block
    while (start >= dbases) {
        blkpos = blkpos2;
        ebases = dbases;
        dbases += sqz_getblkbases(blkbuff + blkpos2, &blkpos2);
    }
    blkoffset = start - ebases; //Position within current block
    dbases  = 0;
    blkpos2 = 0;
    //Main decoding loop. Starting from the beginning of the blk, decode data
    //until number of requested bases is satisfied
    while (dbases < bases) {
        ebases = bases - dbases; //Compute how many bases should be decoded
        //Determine seq blk or n blk
        blklen = *(uint64_t *)(blkbuff + blkpos);
        if (blkoffset >= blklen) {
            //Compute how many Ns have been decoded
            noffset = blkoffset - blklen;
            dbases += nblkdecode(blkbuff + blkpos,
                                 buff + dbases,
                                 noffset,
                                 ebases,
                                 &blkpos,
                                 &blkoffset,
                                 blklen);
            continue;
        }
        //ACGT blk
        uint64_t blkavail = blklen - blkoffset;
        //How much of the block is available
        uint64_t acgt =  ebases < blkavail ? ebases : blkavail;
        //Move blkpos if block is finished
        blkdecode(blkbuff + blkpos + B64,
                  buff + dbases,
                  blkoffset,
                  acgt,
                  blklen);
        dbases += acgt;
        blkoffset += acgt;
    }
    return dbases;
}
