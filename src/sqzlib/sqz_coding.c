#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>

//#define SQZLIB
//#define KLIB
#include "sqz_coding.h"

uint64_t sqz_codeblksize(uint8_t *blkbuff, uint8_t fqflag);

void writens(uint64_t  numn,
             uint8_t *decoded);

uint64_t countnblk(const uint8_t *buff, uint64_t *pos);

uint64_t getblkbytesize(uint64_t blklen);

uint64_t pdecode(uint8_t *blkbuff,
                 uint64_t start,
                 uint64_t bases,
                 uint8_t *buff,
                 uint8_t fqflag);

uint64_t sqz_getblkbytes(uint8_t *blkbuff, uint64_t seqlength);

uint64_t countqbytes(uint8_t *blkbuff, uint64_t bases);


char sqz_fastXencode(sqzfastx_t *sqz, sqzblock_t *blk, uint8_t fqflag)
{
    if (blk->newblk) {
        if (fqflag) return sqz_fastqheadblk(sqz, blk);
        return sqz_fastaheadblk(sqz, blk);
    }
    if (fqflag) return sqz_fastqtailblk(sqz, blk);
    return sqz_fastatailblk(sqz, blk);
}


char sqz_fastqheadblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    //fprintf(stderr, "Encoding from head\n");
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


char sqz_fastqtailblk(sqzfastx_t *sqz, sqzblock_t *blk)
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


char sqz_fastaheadblk(sqzfastx_t *sqz, sqzblock_t *blk)
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


char sqz_fastatailblk(sqzfastx_t *sqz, sqzblock_t *blk)
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


uint64_t sqz_seqencode(const uint8_t *seq, uint64_t lentocode, uint8_t *blkbuff)
{
    uint64_t       blkpos = 0;
    const uint8_t *lstop = seq;            //Track position in sequence
    const uint8_t *nptr = seq + lentocode; //End of string
	  const uint8_t *npos;                   //Track positions where N occurs
    uint64_t       nn;                     //Number of Ns
	  uint64_t       blen = 0;               //Length of segment before first N
    uint64_t       nbases = 0;             //Number of bases encoded
    uint8_t        flag1 = 255;
    uint8_t        flag2 = 0;
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
            nbases += blen;
            //Indicate that a non ACGT block follows
            memcpy(blkbuff + blkpos, &flag1, 1);
            blkpos++;
            //Encode non ACGT block
            blkpos += sqz_nblkcode(blkbuff + blkpos, nn);
            nbases += nn;
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
        nbases += blen;
        /*
          Check that no more sequence needs to be encoded. If more sequence is
          needed write corresponding byte.
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


uint64_t sqz_qualencode(const uint8_t *qual, uint8_t *blkbuff)
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


static uint64_t sqz_pdecode(uint8_t *blkbuff,
                            uint8_t *buff,
                            uint64_t buffsize,
                            uint64_t seqlen,
                            char *name,
                            uint8_t  fqflag)
{
    uint64_t wbytes = 0;
    buff[wbytes++] = fqflag ? FQH : FAH;        //> or @
    memcpy(buff + wbytes, name, strlen(name));  //seq name
    wbytes += strlen(name);
    buff[wbytes++] = NL;                        //newline
    buffsize -= 2 + strlen(name);
    //Decode number of bases that fit in buffer
    uint64_t todecode = seqlen <= buffsize ? seqlen : buffsize;
    pdecode(blkbuff,
            0,
            todecode,
            buff + wbytes,
            fqflag);
    wbytes += todecode;
    buffsize -= todecode;
    //Decode quality if enough space
    if ( (buffsize > 3) && fqflag ) {
        buff[wbytes++] = NL;
        buff[wbytes++] = '+';
        buff[wbytes++] = NL;
        buffsize -= 3;
        sqz_blkqdecode(blkbuff,         //data buffer
                       buff + wbytes,   //output buffer
                       buffsize,        //requested qualities
                       0);               //already decoded qualities
        wbytes += buffsize;
    }
    return wbytes;
}


static uint64_t sqz_rdecode(uint8_t  *blkbuff,
                            uint8_t  *buff,
                            uint64_t buffsize,
                            uint64_t seqlen,
                            uint64_t prevbytes,
                            uint8_t  fqflag)
{
    uint64_t bleft     = 0;
    uint64_t todecode  = 0;
    uint64_t qtodecode = 0;
    uint64_t qdecoded  = 0;
    uint8_t  spacer = 0;
    //Compute how much sequence can be decoded
    todecode = seqlen < prevbytes ? 0 : seqlen - prevbytes;
    qdecoded = seqlen < prevbytes ? prevbytes - seqlen - 3 : 0;
    //fprintf(stderr, "%lu bases need decoding\n", todecode);
    //fprintf(stderr, "%lu qualities have been decoding\n", qdecoded);
    if (todecode) {
        pdecode(blkbuff,
                prevbytes,
                todecode,
                buff,
                fqflag);
        //Check if newline is needed
        bleft = buffsize - todecode;
        if (bleft) {
            if (fqflag) {
                buff[todecode] = NL;
                bleft--;
                buff[todecode + 1] = '+';
                bleft--;
                buff[todecode + 2] = NL;
                bleft--;
                spacer = 3;
            }
        }
    }
    else {
        bleft = buffsize;
        if (!qdecoded) {
            //Add fastq separator if no qualities have been decoded
            buff[0] = NL;
            bleft--;
            buff[1] = '+';
            bleft--;
            buff[2] = NL;
            bleft--;
            spacer = 3;
        }
    }
    if (bleft && fqflag) {
        //qdecoded = seqlen < prevbytes ? prevbytes - seqlen - 3 : 0;
        qtodecode = seqlen - qdecoded < bleft ? seqlen - qdecoded : bleft;
        sqz_blkqdecode(blkbuff,                  //data buffer
                       buff + todecode + spacer, //output buffer
                       qtodecode,                //requested qualities
                       qdecoded);                 //already decoded qualities
    }
    return todecode + qtodecode + spacer;
}


uint64_t sqz_fastXdecode(sqzblock_t *blk,   //Data block
                         uint8_t *buff,     //Array to place decoded data
                         uint64_t buffsize, //Size of buff
                         char fqflag)       //0 - fasta, 1 - fastq
{
    uint8_t  *blkbuff   = blk->blkbuff;
    uint64_t  blkpos    = blk->blkpos;
    uint64_t  namepos   = blk->namepos;
    uint64_t  prevbytes = blk->prevlen;
    //Last 8 bytes store size of name buffer
    uint64_t  namesize  = *(uint64_t *)( blkbuff + ( blk->blksize - B64 ) );
    //Size of sequence data in block
    uint64_t  datasize  = blk->blksize - B64 - namesize;

    uint64_t  wbytes    = 0;
    char     *namebuff  = (char *)( blkbuff + datasize );
    //Read length of first sequence
    uint64_t seqlen     = *(uint64_t *)( blkbuff + blkpos );;
    //uint64_t todecode   = 0;
    //Length of name of first sequence
    uint64_t namelen    = strlen(namebuff + namepos);

    //Header values for fasta or fastq
    uint8_t H = FAH;
    uint8_t E = 3; //Number of new lines (fasta)
    if (fqflag) {  //(fastq)
        H = FQH;
        E = 6;
    }
    /*
      To fit current sequence, length of sequence bytes is needed (2x if fastq
      for quality values). Plus 3 bytes (6 for fastq) for header and 2 new lines
      (4 new lines + header + separator for fastq)
    */
    uint64_t buffneed = ( seqlen * (1 + fqflag) )  + (E + namelen);
    //Check if there is sequence that was not completely decoded make sure to
    //detect if only quality values need to be decoded.
    if (prevbytes) {
        prevbytes -= namelen + 2;
        wbytes = sqz_rdecode(blkbuff + blkpos + B64,
                             buff,
                             buffsize,
                             seqlen,
                             prevbytes,
                             fqflag);
        prevbytes += wbytes;
        if (prevbytes < buffneed - E - namelen)
            goto exit;
        buff[wbytes++] = NL;
        prevbytes = 0;
        blkpos += sqz_codeblksize(blkbuff + blkpos, fqflag);
        //fprintf(stderr, "From incomplete B\n");
        //Update buffer size
        //Get new sequence length, sequence name, sequence name lenght, needed
        //buffer size, and name position
        buffsize -= wbytes;
        //TODO: Overflow read when end of block has been reached
        //New sequence
        namepos += namelen + 1;
        namelen  = strlen(namebuff + namepos);
        seqlen   = *(uint64_t *)( blkbuff + blkpos );
        buffneed = ( seqlen * (1 + fqflag) )  + (E + namelen);
    }
    //While there is data to decode
    while ( blkpos < datasize ) {
        //Test if there is enough space for current sequence
        if ( buffsize < buffneed ) {
            //Test if at least: ">"|"@" + name length + '\n' + 1 base fits
            if (buffsize < namelen + 3) {
                prevbytes    = 0;
                blk->blkpos  = blkpos;
                blk->namepos = namepos;
                goto exit;
            }
            prevbytes = sqz_pdecode(blkbuff + blkpos + B64,
                                    buff + wbytes,
                                    buffsize,
                                    seqlen,
                                    namebuff + namepos,
                                    fqflag);
            wbytes += prevbytes;
            blk->blkpos  = blkpos;
            blk->namepos = namepos;
            goto exit;
        }
        blkpos += B64;
        buff[wbytes++] = H; //Add header symbol: ">" or "@"
        memcpy(buff + wbytes, namebuff + namepos, namelen); //Copy sequence name
        wbytes += namelen;
        buff[wbytes++] = NL;
        //Decode sequence and quality if fastq
        blkpos += sqz_seqdecode( blkbuff + blkpos,
                                 buff + wbytes,
                                 seqlen,
                                 fqflag,
                                 &wbytes);
        //fprintf(stderr, "From complete B\n");
        buffsize -= buffneed;
        namepos += namelen + 1;
        //New sequence
        seqlen   = *(uint64_t *)( blkbuff + blkpos );
        namelen  = strlen(namebuff + namepos);
        buffneed = ( seqlen * (1 + fqflag) )  + (E + namelen);
    }
    //fprintf(stderr, "Done with block\n");
    blk->newblk  = 0;
    blk->blkpos  = 0;
    blk->namepos = 0;
    exit:
        blk->prevlen = prevbytes;
        return wbytes;
}


uint64_t sqz_seqdecode(const uint8_t *codebuff,
                       uint8_t       *decodebuff,
                       uint64_t       length,
                       char           qflag,
                       uint64_t      *wbytes)
{
    uint64_t seqlen     = length;
    uint64_t buffpos    = 0;
    uint64_t blklen     = 0;
    uint64_t prevblk    = 0;
    uint64_t seqpos     = 0;
    uint8_t nflag;
    uint64_t nnum;
    uint64_t nbytes;
    while (length > 0) {
        blklen = *(uint64_t *)(codebuff + buffpos);
        buffpos += B64;
        buffpos += sqz_blkdecode(codebuff + buffpos,
                                 decodebuff + seqpos,
                                 &seqpos,
                                 blklen);
        length -= blklen;
        if (length) {     //Partial sequence decoding
            nflag = *(codebuff + buffpos);
            buffpos++;
            if (nflag) {  //Check if it is an N block
                nnum = countnblk(codebuff + buffpos, &nbytes);
                buffpos += nbytes;
                writens(nnum, decodebuff + seqpos);
                seqpos += nnum;
                length -= nnum;
            }
            else {
                if (qflag) {
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
                    buffpos += sqz_qualdecode(codebuff + buffpos,
                                              decodebuff + seqlen + 3,
                                              blklen);
                    //Change sequence length to reflect part of the quality
                    //string has been decoded
                    seqlen -= blklen;
                    //Store this blk length for proper placement of decoded
                    //qualities
                    prevblk = blklen;
                }
            }
            continue;
        }
        nflag = *(codebuff + buffpos);
        buffpos++;
    }
    //Decode quality strings for fastq data
    if (qflag) {
        decodebuff[seqpos++] = NL;
        decodebuff[seqpos++] = '+';
        decodebuff[seqpos++] = NL;
        buffpos += sqz_qualdecode(codebuff + buffpos,
                                  decodebuff + seqpos + prevblk,
                                  seqlen);
        seqpos += seqlen + prevblk;
    }
    decodebuff[seqpos++] = NL;
    *wbytes += seqpos;
    return buffpos;
}


size_t sqz_qualdecode(const uint8_t *codebuff,
                      uint8_t *qualstr,
                      uint64_t length)
{
    size_t byte = 0;
    size_t offset = 0;
    size_t decoded = 0;
    uint8_t code;
    unsigned char count;
    unsigned char q;
    while (decoded != length) {
        code = *(codebuff + byte); //get byte value
        count = code & 31;         //get current symbol count
        q = code >> 5;             //get symbol index value
        for (int i = 0; i <= count; i++) {
            qualstr[offset] = qual_val_table[q];
            offset++;
        }
        decoded += ++count;
        byte++;
    }
    return byte;
}


uint8_t sqz_8binqual(uint8_t q)
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


const uint8_t *sqz_findn(const uint8_t *seq)
{
    while (seq_nt4_tableSQZ[*seq] < 4) {
        seq++;
    }
    return seq;
}


uint64_t sqz_bit2encode(const uint8_t *seq,
                        size_t seqlen)
{
    uint64_t result = 0;
    for (uint64_t i = 0; i < seqlen; i++) {
        result = (result << 2) | seq_nt4_tableSQZ[seq[i]];
    }
    return result;
}



unsigned char sqz_bit2decode(const uint64_t *mer,
                             uint8_t *decoded,
                             unsigned char len)
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


unsigned char sqz_writens(uint8_t  numn,
                          uint8_t *decoded)
{
    unsigned char nwritten = 0;
    while (numn-- > 0) {
        decoded[numn] = 'N';
        nwritten++;
    }
    return nwritten;
}


static uint64_t sqz_blkcode(uint64_t *buff, const uint8_t *seq, uint64_t len)
{
    if (!len) return 0;
    uint64_t nblks = ( (len / 32) + ( (len % 32) > 0 ) ) - 1;
    uint64_t lenleft;
    uint64_t code;
    uint64_t buffpos;
    uint64_t nbytes = 0;

    //1 less than the total number of iterations needed.
    for (buffpos = 0; buffpos < nblks; buffpos++) {
        //Encode 32mer or what is left of sequence block
        code = sqz_bit2encode(seq, 32);
        memcpy(buff + buffpos, &code, B64);
        nbytes += B64;
        seq += 32;
    }
    lenleft = ( len - ( nblks * 32) );
    code = sqz_bit2encode(seq, lenleft);
    memcpy(buff + buffpos, &code, B64);
    buffpos++;
    nbytes += B64;
    return nbytes;
}


static uint64_t sqz_nblkcode(uint8_t *buff, uint64_t len)
{
    if (!len) return 0;
    //1 less than the total number of iterations needed.
    uint64_t nblks = ( (len / 127) + ( (len % 127) > 0 ) ) - 1;
    uint64_t buffpos = 0;
    uint8_t  l;
    for (buffpos = 0; buffpos < nblks; buffpos++)
        buff[buffpos] = 127; // 127 consecutive non ACGT bases

    //Last itteration needs special treatment
    //Set bit 7 if last N block before next base
    l = (len - (nblks * 127)) | NEND;
    //TODO this memcpy is dum buff[buffpos] = l;
    memcpy(buff + buffpos, &l, 1);
    buffpos++;
    return buffpos;
}


static uint64_t sqz_blkdecode(const uint8_t *codebuff,
                              uint8_t       *decodebuff,
                              uint64_t      *wbytes,
                              uint64_t      blklen)
{
    uint64_t codepos = 0;
    uint64_t decodepos = 0;
    uint64_t blknum = blklen / 32;;
    uint64_t blkleft  = blklen % 32;
    uint64_t *mer;
    uint64_t blkidx;
    mer = (uint64_t *)codebuff;
    for (blkidx = 0; blkidx < blknum; blkidx++) {
        sqz_bit2decode(mer + blkidx, decodebuff + decodepos, 32);
        //Keep traked of amount of decoded data: 64bit per mer
        codepos += B64;
        //Move sequence pointer by amount of bases decoded
        decodepos += 32;
    }
    if (blkleft) {
        sqz_bit2decode(mer + blkidx, decodebuff + decodepos, blkleft);
        codepos += B64;
        decodepos += blkleft;
    }
    *wbytes += decodepos;
    return codepos;
}
