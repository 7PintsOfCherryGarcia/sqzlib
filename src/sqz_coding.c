#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>

#define SQZLIB
#define KLIB
#include "sqz_coding.h"


char sqz_fastqencode(sqzfastx_t *sqz,
                     sqzblock_t *blk)
{
    if (blk->newblk) {
        return sqz_headblk(sqz, blk);
    }
    else {
        return sqz_tailblk(sqz, blk);
    }
}


char sqz_headblk(sqzfastx_t *sqz,
                 sqzblock_t *blk)
{
    uint8_t  *blkbuff    = blk->blkbuff;
    uint64_t  blkpos     = blk->blkpos;
    uint8_t  *seqbuffer  = sqz->seqbuffer;
    uint8_t  *qualbuffer = sqz->qualbuffer;
    uint64_t  sqzsize    = sqz->offset;

    uint8_t *seq    = NULL;
    uint8_t *qual   = NULL;
    uint64_t seqlen = 0;
    uint64_t seqread;
    uint64_t k = 0;

    while ( k < sqzsize ) {
        //Extract sequence length
        seqlen = *(uint64_t *)( seqbuffer + k );
        k += B64;
        //Store sequence length
        memcpy(blkbuff + blkpos, &seqlen, B64);
        blkpos += B64;
        //Position seq and qual strings to corresponding place in buffers
        seq  = seqbuffer + k;
        qual = qualbuffer + k;
        //Determine how much sequence is contained in the buffer
        //This happens because not all the sequence may be loaded in buffer
        seqread = seqlen < (sqzsize - k) ? seqlen : sqzsize - k - 1;
        k += seqread + 1;
        //encode sequence
        blkpos += sqz_seqencode(seq, seqread, blkbuff + blkpos, seqlen);
        //encode qualities
        blkpos += sqz_qualencode(qual, seqread, blkbuff + blkpos, seqlen);
    }
    if (k != sqzsize) {
        fprintf(stderr, "[sqzlib ERROR]: Failed to unload raw data buffer\n");
        return 0;
    }
    //Indicate if there is more sequence to read
    if (seqread < seqlen) {
        //Unset new block flag. The rest of the sequence needs to be encoded
        //before a new block can be started.
        blk->newblk = 0;
        //Indicate how much sequence has been read
        sqz->toread = seqread;
    }
    if (!sqz->endflag) {
        memcpy(blkbuff + blkpos, sqz->namebuffer, sqz->namepos);
        blkpos += sqz->namepos;
        memcpy(blkbuff + blkpos, &(sqz->namepos), B64);
        blkpos += B64;
    }
    blk->blkpos = blkpos;
    return 1;
}


char sqz_tailblk(sqzfastx_t *sqz,
                 sqzblock_t *blk)
{
    uint64_t seqlen = sqz->prevlen;
    uint8_t *seq    =  sqz->seqbuffer;
    uint8_t *qual   = sqz->qualbuffer;

    uint8_t *blkbuff = blk->blkbuff;
    uint64_t blksize = blk->blksize;
    uint64_t blkpos  = blk->blkpos;
    uint64_t seqread = sqz->rem;
    uint64_t seqleft = seqlen - sqz->toread;
    //encode sequence
    blkpos += sqz_seqencode(seq, seqread, blkbuff + blkpos, seqleft);
    //encode quality
    blkpos += sqz_qualencode(qual, seqread, blkbuff + blkpos, seqleft);
    //Update how much sequence has been read
    sqz->toread += sqz->offset;

    //Indicate if there is more sequence to read
    if (sqz->toread < seqlen) {
        blk->newblk = 0;
    }
    else {
        blk->newblk = 1;
    }
    if (!sqz->endflag) {
        //Copy name data making sure it fits
        if ( (blkpos + sqz->namepos) >= blksize ) {
            fprintf(stderr, "REALLOCING\n");
            blksize += sqz->namepos + B64;
            blkbuff = realloc(blk->blkbuff, blksize);
            blk->blkbuff = blkbuff;
            blk->blksize = blksize;
        }
        memcpy(blkbuff + blkpos, sqz->namebuffer, sqz->namepos);
        blkpos += sqz->namepos;
        memcpy(blkbuff + blkpos, &(sqz->namepos), B64);
        blkpos += B64;
    }
    blk->blkpos = blkpos;
    return 1;
}


char sqz_fastaencode(sqzfastx_t *sqz, sqzblock_t *blk)
{
    if (blk->newblk) {
        return sqz_fastaheadblk(sqz, blk);
    }
    else {
        return sqz_fastatailblk(sqz, blk);
    }
}


char sqz_fastaheadblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    uint8_t *blkbuff   = blk->blkbuff;
    uint64_t blkpos    = blk->blkpos;
    uint8_t *seqbuffer = sqz->seqbuffer;
    uint64_t sqzsize   = sqz->offset;

    uint8_t *seq = NULL;
    uint64_t seqlen = 0;
    uint64_t seqread;
    uint64_t k = 0;

    //
    char *namebuff = (char *)(sqz->namebuffer);
    //

    while ( k < sqzsize ) {
        //Extract sequence length
        seqlen = *(uint64_t *)( seqbuffer + k );

        //Keep track of number of bytes read from seqbuffer
        k += B64;
        //Store sequence length in code buffer
        memcpy(blkbuff + blkpos, &seqlen, B64);
        blkpos += B64;
        //Position seq and qual strings to corresponding place in buffers
        seq = seqbuffer + k;
        namebuff += strlen(namebuff) + 1;
        //Determine how much sequence is contained in the buffer
        //This happens because not all the sequence may be loaded in buffer
        seqread = seqlen < (sqzsize - k) ? seqlen : sqzsize - k - 1;
        k += seqread + 1;
        //encode sequence
        blkpos += sqz_seqencode(seq, seqread, blkbuff + blkpos, seqlen);
    }
    if (k != sqzsize) {
        fprintf(stderr, "[sqzlib ERROR]: Failed to unload raw data buffer\n");
        return 0;
    }
    //Indicate if there is more sequence to read
    if (seqread < seqlen) {
        //Unset new block flag. The rest of the sequence needs to be encoded
        //before a new block can be started.
        blk->newblk = 0;
        //Indicate how much sequence has been read
        sqz->toread = seqread;
    }
    if (!sqz->endflag) {
        memcpy(blkbuff + blkpos, sqz->namebuffer, sqz->namepos);
        blkpos += sqz->namepos;
        memcpy(blkbuff + blkpos, &(sqz->namepos), B64);
        blkpos += B64;
    }
    blk->blkpos = blkpos;
    return 1;
}


char sqz_fastatailblk(sqzfastx_t *sqz,
                      sqzblock_t *blk)
{
    uint64_t seqlen = sqz->prevlen;
    uint8_t *seq    = sqz->seqbuffer;

    uint8_t *blkbuff = blk->blkbuff;
    uint64_t blksize = blk->blksize;
    uint64_t blkpos  = blk->blkpos;
    uint64_t seqread = sqz->rem;
    uint64_t seqleft = seqlen - sqz->toread;

    //encode sequence
    blkpos += sqz_seqencode(seq, seqread, blkbuff + blkpos, seqleft);
    //Update how much sequence has been read
    sqz->toread += sqz->offset;

    //Indicate if there is more sequence to read
    if (sqz->toread < seqlen) {
        blk->newblk = 0;
    }
    else {
        blk->newblk = 1;
    }
    if (!sqz->endflag) {
        if ( (blkpos + sqz->namepos) >= blksize ) {
            fprintf(stderr, "REALLOCING\n");
            blksize += sqz->namepos + B64;
            blkbuff = realloc(blk->blkbuff, blksize);
            blk->blkbuff = blkbuff;
            blk->blksize = blksize;
        }
        memcpy(blkbuff + blkpos, sqz->namebuffer, sqz->namepos);
        blkpos += sqz->namepos;
        memcpy(blkbuff + blkpos, &(sqz->namepos), B64);
        blkpos += B64;
    }
    blk->blkpos = blkpos;
    return 1;
}


uint64_t sqz_seqencode(const uint8_t *seq,
                       uint64_t seqlen,
                       uint8_t *blkbuff,
                       uint64_t seqlenOG)
{
    //TODO function to long and there is code repetition. Create new functions
    //for reduntsnt code
    //Set coding buffer to the next unused position in the block
    uint64_t blkpos = 0;
    //Track position in sequence
    const uint8_t *lstop = seq;
    const uint8_t *nptr = seq + seqlen; //End of string
    size_t nn;                          //Number of Ns
    unsigned char wn;                   //127 N block. 1 bit flag 7 bit count
	  const uint8_t *npos;                //Track positions where N occurs
	  uint64_t blen = 0;                    //Length of segment before first N
    uint64_t code;                      //2 bit encoded sequence
    uint64_t nbases = 0;                //Number of bases encoded
    unsigned char flag1 = 255;
    unsigned char flag2 = 0;
    //Main encoding loop
    do {
        //Get position of first non ACGTacgt base
	      npos = sqz_findn(lstop);
        if (*npos) {
            //Determine block length up to found N
            blen = npos - lstop;
            //Determine number of consecutive Ns until next base
            nn = 0;
	          while ( seq_nt4_table[*npos] == 4) {
	              nn++;
		            npos++;
	          }
            //Write block length [blen]
            memcpy(blkbuff + blkpos, &blen, B64);
            blkpos += B64;
            //Loop over sequence block in chunks of 32 bases and encode
            //TODO Change loop to always fit 32mers, add final encoding of
            //whatever is left outside of loop (hint use modulo 32)
            for (uint64_t i = 0; i < blen; i+=32) {
                //Encode 32mer or what is left of sequence block
                code = sqz_bit2encode(lstop+i,((i+32)<= blen)?32:blen-i);
                //TODO change sizeof(uint64_t), by definition it's 8 bytes
                //so just use 8
                memcpy(blkbuff + blkpos, &code, sizeof(uint64_t));
                blkpos += B64;
                nbases += ( (i+32) <= blen ) ? 32 : blen - i;
            }
            //Indicate that an N block follows
            memcpy(blkbuff + blkpos, &flag1, 1);
            blkpos++;
            //Loop over Ns in chunks of 127 and encode
            //TODO move to function
            for (size_t i = 0; i < nn; i +=127) {
                //Encode 127 Ns or what is left
                wn = ((i+127) <= nn)?127:nn-i;
                //Set bit 7 if last N block before next base
                if (i+127 >= nn) wn = wn | NEND;
                memcpy(blkbuff + blkpos, &wn, 1);
                blkpos++;
                nbases += ((i+127) <= nn)?127:nn-i;
            }
            //Trace next base after sequence of Ns
            lstop = npos;
	      }
    } while (*npos);
    //Detect and encode trailing bases
    blen = nptr - lstop;
    if (blen) {
        memcpy(blkbuff + blkpos, &blen, B64);
        blkpos += B64;
        for ( uint64_t i = 0; i < blen; i+=32 ) {
            code = sqz_bit2encode(lstop + i, ((i+32) <= blen) ? 32 : blen - i);
            memcpy(blkbuff + blkpos, &code, B64);
            blkpos += B64;
            nbases += ((i+32) <= blen) ? 32 : blen - i;
        }
    }
    //Check that no more sequence needs to be loaded. If more sequence is needed
    //write corresponding byte
    if (seqlenOG - nbases) {
        memcpy(blkbuff + blkpos, &flag2, 1);
        blkpos++;
    }
    //Move offset by number of bytes written
    //blk->blksize += wbytes;
    return blkpos;
}


uint64_t sqz_qualencode(const uint8_t *qual,
                        uint64_t quallen,
                        uint8_t *blkbuff,
                        uint64_t seqlen)
{
    //uint8_t *codebuff = blk->blkbuff + blk->blksize;
    uint64_t blkpos = 0;
    uint8_t q       = sqz_8binqual(*qual);
    uint8_t c       = 0;
    uint8_t code    = 0;
    while(*qual) {
        if ( sqz_8binqual(*(qual+1)) == q) {
            c++;
            if (c == 31) {
                //Encode
                code = code | q;
                code = code | c;
                memcpy(blkbuff + blkpos, &code, 1);
                blkpos += 1;
                c = 0;
                code = 0;
                q = sqz_8binqual(*(++qual));
            }
            qual++;
            continue;
        }
        //Encode
        code = code | q;
        code = code | c;
        memcpy(blkbuff + blkpos, &code, 1);
        qual++;
        q = sqz_8binqual(*qual);
        c = 0;
        blkpos += 1;
        code = 0;
    }
    return blkpos;
}


size_t sqz_fastqdecode(sqzblock_t *blk)
{
    uint64_t size = blk->blksize;
    uint8_t *buff = blk->blkbuff;
    uint64_t blkpos = 0;
    size_t seqlen;
    char *seqstr;
    char *qualstr;
    uint64_t n = 0;
    uint64_t namelen = *(uint64_t *)(buff + (size - B64));
    char *namebuff = (char *)(buff + (size - B64 - namelen));
    while ( blkpos < (size - B64 - namelen) ) {
        seqlen = *(size_t *)(buff + blkpos);
        blkpos += B64;
        //TODO get rid of these mallocs
        seqstr = malloc(seqlen + 1);
        qualstr = malloc(seqlen + 1);
        seqstr[seqlen] = '\0';
        qualstr[seqlen] = '\0';
        blkpos += sqz_seqdecode(buff + blkpos, seqstr, qualstr, seqlen);
        fprintf(stdout, "@%s\n%s\n+\n%s\n", namebuff, seqstr, qualstr);
        namebuff += strlen(namebuff) + 1;
        free(seqstr);
        free(qualstr);
        n++;
    }
    return 1;
}


size_t sqz_seqdecode(const uint8_t *codebuff,
                     char *seqstr,
                     char *qualstr,
                     size_t length)
{
    uint64_t buffpos = 0;
    uint64_t blklen;
    uint64_t *mer;
    size_t blkidx;
    char merlen;
    size_t seqpos = 0;
    size_t qualpos = 0;
    const uint8_t *nstr;
    //tmp
    uint64_t wbytes = 0;
    while (length > 0) {
        //blen - length of block (within sequence up to first non acgt base)
        blklen = *(uint64_t *)(codebuff + buffpos);
        buffpos += B64;
        //buffpos records number of bytes read, but when used as an array offset,
        //it moves as many positions as the type. So if type is uint64_t(8bytes)
        //adding buffpos as ofset moves buffpos*8 bytes
        mer = (uint64_t *)(codebuff + buffpos);
        blkidx = 0;
        //TODO, compute number of iterations instead of increments of 32
        //TODO move to function
        for (size_t i = 0; i < blklen; i+=32) {
            //Blocks code 32mers except last block which may code a shorter kmer
            merlen = ((i+32)<= blklen)?32:blklen-i;
            sqz_bit2decode(mer + blkidx, (uint8_t *)(seqstr + seqpos), merlen);
            //Keep traked of amount of decoded data: 64bit per mer
            //Same as before
            buffpos += B64;
            //Move sequencepointerbyabount of bases decoded (32 except last mer)
            seqpos += merlen;
            //Track next block to decode
            blkidx++;
        }
        length -= blklen;

        //Check how much sequence has been decoded
        if (length) {
            nstr = codebuff + buffpos;
            buffpos++;
            //Check if it is an N block
            if (*nstr) {
                nstr = codebuff + buffpos;
                while (1) {
                    unsigned char numn = *nstr & ~(1<<7);
                    sqz_writens(*nstr & ~(1<<7), (uint8_t *)(seqstr + seqpos));
                    seqpos += numn;
                    length -= numn;
                    buffpos++;
                    if (*nstr & 128)
                        break;
                    nstr++;
                }
            }
            else {
                //Otherwise, it's a quality block of blklen
                if (qualstr) {
                    buffpos += sqz_qualdecode(codebuff + buffpos,
                                              (uint8_t *)(qualstr + qualpos),
                                              blklen,
                                              &wbytes);
                    qualpos += blklen;
                }
            }
        }
    }
    if (qualstr)
        buffpos += sqz_qualdecode(codebuff + buffpos,
                                  (uint8_t *)(qualstr + qualpos),
                                  seqpos  - qualpos,
                                  &wbytes);
    return buffpos;
}


uint64_t sqz_fastXdecode(sqzblock_t *blk,
                         uint8_t *klibbuff,
                         uint64_t size,
                         char fastqflag)
{
    uint8_t *blkbuff  = blk->blkbuff;
    uint64_t blksize  = blk->blksize;
    uint64_t blkpos   = blk->blkpos;
    uint64_t namepos  = blk->namepos;
    uint64_t namesize = *(uint64_t *)( blkbuff + ( blksize - B64 ) );
    uint64_t datasize = blksize - B64 - namesize;
    uint64_t wbytes   = 0;
    uint64_t seqlen   = 0;
    uint64_t namelen  = 0;
    char *namebuff = (char *)( blkbuff + datasize );
    while ( blkpos < datasize ) {
        //Check how much sequence can be decoded
        //TODO add this computation to the while condition
        seqlen = *(uint64_t *)( blkbuff + blkpos );
        if ( (size - wbytes)  < ( seqlen * 2 )  + (6 + namelen ) ) {
            blk->blkpos  = blkpos;
            blk->namepos = namepos;
            goto exit;
        }
        namelen = strlen(namebuff + namepos);
        blkpos += B64;

        //Add fastq header
        klibbuff[wbytes++] = FQH;
        //Copy sequence name
        memcpy(klibbuff + wbytes, namebuff + namepos, namelen);
        wbytes += namelen ;
        //Add newline
        klibbuff[wbytes++] = NL;
        //Decode sequence and quality if fastq
        blkpos += sqz_seqdecode2( blkbuff + blkpos,
                                  klibbuff + wbytes,
                                  seqlen,
                                  fastqflag,
                                  &wbytes);
        namepos += namelen + 1;
    }
    //Update blk on buffpos
    blk->blkpos  = 0;
    blk->blksize = 0;
    blk->namepos = 0;
    exit:
        return wbytes;
}


uint64_t sqz_seqdecode2(const uint8_t *codebuff,
                        uint8_t       *decodebuff,
                        uint64_t       length,
                        char           qflag,
                        uint64_t      *wbytes)
{
    uint64_t seqlen     = length;
    uint64_t buffpos    = 0;
    uint64_t blklen     = 0;
    uint64_t prevblk    = 0;
    uint64_t *mer       = NULL;
    uint64_t blkidx     = 0;
    char merlen         = 0;
    uint64_t seqpos     = 0;
    const uint8_t *nstr = NULL;
    while (length > 0) {
        blklen = *(uint64_t *)(codebuff + buffpos);
        buffpos += B64;
        mer = (uint64_t *)(codebuff + buffpos);
        blkidx = 0;
        //TODO, compute number of iterations instead of increments of 32
        //TODO move to function
        for (size_t i = 0; i < blklen; i+=32) {
            //Blocks code 32mers except last block which may code a shorter kmer
            merlen = ((i+32)<= blklen)?32:blklen-i;
            sqz_bit2decode(mer + blkidx, decodebuff + seqpos, merlen);
            //Keep traked of amount of decoded data: 64bit per mer
            //Same as before
            buffpos += B64;
            //Move sequencepointerbyabount of bases decoded (32 except last mer)
            seqpos += merlen;
            //Track next block to decode
            blkidx++;
        }
        length -= blklen;
        //Check how much sequence has been decoded
        if (length) {
            nstr = codebuff + buffpos;
            buffpos++;
            //Check if it is an N block
            if (*nstr) {
                nstr = codebuff + buffpos;
                while (1) {
                    unsigned char numn = *nstr & ~(1<<7);
                    sqz_writens(*nstr & ~(1<<7), decodebuff + seqpos);
                    seqpos += numn;
                    length -= numn;
                    buffpos++;
                    if (*nstr & 128)
                        break;
                    nstr++;
                }
            }
            else {
                //Otherwise, it's a quality block of blklen
                if (qflag) {
                    /*
                    Even though entore sequence has not been decoded completely,
                    there is quality data to be decoded. This happens because in
                    this particular sequence, a buffer was filled before entire
                    sequence could be loaded. In this case, what was loaded of
                    the sequence is encoded, and then the reast of the sequence
                    finishes loading. Resulting in the following patter in the
                    code block:
                        seqcode:qualcode:seqcode:qualcode
                    Instead of:
                        seqcode:qualcode
                    */
                    buffpos += sqz_qualdecode(codebuff + buffpos,
                                              decodebuff + seqlen + 3,
                                              blklen,
                                              wbytes);
                    seqlen -= blklen;
                    prevblk = blklen;
                }
            }
        }
    }

    if (qflag) {
        decodebuff[seqpos++] = NL;
        decodebuff[seqpos++] = '+';
        decodebuff[seqpos++] = NL;
        buffpos += sqz_qualdecode(codebuff + buffpos,
                                  decodebuff + seqpos + prevblk,
                                  seqlen,
                                  wbytes);
        seqpos += seqlen + prevblk;
    }
    decodebuff[seqpos++] = NL;
    *wbytes += seqpos;
    return buffpos;
}


size_t sqz_fastadecode(const uint8_t *buff,
                       size_t size)
{
    char ret = 0;
    size_t buffpos = 0;
    size_t seqlen;
    char *seqstr;
    uint64_t namelen = *(uint64_t *)(buff + (size - 8));
    char *namebuff = (char *)(buff + (size - 8 - namelen));
    while (buffpos < (size - sizeof(uint64_t) - namelen) ) {
        seqlen = *(size_t *)(buff + buffpos);
        //TODO change to constant
        buffpos += sizeof(size_t);
        seqstr = malloc(seqlen + 1);
        if (!seqstr) {
            fprintf(stderr, "[sqzlib ERROR]: Insufficient memory.\n");
            goto exit;
        }
        seqstr[seqlen] = '\0';
        buffpos += sqz_seqdecode(buff + buffpos, seqstr, NULL, seqlen);
        fprintf(stdout, ">%s\n%s\n",namebuff, seqstr);
        namebuff += strlen(namebuff) + 1;
        free(seqstr);
    }
    ret = 1;
    exit:
        return ret;
}


size_t sqz_qualdecode(const uint8_t *codebuff,
                      uint8_t *qualstr,
                      uint64_t length,
                      uint64_t *wbytes)
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
    //TODO bottom bound bin instead of upper bin
    if (q >= 73) return 7<<5;
    if ( (q < 73) && (q >= 68) ) return 6<<5;
    if ( (q < 68) && (q >= 63) ) return 5<<5;
    if ( (q < 63) && (q >= 58) ) return 4<<5;
    if ( (q < 58) && (q >= 53) ) return 3<<5;
    if ( (q < 53) && (q >= 43) ) return 2<<5;
    if ( (q < 43) && (q >= 35) ) return 1<<5;
    return 0;
}


const uint8_t *sqz_findn(const uint8_t *seq)
{
    while (seq_nt4_table[*seq] < 4) {
        seq++;
    }
    return seq;
}


uint64_t sqz_bit2encode(const uint8_t *seq,
                        size_t seqlen)
{
    uint64_t result = 0;
    for (size_t i = 0; i < seqlen; i++) {
        result = (result << 2) | seq_nt4_table[seq[i]];
    }
    return result;
}


void sqz_blkdestroy(sqzblock_t *blk)
{
    if (blk) {
        free(blk->blkbuff);
        free(blk->cmpbuff);
        free(blk);
    }
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
        decoded[len] = seq_dec_table[byte];
        nbase++;
        code >>= 2;
    } while (len-- != 0);
    return nbase;
}


unsigned char sqz_writens(unsigned char numn,
                          uint8_t *decoded)
{
    unsigned char nwritten = 0;
    while (numn-- > 0) {
        decoded[numn] = 'N';
        nwritten++;
    }
    return nwritten;
}


