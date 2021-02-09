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
    if (blk->newblk) return sqz_headblk(sqz, blk);
    return sqz_tailblk(sqz, blk);
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
    //TODO test if we can make due without this
    if (k != sqzsize) {
        fprintf(stderr, "[sqzlib ERROR]: Failed to unload raw data buffer\n");
        return 0;
    }
    //Indicate if there is more sequence to read
    //if (seqread < seqlen) {
    if (sqz->endflag) {
        //Unset new block flag. The rest of the sequence needs to be encoded
        //before a new block can be started.
        blk->newblk = 0;
        //Indicate how much sequence has been read
        //TODO check if this is usefull
        sqz->toread = seqread;
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


char sqz_tailblk(sqzfastx_t *sqz,
                 sqzblock_t *blk)
{
    uint64_t seqlen = sqz->prevlen;
    uint8_t *seq    = sqz->seqbuffer;
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


char sqz_fastaencode(sqzfastx_t *sqz, sqzblock_t *blk)
{
    if (blk->newblk) return sqz_fastaheadblk(sqz, blk);
    return sqz_fastatailblk(sqz, blk);
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
    while ( k < sqzsize ) {
        //Extract sequence length
        seqlen = *(uint64_t *)( seqbuffer + k );
        k += B64;
        //Store sequence length in code buffer
        memcpy(blkbuff + blkpos, &seqlen, B64);
        blkpos += B64;
        //Position seq and qual strings to corresponding place in buffers
        seq = seqbuffer + k;
        //Determine how much sequence is contained in the buffer
        //This happens because not all the sequence may be loaded in buffer
        seqread = seqlen < (sqzsize - k) ? seqlen : sqzsize - k - 1;
        k += seqread + 1;
        //encode sequence
        blkpos += sqz_seqencode(seq, seqread, blkbuff + blkpos, seqlen);
    }
    //TODO test if we can make due without this
    if (k != sqzsize) {
        fprintf(stderr, "[sqzlib ERROR]: Failed to unload raw data buffer\n");
        return 0;
    }
    //Indicate if there is more sequence to read
    if (sqz->endflag) {
        //Unset new block flag. The rest of the sequence needs to be encoded
        //before a new block can be started.
        blk->newblk = 0;
        //Indicate how much sequence has been read
        //TODO check if this is usefull
        sqz->toread = seqread;
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
	          while ( seq_nt4_tableSQZ[*npos] == 4) {
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


uint64_t sqz_fastXdecode(sqzblock_t *blk,
                         uint8_t *buff,
                         uint64_t size,
                         char fqflag)
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
    uint8_t H = FAH;
    if (fqflag) H = FQH;
    while ( blkpos < datasize ) {
        //Check how much sequence can be decoded
        seqlen = *(uint64_t *)( blkbuff + blkpos );
        namelen = strlen(namebuff + namepos);
        ///TODO awkward way of handling fq file
        /*
        fix bug when sequence is larger than buffer size
        When a sequence can not be completely decoded into a buffer, decoding
        is interupted so that buffers can be offloaded. Afterwards decoding
        continues from where it left. The problem arrises when a single sequence
        is larger than the provided buffer. In this scenario, the sequence is
        never loaded. This is not fatal during library execution as the decoding
        function exist as if it had decoded an entore block.
        */
        if ( (size - wbytes)  < ( seqlen * (1 + fqflag) )  + (6 + namelen ) ) {
            //if (seq)
            blk->blkpos  = blkpos;
            blk->namepos = namepos;
            goto exit;
        }
        blkpos += B64;
        //Add fastq header
        buff[wbytes++] = H;
        //Copy sequence name
        memcpy(buff + wbytes, namebuff + namepos, namelen);
        wbytes += namelen ;
        //Add newline
        buff[wbytes++] = NL;
        //Decode sequence and quality if fastq
        blkpos += sqz_seqdecode( blkbuff + blkpos,
                                 buff + wbytes,
                                 seqlen,
                                 fqflag,
                                 &wbytes);
        namepos += namelen + 1;
    }
    //Update blk on buffpos
    blk->blkpos  = 0;
    blk->namepos = 0;
    exit:
        return wbytes;
}


uint64_t sqz_seqdecode(const uint8_t *codebuff,
                       uint8_t       *decodebuff,
                       uint64_t       length,
                       char           qflag,
                       uint64_t      *wbytes)
{
    //TODO: Add feature to decode at most N bytes
    /*
      Currently, output buffer must be large enough to store entire fastX
      sequence. This has the limitation of allowing a maximum sequence length of
      LOAD_SIZE - 6 - sequence name size for fastA and
      (LOAD_SIZE - 6 - sequence name) / 2 for fastQ sequences
    */
    uint64_t seqlen     = length;
    uint64_t buffpos    = 0;
    uint64_t blklen     = 0;
    uint64_t prevblk    = 0;
    uint64_t seqpos     = 0;
    const uint8_t *nstr = NULL;
    while (length > 0) {
        blklen = *(uint64_t *)(codebuff + buffpos);
        buffpos += B64;
        if (blklen > length) blklen = length;
        buffpos += sqz_blkdecode(codebuff + buffpos,
                                 decodebuff + seqpos,
                                 &seqpos,
                                 blklen);
        length -= blklen;
        //Check how much sequence has been decoded
        if (length) {
            nstr = codebuff + buffpos;
            buffpos++;
            //Check if it is an N block
            if (*nstr) {
                //TODO move to function
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
                //TODO: This if might be unnnecessary
                if (qflag) {
                    /*
                    Even though entire sequence has not been decoded completely,
                    there is quality data to be decoded. This happens because in
                    this particular sequence, during the encoding process a buffer
                    was filled before entire sequence could be loaded. In this
                    case, what was loaded of the sequence is encoded, and then
                    the reast of the sequence finishes loading. Resulting in the
                    following patter in the code block:
                        blklen:seqcode:qualcode:blklen:seqcode:qualcode
                    Instead of:
                        blklen:seqcode:qualcode
                    */
                    buffpos += sqz_qualdecode(codebuff + buffpos,
                                              decodebuff + seqlen + 3,
                                              blklen,
                                              wbytes);
                    //Change sequence length to reflect part of the quality
                    //string has been decoded
                    seqlen -= blklen;
                    //Store this blk length for proper placement of decoded
                    //qualities
                    prevblk = blklen;
                }
            }
        }
    }
    //Decode quality strings for fastq data
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
        decoded[len] = seq_dec_tableSQZ[byte];
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


static uint64_t sqz_blkdecode(const uint8_t *codebuff,
                              uint8_t       *decodebuff,
                              uint64_t      *wbytes,
                              uint64_t      blklen)
{
    uint64_t codepos = 0;
    uint64_t decodepos = 0;
    uint64_t blknum;
    uint64_t *mer;
    uint64_t blkidx;
    blknum = blklen / 32;
    mer = (uint64_t *)codebuff;
    for (blkidx = 0; blkidx < blknum; blkidx++) {
        sqz_bit2decode(mer + blkidx, decodebuff + decodepos, 32);
        //Keep traked of amount of decoded data: 64bit per mer
        codepos += B64;
        //Move sequence pointer by amount of bases decoded
        decodepos += 32;
    }
    sqz_bit2decode(mer + blkidx, decodebuff + decodepos, blklen % 32);
    codepos += B64;
    decodepos += blklen % 32;
    *wbytes += decodepos;
    return codepos;
}
