#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include "sqz_coding.h"


char sqz_encode(sqzfastx_t *sqz, sqzblock_t *blk)
{
    if (blk->newblk) {
        fprintf(stderr, "[sqzlib INFO]: New block - %lu sequences.\n", sqz->n);
        return sqz_headblk(sqz, blk);
    }
    else {
         return sqz_tailblk(sqz, blk);
    }
}


char sqz_headblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    uint8_t *seq = NULL;
    uint8_t *qual = NULL;
    uint64_t seqlen = 0;
    size_t seqread;
    size_t lenbytes = sizeof(size_t);
    size_t k = 0;
    size_t n = 0;
    size_t sqzsize = sqz->offset;
    while ( k < sqzsize) {
        //Keep track of how many sequences have been encoded
        n++;
        //Extract sequence length
        seqlen = *(uint64_t *)(sqz->seqbuffer + k);
        //Keep track of number of bytes read from seqbuffer
        k += lenbytes;
        //Store sequence length in code buffer
        memcpy(blk->codebuff + blk->blksize, &seqlen, sizeof(uint64_t));
        //Keep track of number of bytes written to codebuff
        blk->blksize += sizeof(uint64_t);
        //Position seq and qual strings to corresponding place in buffers
        seq =  sqz->seqbuffer + k;
        qual = sqz->qualbuffer + k;
        //Determine how much sequence is contained in the buffer
        //This happens because not all the sequence may be loaded in buffer
        seqread = seqlen < (sqzsize - k) ? seqlen : sqzsize - k - 1;
        /*
        TODO - Problem, if only a fraction of the sequence was loaded, that
               number of bytes should be kept track of. If the entire sequence
               was extracted, and additional byte should be counted
               corresponding to the end of sequence NULL
        */
        k += seqread + 1;
        //encode sequence
        sqz_seqencode(seq, seqread, blk, seqlen);
        //do something with the qualities
        sqz_qualencode(qual, seqread, blk, seqlen);
    }
    //Indicate if there is more sequence to read
    if (seqread < seqlen) {
        //Unset new block flag. The rest of the sequence needs to be encoded
        //before a new block can be started.
        blk->newblk = 0;
        //Indicate how much sequence has been read
        sqz->toread = seqread;
        //Indicate length of sequence still needing loading
        //sqz->prevlen = seqlen;
    }
    if (k != sqzsize) {
        fprintf(stderr, "[sqzlib ERROR]: Failed to unload raw data buffer\n");
        fprintf(stderr, "\t\t%lu : %lu\n", k, sqzsize);
        return 0;
    }
    return 1;
}


char sqz_tailblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    //Check if reading leftover sequence
    size_t seqlen = sqz->prevlen;
    uint8_t *seq =  sqz->seqbuffer;
    uint8_t *qual = sqz->qualbuffer;
    //encode sequence
    sqz_seqencode(seq, sqz->rem, blk, seqlen - sqz->toread);
    //encode quality
    sqz_qualencode(qual, sqz->rem, blk, seqlen - sqz->toread);
    //Update how much sequence has been read
    sqz->toread += sqz->offset;

    //Indicate if there is more sequence to read
    if (sqz->toread < seqlen) {
        blk->newblk = 0;
    }
    else {
        blk->newblk = 1;
    }
    return 1;
}


char sqz_fastaencode(sqzfastx_t *sqz, sqzblock_t *blk)
{
    if (blk->newblk) {
        fprintf(stderr, "[sqzlib INFO]: New block - %lu sequences.\n", sqz->n);
        fprintf(stderr, "HEAD||||%d\n", sqz->endflag);
        return sqz_fastaheadblk(sqz, blk);
    }
    else {
        fprintf(stderr, "TAIL||||%d\n", sqz->endflag);
        return sqz_fastatailblk(sqz, blk);
    }
}


char sqz_fastaheadblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    uint8_t *seq = NULL;
    uint64_t seqlen = 0;
    size_t seqread;
    size_t lenbytes = sizeof(size_t);
    size_t k = 0;
    size_t n = 0;
    size_t sqzsize = sqz->offset;
    while ( k < sqzsize) {
        //Keep track of how many sequences have been encoded
        n++;
        //Extract sequence length
        seqlen = *(uint64_t *)(sqz->seqbuffer + k);
        //fprintf(stderr, "SEQLEN: %lu %lu\n", seqlen, n);
        //Keep track of number of bytes read from seqbuffer
        k += lenbytes;
        //Store sequence length in code buffer
        memcpy(blk->codebuff + blk->blksize, &seqlen, sizeof(uint64_t));
        //Keep track of number of bytes written to codebuff
        blk->blksize += sizeof(uint64_t);
        //Position seq and qual strings to corresponding place in buffers
        seq =  sqz->seqbuffer + k;
        //Determine how much sequence is contained in the buffer
        //This happens because not all the sequence may be loaded in buffer
        seqread = seqlen < (sqzsize - k) ? seqlen : sqzsize - k - 1;
        /*
        TODO - Problem, if only a fraction of the sequence was loaded, that
               number of bytes should be kept track of. If the entire sequence
               was extracted, and additional byte should be counted
               corresponding to the end of sequence NULL
        */
        k += seqread + 1;
        //encode sequence
        sqz_seqencode(seq, seqread, blk, seqlen);
    }
    //Indicate if there is more sequence to read
    if (seqread < seqlen) {
        //Unset new block flag. The rest of the sequence needs to be encoded
        //before a new block can be started.
        blk->newblk = 0;
        //Indicate how much sequence has been read
        sqz->toread = seqread;
        //Indicate length of sequence still needing loading
        //sqz->prevlen = seqlen;
    }
    if (k != sqzsize) {
        fprintf(stderr, "[sqzlib ERROR]: Failed to unload raw data buffer\n");
        fprintf(stderr, "\t\t%lu : %lu\n", k, sqzsize);
        return 0;
    }
    if (!sqz->endflag) {
        fprintf(stderr, "Adding names1 %lu\n", sqz->namelen);
        //memcpy(blk->codebuff + blk->blksize, sqz->namebuffer, sqz->namelen);
        //blk->blksize += sqz->namelen;
        //memcpy(blk->codebuff + blk->blksize, &(sqz->namelen), sizeof(size_t));
        //blk->blksize += sizeof(size_t);
        fprintf(stderr, "Adding names1 %lu\n", blk->blksize);
    }
    return 1;
}


char sqz_fastatailblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    //Check if reading leftover sequence
    size_t seqlen = sqz->prevlen;
    uint8_t *seq =  sqz->seqbuffer;
    //encode sequence
    sqz_seqencode(seq, sqz->rem, blk, seqlen - sqz->toread);
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
        fprintf(stderr, "Adding names2 %lu\n", sqz->namelen);
        //memcpy(blk->codebuff + blk->blksize, sqz->namebuffer, sqz->namelen);
        //blk->blksize += sqz->namelen;
        //memcpy(blk->codebuff + blk->blksize, &(sqz->namelen), sizeof(size_t));
        //blk->blksize += sizeof(size_t);
        fprintf(stderr, "Adding names2 %lu\n", blk->blksize);
    }
    return 1;
}


void sqz_seqencode(const uint8_t *seq,
                   uint64_t seqlen,
                   sqzblock_t *blk,
                   uint64_t seqlenOG)
{
    //TODO function to long and there is code repetition. Create new functions
    //for reduntsnt code
    //Set coding buffer to the next unused position in the block
    uint8_t *codebuff = blk->codebuff + blk->blksize;
    //Track position in sequence
    const uint8_t *lstop = seq;
    const uint8_t *nptr = seq + seqlen; //End of string
    size_t nn;                          //Number of Ns
    unsigned char wn;                   //127 N block. 1 bit flag 7 bit count
	  const uint8_t *npos;                //Track positions where N occurs
	  size_t blen = 0;                    //Length of segment before first N
    uint64_t code;                      //2 bit encoded sequence
    size_t wbytes = 0;                  //Number of bytes used to encode
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
            //fprintf(stderr, "Firstblock %s: %lu\n", npos, blen);
            //Determine number of consecutive Ns until next base
            nn = 0;
	          while ( seq_nt4_table[*npos] == 4) {
	              nn++;
                //fprintf(stderr, "inside ncount %s nn: %lu %u:\n",
                //        npos, nn, seq_nt4_table[*npos]);
		            npos++;
	          }
            //fprintf(stderr, "after ncount %s:\n", npos);
            //Write block length [blen]
            //fprintf(stderr, "||\t\tblklen: %lu\n", blen);
            memcpy(codebuff + wbytes, &blen, sizeof(size_t));
            wbytes += sizeof(size_t);
            //Loop over sequence block in chunks of 32 bases and encode
            //TODO Change loop to always fit 32mers, add final encoding of
            //whatever is left outside of loop (hint use modulo 32)
            for (size_t i = 0; i < blen; i+=32) {
                //Encode 32mer or what is left of sequence block
                code = sqz_bit2encode(lstop+i,((i+32)<= blen)?32:blen-i);
                //TODO change sizeof(uint64_t), by definition it's 8 bytes
                //so just use 8
                memcpy(codebuff + wbytes, &code, sizeof(uint64_t));
                wbytes += sizeof(uint64_t);
                nbases += ((i+32)<= blen) ? 32 : blen - i;
            }
            //Indicate that an N block follows
            memcpy(codebuff + wbytes, &flag1, 1);
            wbytes++;
            //Loop over Ns in chunks of 127 and encode
            //TODO move to function
            for (size_t i = 0; i < nn; i +=127) {
                //Encode 127 Ns or what is left
                wn = ((i+127) <= nn)?127:nn-i;
                //Set bit 7 if last N block before next base
                if (i+127 >= nn) wn = wn | NEND;
                memcpy(codebuff + wbytes, &wn, 1);
                wbytes++;
                nbases += ((i+127) <= nn)?127:nn-i;
            }
            //Trace next base after sequence of Ns
            lstop = npos;
	      }
    } while (*npos);
    //Detect and encode trailing bases
    blen = nptr - lstop;
    if (blen) {
        memcpy(codebuff + wbytes, &blen, sizeof(size_t));
        wbytes += sizeof(size_t);
        for (uint32_t i = 0; i < blen; i+=32) {
            code = sqz_bit2encode(lstop + i, ((i+32) <= blen) ? 32 : blen - i);
            memcpy(codebuff + wbytes, &code, sizeof(uint64_t));
            wbytes += sizeof(uint64_t);
            nbases += ((i+32) <= blen) ? 32 : blen - i;
        }
    }
    //Check that no more sequence needs to be loaded. If more sequence is needed
    //write corresponding byte
    if (seqlenOG - nbases) {
        memcpy(codebuff + wbytes, &flag2, 1);
        wbytes++;
    }
    //Move offset by number of bytes written
    blk->blksize += wbytes;

}


size_t sqz_qualencode(const uint8_t *qual,
                      size_t quallen,
                      sqzblock_t *blk,
                      uint64_t seqlen)
{
    uint8_t *codebuff = blk->codebuff + blk->blksize;
    size_t bytes = 0;
    uint8_t q = sqz_8binqual(*qual);
    uint8_t c = 0;
    uint8_t code = 0;
    while(*qual) {
        if ( sqz_8binqual(*(qual+1)) == q) {
            c++;
            if (c == 31) {
                //Encode
                code = code | q;
                code = code | c;
                memcpy(codebuff + bytes, &code, 1);
                bytes += 1;
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
        memcpy(codebuff + bytes, &code, 1);
        qual++;
        q = sqz_8binqual(*qual);
        c = 0;
        bytes += 1;
        code = 0;
    }
    blk->blksize += bytes;
    return bytes;
}


size_t sqz_fastqdecode(const uint8_t *buff, size_t size)
{
    char ret = 0;
    size_t buffpos = 0;
    size_t seqlen;
    char *seqstr;
    char *qualstr;
    uint64_t n = 0;
    while (buffpos < size) {
        seqlen = *(size_t *)(buff + buffpos);
        //TODO change to constant
        buffpos += sizeof(size_t);
        seqstr = malloc(seqlen + 1);
        qualstr = malloc(seqlen + 1);
        if (!seqstr | !qualstr) {
            fprintf(stderr, "[sqzlib ERROR]: Insufficient memory.\n");
            goto exit;
        }
        seqstr[seqlen] = '\0';
        qualstr[seqlen] = '\0';
        /*
          TODO - Problem, when a sequence is only partially loaded. The
          sequence code is followed by the quality code bedofre the the entire
          sequence information is stored. This is not taken into account in the
          decoding logic
        */

        buffpos += sqz_seqdecode(buff + buffpos, seqstr, qualstr, seqlen);
        //buffpos += sqz_qualdecode(buff + buffpos, qualstr, seqlen);
        fprintf(stdout, "@\n%s\n+\n", seqstr);
        fprintf(stdout, "%s\n", qualstr);
        free(seqstr);
        free(qualstr);
        n++;
    }
    ret = 1;
    fprintf(stderr, "\t%lu sequences in block\n", n);
    exit:
        return ret;
}


size_t sqz_seqdecode(const uint8_t *codebuff,
                     char *seqstr,
                     char *qualstr,
                     size_t length)
{
    size_t buffpos = 0;
    uint64_t blklen;
    uint64_t *mer;
    size_t blkidx;
    char merlen;
    size_t seqpos = 0;
    size_t qualpos = 0;
    const uint8_t *nstr;
    while (length > 0) {
        //blen - length of block (within sequence up to first non acgt base)
        blklen = *(uint64_t *)(codebuff + buffpos);
        buffpos += sizeof(uint64_t);
        //TODO Define a constante equal to 8
        //buffpos records number of bytes read, but when used as an array offset,
        //it moves as many positions as the type. So if type is uint64_t(8bytes)
        //adding buffpos as ofset moves buffpos*8 bytes
        //TODO try to remove this cast
        mer = (uint64_t *)(codebuff + buffpos);
        blkidx = 0;
        //TODO, compute number of iterations instead of increments of 32
        //TODO move to function
        for (size_t i = 0; i < blklen; i+=32) {
            //Blocks code 32mers except last block which may code a shorter kmer
            merlen = ((i+32)<= blklen)?32:blklen-i;
            sqz_bit2decode(mer + blkidx, seqstr + seqpos, merlen);
            //Keep traked of amount of decoded data: 64bit per mer
            //Same as before
            buffpos += sizeof(uint64_t);
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
                    sqz_writens(*nstr & ~(1<<7), seqstr + seqpos);
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
                                              qualstr + qualpos,
                                              blklen);
                    qualpos += blklen;
                }
            }
        }

    }
    if (qualstr)
        buffpos += sqz_qualdecode(codebuff + buffpos,
                                  qualstr + qualpos,
                                  seqpos  - qualpos);
    return buffpos;
}


size_t sqz_fastadecode(const uint8_t *buff, size_t size)
{
    char ret = 0;
    size_t buffpos = 0;
    size_t seqlen;
    char *seqstr;
    uint64_t n = 0;
    while (buffpos < size) {
        seqlen = *(size_t *)(buff + buffpos);
        //TODO change to constant
        buffpos += sizeof(size_t);
        seqstr = malloc(seqlen + 1);
        if (!seqstr) {
            fprintf(stderr, "[sqzlib ERROR]: Insufficient memory.\n");
            goto exit;
        }
        seqstr[seqlen] = '\0';
        /*
          TODO - Problem, when a sequence is only partially loaded. The
          sequence code is followed by the quality code bedofre the the entire
          sequence information is stored. This is not taken into account in the
          decoding logic
        */

        buffpos += sqz_seqdecode(buff + buffpos, seqstr, NULL, seqlen);
        //buffpos += sqz_qualdecode(buff + buffpos, qualstr, seqlen);
        fprintf(stdout, ">\n%s\n", seqstr);
        //fprintf(stdout, "%s\n", qualstr);
        free(seqstr);
        n++;
    }
    ret = 1;
    fprintf(stderr, "\t%lu sequences in block\n", n);
    exit:
        return ret;
}





size_t sqz_qualdecode(const uint8_t *codebuff, char *qualstr, size_t length)
{
    size_t byte = 0;
    size_t offset = 0;
    size_t decoded = 0;
    uint8_t code;
    unsigned char count;
    unsigned char q;
    uint64_t total = 0;
    while (decoded != length) {
        code = *(codebuff + byte); //get byte value
        count = code & 31;         //get current symbol count
        q = code >> 5;             //get symbol index value
        total += count + 1;
        //fprintf(stderr, "\tsymbol: %d count: %d\n", q, count);
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


uint64_t sqz_bit2encode(const uint8_t *seq, size_t seqlen)
{
    uint64_t result = 0;
    for (size_t i = 0; i < seqlen; i++) {
        result = (result << 2) | seq_nt4_table[seq[i]];
    }
    return result;
}


sqzblock_t *sqz_sqzblkinit(size_t size)
{
    sqzblock_t *blk = malloc(sizeof(sqzblock_t));
    if (!blk) return NULL;
    //Code buffer array
    blk->codebuff = malloc(2*size);
    if (!blk->codebuff) {
        fprintf(stderr, "[libsqz ERROR]: memory error.\n");
        free(blk);
        return NULL;
    }
    blk->blksize = 0;
    blk->newblk = 1;
    //Compression buffer array
    blk->cmpbuff = malloc(2*size);
    if (!blk->cmpbuff) {
        fprintf(stderr, "[libsqz ERROR]: memory error.\n");
        free(blk->codebuff);
        free(blk);
        return NULL;
    }
    blk->cmpsize = 2*size;
    return blk;
}


void sqz_blkdestroy(sqzblock_t *blk)
{
    if (blk) {
        free(blk->codebuff);
        free(blk->cmpbuff);
        free(blk);
    }
}


unsigned char sqz_bit2decode(const uint64_t *mer,
                             char *decoded,
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


unsigned char sqz_writens(unsigned char numn, char *decoded)
{
    unsigned char nwritten = 0;
    while (numn-- > 0) {
        decoded[numn] = 'N';
        nwritten++;
    }
    return nwritten;
}


