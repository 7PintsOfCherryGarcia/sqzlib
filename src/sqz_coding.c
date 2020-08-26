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
    size_t *seqlen = NULL;
    size_t seqread;
    size_t lenbytes = sizeof(size_t);
    size_t k = 0;
    size_t n = 0;
    size_t sqzsize = sqz->offset;
    while ( k < sqzsize) {
        n++;
        seqlen = (size_t *)(sqz->seqbuffer + k);
        k += lenbytes;
        seq =  sqz->seqbuffer + k;
        qual = sqz->qualbuffer + k;
        //Determine how much sequence is contained in the buffer
        /*This happens because not all the sequence can be stored in buffer*/
        seqread = *seqlen < (sqzsize - k)?*seqlen + 1:sqzsize - k;
        k += seqread;
        //encode sequence
        sqz_seqencode(seq, seqread, blk);
        //do something with the qualities
        sqz_qualencode(qual, seqread, blk);
    }
    //Indicate if there is more sequence to read
    if (seqread < *seqlen) {
        blk->newblk = 0;
        //Indicate how much sequence has been read
        sqz->toread = seqread;
        //Indicate length of sequence still needing loading
        sqz->prevlen = *seqlen;
    }
    if (k != sqzsize) {
        fprintf(stderr, "[sqzlib ERROR]: Failed to unloading raw data buffer\n");
        return 0;
    }
    return 1;
}


char sqz_tailblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    //Check if reading leftover sequence
    size_t *seqlen = &(sqz->prevlen);
    uint8_t *seq =  sqz->seqbuffer;
    uint8_t *qual = sqz->qualbuffer;
    //encode sequence
    sqz_seqencode(seq, sqz->offset, blk);
    //encode quality
    sqz_qualencode(qual, sqz->offset, blk);
    //Update how much sequence has been read
    sqz->toread += sqz->offset;
    blk->newblk = 1;
    //Indicate if there is more sequence to read
    if (sqz->toread < *seqlen) {
        fprintf(stderr, "More seq please!! %lu %lu\n",
                sqz->toread, sqz->prevlen);
        blk->newblk = 0;
    }
    return 1;
}


void sqz_seqencode(const uint8_t *seq, size_t seqlen, sqzblock_t *blk)
{
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
    //Write sequence length
    memcpy(codebuff + wbytes, &seqlen, sizeof(size_t));
    wbytes += sizeof(size_t);
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
            }
            //Loop over Ns in chunks of 127 and encode
            for (size_t i = 0; i < nn; i +=127) {
                //Encode 127 Ns or what is left
                wn = ((i+127) <= nn)?127:nn-i;
                //Set first bit if last N block before next base
                if (i+127 >= nn) wn = wn | NEND;
                memcpy(codebuff + wbytes, &wn, 1);
                wbytes++;
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
            code = sqz_bit2encode(lstop + i, ((i+32) <= blen)?32:blen-i);
            memcpy(codebuff + wbytes, &code, sizeof(uint64_t));
            wbytes += sizeof(uint64_t);
        }
    }
    //Move offset by number of bytes written
    blk->blksize += wbytes;
}


size_t sqz_seqdecode(const uint8_t *buff)
{
    size_t cbytes = 0;
    size_t seqlen;
    char *seqstr;
    char *qualstr;
    const unsigned char *nstr;
    size_t seqpos;
    uint64_t *mer;
    size_t blocklen;
    char merlen;
    size_t blockidx;
    size_t decoded;
    //slen
    seqlen = *(size_t *)buff;
    fprintf(stderr, "\t\t%lu\n", seqlen);
    cbytes += sizeof(size_t);
    seqstr = malloc(seqlen);
    qualstr = malloc(seqlen);
    if (!seqstr | !qualstr) {
        fprintf(stderr, "[sqzlib ERROR]: Insufficient memory.\n");
        free(qualstr);
        free(seqstr);
        return 0;
    }
    seqstr[seqlen - 1] = '\0';
    qualstr[seqlen - 1] = '\0';
    seqpos = 0;
    fprintf(stderr, ">>>>%p\n", seqstr) ;
    while (seqlen > 0) {
        //blen - length of block (within sequence up to first non acgt base)
        blocklen = *(size_t *)(buff + cbytes);
        fprintf(stderr, "^^block length %lu\n", blocklen);
        cbytes += sizeof(size_t);
        //TODO try to remove this cast
        mer = (uint64_t *)(buff + cbytes);
        blockidx = 0;
        for (size_t i = 0; i < blocklen; i+=32) {
            //Blocks code 32mers except last block which may code a shorter kmer
            merlen = ((i+32) <= blocklen)?32:blocklen - i - 1;
            decoded += sqz_bit2decode(mer + blockidx, seqstr + seqpos, merlen);
            //Keep traked of amount of decoded data: 64bit per mer
            cbytes += sizeof(uint64_t);
            //Move sequencepointerbyabount of bases decoded (32 except last mer)
            seqpos += merlen;
            //Track next block to decompress
            blockidx++;
        }
        seqlen -= blocklen;
        if (seqlen) {
            fprintf(stderr, "There are NS!!!!!!!\n");
            nstr = buff + cbytes;
            cbytes++;
            while (1) {
                unsigned char numn = *nstr & ~(1<<7);
                decoded += sqz_writens(*nstr & ~(1<<7), seqstr);
                seqpos += numn;
                seqlen -= numn;
                if (*nstr & 128)
                    break;
                nstr++;
                cbytes++;
            }
        }
    }

    fprintf(stderr, "bytes decoded: %lu\n", cbytes);
    sqz_qualdecode(buff + cbytes, qualstr, seqpos - 1);
    fprintf(stderr, "<<<<%p\n", seqstr);
    fprintf(stdout, "||%s\n", seqstr);
    free(seqstr);
    free(qualstr);
    return 0;
}

size_t sqz_qualencode(const uint8_t *qual, size_t quallen, sqzblock_t *blk)
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


uint8_t sqz_8binqual(uint8_t q)
{
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
    do {
    }
    while(seq_nt4_table[*(++seq)] < 4);
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


size_t sqz_qualdecode(const uint8_t *buff, char *uncode, size_t length)
{
    size_t byte = 0;
    size_t offset = 0;
    size_t decoded = 0;
    uint8_t code;
    unsigned char count;
    unsigned char q;
    while (decoded != length) {
        code = *(buff + byte); //get byte value
        count = code & 31;         //get current symbol count
        q = (code & 224) >> 5;     //get symbol index value
        for (int i = 0; i <= count; i++) {
            uncode[offset] = qual_val_table[q];
            offset++;
        }
        decoded += ++count;
        byte++;
    }
    uncode[offset] = '\0';
    return decoded;
}


unsigned char sqz_bit2decode(const uint64_t *mer, char *decoded, uint32_t len)
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
