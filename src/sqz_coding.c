#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include "sqz_coding.h"


char sqz_encode(sqzfastx_t *sqz, sqzblock_t *blk)
{
    fprintf(stderr, "|||%u\n", LOAD_SIZE);
    if (blk->newblk) {
        return sqz_headblk(sqz, blk);
    }
    else {
         return sqz_tailblk(sqz, blk);
    }
}


char sqz_headblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    uint8_t *seq = NULL;
    //char *qual = NULL;
    size_t *seqlen = NULL;
    size_t seqread;
    size_t lenbytes = sizeof(size_t);
    size_t k = 0;
    size_t n = 0;
    size_t sqzsize = sqz->offset;
    //uint8_t *qualbuff = NULL;
    //Process sequence data in buffer
    while ( k < sqzsize) {
        n++;
        //Extract length of sequences in buffer
        seqlen = (size_t *)(sqz->seqbuffer + k);
        k += lenbytes;
        seq = sqz->seqbuffer + k;
        //qual = (char *)sqz->qualbuffer + k;
        //Determine how much sequence is contained in the buffer
        seqread = *seqlen < (sqzsize - k)?*seqlen + 1:sqzsize - k;
        k += seqread;
        //encode sequence
        sqz_seqencode(seq, seqread, blk);
        //do something with the qualities
        //sqz_qualencode(qual, qualbuff);
    }
    //Indicate if there is more sequence to read
    if (seqread < *seqlen) {
        //fprintf(stderr, "More seq please!! %lu %lu\n", seqread, *seqlen);
        blk->newblk = 0;
        //Indicate how much sequence has been read
        sqz->toread = seqread;
        //Indicate length of sequence still needing loading
        sqz->prevlen = *seqlen;
    }
    if (k != sqzsize) {
        fprintf(stderr, "[ERROR] Unloading buffer\n");
        return 0;
    }
    return 1;
}


char sqz_tailblk(sqzfastx_t *sqz, sqzblock_t *blk)
{
    fprintf(stderr, "|| from seqlen %lu\n", sqz->prevlen);
    //Check if reading leftover sequence
    size_t *seqlen = &(sqz->prevlen);
    uint8_t *seq = sqz->seqbuffer;
    //char *qual = (char *)sqz->qualbuffer;
    //encode sequence
    sqz_seqencode(seq, sqz->offset, blk);
    //encode quality
    //sqz_qualencode(qual, qualbuff);
    //Update how much sequence has been read
    sqz->toread += sqz->offset;
    blk->newblk = 1;
    //Indicate if there is more sequence to read
    if (sqz->toread < *seqlen) {
        fprintf(stderr, "More seq please!! %lu %lu\n",
                sqz->toread, sqz->prevlen);
        blk->newblk = 0;
        //Indicate length of sequence still needing loading
        //sqz->prevlen = *seqlen;
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
	  size_t nn;                     //Number of Ns
    unsigned char wn;                //127 N block. 1 bit flag 7 bit count
	  const uint8_t *npos;       //Track positions where N occurs
	  size_t blen = 0;                 //Length of segment before first N
    uint64_t code;                   //2 bit encoded sequence
    size_t wbytes = 0;               //Number of bytes used to encode
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


size_t sqz_qualencode(const unsigned char *strqual, uint8_t *codebuff)
{
    size_t bytes = 0;
    unsigned char q = sqz_8binqual(*strqual);
    unsigned char  c = 0;
    unsigned char code = 0;
    while(*strqual) {
        if ( sqz_8binqual(*(strqual+1)) == q) {
            c++;
            if (c == 31) {
                //Encode
                code = code | q;
                code = code | c;
                memcpy(codebuff + bytes, &code, 1);
                bytes += 1;
                c = 0;
                code = 0;
                q = sqz_8binqual(*(++strqual));
            }
            strqual++;
            continue;
        }
        //Encode
        code = code | q;
        code = code | c;
        memcpy(codebuff + bytes, &code, 1);
        strqual++;
        q = sqz_8binqual(*strqual);
        c = 0;
        bytes += 1;
        code = 0;
    }
    return bytes;
}


unsigned char sqz_8binqual(char q)
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
    blk->codebuff = malloc(size);
    if (!blk->codebuff) {
        fprintf(stderr, "[libsqz ERROR]: memory error.\n");
        free(blk);
        return NULL;
    }
    blk->blksize = 0;
    blk->newblk = 1;
    return blk;
}


void sqz_blkdestroy(sqzblock_t *blk)
{
    free(blk->codebuff);
    free(blk);
}
