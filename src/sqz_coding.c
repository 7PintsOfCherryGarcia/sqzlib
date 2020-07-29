#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include "sqz_coding.h"


int sqz_encode(sqzfastx_t *sqz, size_t sqzsize, sqzcodeblock_t *blk)
{
    unsigned char *seq = NULL;
    unsigned char *qual = NULL;
    size_t *seqlen = NULL;
    size_t seqread;
    size_t lenbytes = sizeof(size_t);
    size_t k = 0;
    size_t n = 0;
    uint8_t *qualbuff = NULL;
    //Check if reading leftover sequence
    if (sqz->endflag) {
        seqlen = &(sqz->prevlen);
        seq = (unsigned char *)sqz->seqbuffer + k;
        qual = (unsigned char *)sqz->qualbuffer + k;
        //encode sequence
        sqz_seqencode(seq, sqzsize, blk);
        //encode quality
        sqz_qualencode(qual, qualbuff);
        //If this happens the buffer will only contain data from the sequence being read
        k = sqzsize;
        //Update how much sequence has been read
        sqz->toread += sqzsize;
        seqread = sqz->toread;
        sqz->endflag = 0;
    }
    //Process sequence data in buffer
    while ( k < sqzsize) {
        n++;
        //Extract length of sequences in buffer
        seqlen = (size_t *)(sqz->seqbuffer + k);
        k += lenbytes;
        seq = (unsigned char *)sqz->seqbuffer + k;
        qual = (unsigned char *)sqz->qualbuffer + k;
        //Determine how much sequence is contained in the buffer
        seqread = *seqlen < (sqzsize - k)?*seqlen + 1:sqzsize - k;
        //fprintf(stderr, "In sequence: %lu read %lu len %lu\n", n, seqread, *seqlen);
        k += seqread;
        //encode sequence
        sqz_seqencode(seq, seqread, blk);
        //do something with the qualities
        sqz_qualencode(qual, qualbuff);
    }
    //Indicate if there is more sequence to read
    if (seqread < *seqlen) {
        //fprintf(stderr, "More seq please!! %lu %lu\n", seqread, *seqlen);
        sqz->endflag = 1;
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


size_t sqz_seqencode(const unsigned char *str, uint32_t strlen, sqzcodeblock_t *codeblk)
{
    //TODO control for cmpbuff length
    //Set coding buffer to the next unused position in the block
    uint8_t *cmpbuff = codeblk->codebuff + codeblk->offset;
    const unsigned char *lstop = str;         //Track position within sequence
	  const unsigned char *nptr = str + strlen; //End of string
	  uint32_t nn;                              //Number of Ns
    unsigned char wn;                         //127 N block. 1 bit flag 7 bit count
	  const unsigned char *npos;                //Track positions where N occurs
	  uint32_t blocklen = 0;                    //Length of segment before first N
    uint64_t code;                            //2 bit encoded sequence
    size_t wbytes = 0;                        //Number of bytes written
    //Write sequence length
    memcpy(cmpbuff + wbytes, &strlen, sizeof(uint32_t));
    wbytes += sizeof(uint32_t);
    //Search for stretches of Ns
    do {
        //Get position of first non ACGTacgt base
	      npos = sqz_findn(lstop);
	      if (*npos) {
            //Determine block length up to found N
            blocklen = npos - lstop;
            //Determine number of consecutive Ns until next base
            nn = 0;
	          while ( seq_nt4_table[*npos] == 4) {
	              nn++;
		            npos++;
	          }
            //Write block length [blen]
            memcpy(cmpbuff + wbytes, &blocklen, sizeof(uint32_t));
            wbytes += sizeof(uint32_t);
            //Loop over sequence block in chunks of 32 bases
            for (uint32_t i = 0; i < blocklen; i+=32) {
                //Encode 32mer or what is left of sequence block
                code = bit2encode(lstop + i, ((i+32) <= blocklen)?32:blocklen-i);
                memcpy(cmpbuff + wbytes, &code, sizeof(uint64_t));
                wbytes += sizeof(uint64_t);
            }
            //Loop over Ns in chunks of 127
            for (uint32_t i = 0; i < nn; i +=127) {
                //Encode 127 Ns or what is left
                wn = ((i+127) <= nn)?127:nn-i;
                //Set first bit if last N block before next base
                if (i+127 >= nn) wn = wn | NEND;
                memcpy(cmpbuff + wbytes, &wn, 1);
                wbytes++;
            }
            //Trace next base after sequence of Ns
            lstop = npos;
	      }
    } while (*npos);
    //Detect and encode trailing bases
    blocklen = nptr - lstop;
    if (blocklen) {
        memcpy(cmpbuff + wbytes, &blocklen, sizeof(uint32_t));
        wbytes += sizeof(uint32_t);
        for (uint32_t i = 0; i < blocklen; i+=32) {
            code = bit2encode(lstop + i, ((i+32) <= blocklen)?32:blocklen-i);
            memcpy(cmpbuff + wbytes, &code, sizeof(uint64_t));
            wbytes += sizeof(uint64_t);
        }
    }
    //Move offset by number of bytes written
    codeblk->offset += wbytes;
    return wbytes;
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


const unsigned char *sqz_findn(const unsigned char *strseq)
{
    do {
    }
    while(seq_nt4_table[*(++strseq)] < 4);
    return strseq;
}


uint64_t bit2encode(const unsigned char *str, uint32_t strlen)
{
    uint64_t result = 0;
    for (uint32_t i = 0; i < strlen; i++) {
        result = (result << 2) | seq_nt4_table[str[i]];
    }
    return result;
}
