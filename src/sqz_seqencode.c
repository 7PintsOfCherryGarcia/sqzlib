#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "sqz_data.h"
#include "sqz_seqencode.h"


unsigned char bit2decode(const uint64_t *mer, char *decoded, uint32_t len)
{
    unsigned char byte;
    unsigned char nbase = 0;
    uint64_t code = *mer;
    --len;
    do {
        byte = code & TWO_BIT_MASK;
        decoded[len] = seq_dec_table[byte];
        nbase += 1;
        code >>= 2;
    } while (len-- != 0);
    return nbase;
}


uint64_t bit2encode(const unsigned char *str, uint32_t strlen)
{
    uint64_t result = 0;
    for (uint32_t i = 0; i < strlen; i++) {
        result = (result << 2) | seq_nt4_table[str[i]];
    }
    return result;
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


const unsigned char *sqz_findn(const unsigned char *strseq)
{
    do {
    }
    while(seq_nt4_table[*(++strseq)] < 4);
    return strseq;
}


size_t sqz_seqdecode(const uint8_t *buff)
{
    uint32_t strlen;
    char *seqpos;
    uint32_t blocklen;
    const uint64_t *mer;
    unsigned char merlen;
    uint32_t blockidx;
    const unsigned char *nstr;
    size_t bytes = 0;
    uint32_t decoded = 0;
    char *seqstr;
    //Get sequence length of current block
    memcpy(&strlen, buff, sizeof(uint32_t));
    bytes += sizeof(uint32_t);
    //Allocate memory for current sequence
    seqstr = malloc(strlen + 1);
    if (!seqstr) {
        fprintf(stderr, "MEM ERROR\n");
    }
    uint32_t strlen2 = strlen;
    seqpos = seqstr;
    //Decode for as long as there are bases to decode
    while (strlen > 0) {
        //Get block length
        memcpy(&blocklen, buff + bytes, sizeof(uint32_t));
        bytes += sizeof(uint32_t);
        //Get starting porition of 64bit integers
        mer = (uint64_t*)(buff + bytes);
        //Track block position
        blockidx = 0;
        //Loop over blocks
        //fprintf(stderr, "DECODING block of length: %u\n", blocklen);
        for (uint32_t i = 0; i < blocklen; i+=32) {
            //Blocks code 32mers except last block which may code a shorter kmer
            merlen = ((i+32) <= blocklen)?32:blocklen - i;
            decoded += bit2decode(mer + blockidx, seqpos, merlen);
            //Keep traked of amount of decoded data: 64bit per mer
            bytes += sizeof(uint64_t);
            //Move sequence pointer by abount of bases decoded (32 except last mer)
            seqpos += merlen;
            //Track next block to decompress
            blockidx++;
        }
        //TODO make sure decoded bases equals block length
        //Update amount of data left to decompress
        strlen -= blocklen;
        //strlen is != 0 when blocks are followed by an Nblock
        if (strlen) {
            //Read nbyte
            nstr = buff + bytes;
            bytes++;
            while (1) {
                unsigned char numn = *nstr & ~(1<<7);
                //fprintf(stderr, "N value of: %u\n", numn);
                decoded += sqz_writens(*nstr & ~(1<<7), seqpos);
                seqpos += numn;
                strlen-= numn;
                if (*nstr & 128)
                    break;
                nstr++;
                bytes++;
            }
        }
    }

    FILE *raw = fopen("decomseq", "wb");
    fwrite(seqstr, 1, strlen2+1, raw);
    fclose(raw);


    free(seqstr);
    return strlen;
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
