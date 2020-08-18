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
        //Loop over blocks and decode
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
