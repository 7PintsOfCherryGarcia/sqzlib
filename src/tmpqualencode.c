#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)


unsigned char qual_val_table[8] = {33,39,48,55,60,66,70,73};

char bin8qual(char q) {
    if (q >= 73) return 7<<5;
    if ( (q < 73) && (q >= 68) ) return 6<<5;
    if ( (q < 68) && (q >= 63) ) return 5<<5;
    if ( (q < 63) && (q >= 58) ) return 4<<5;
    if ( (q < 58) && (q >= 53) ) return 3<<5;
    if ( (q < 53) && (q >= 43) ) return 2<<5;
    if ( (q < 43) && (q >= 35) ) return 1<<5;
    return 0;
}


size_t qdecode0bitlength(uint8_t *codebuff, char *uncode, size_t length)
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


size_t qencode0bitlength(const char *strqual, uint8_t *codebuff)
{
    size_t bytes = 0;
    char q = bin8qual(*strqual);
    unsigned char  c = 0;
    unsigned char code = 0;
    while(*strqual) {
        if ( bin8qual(*(strqual+1)) == q) {
            c++;
            if (c == 31) {
                //Encode
                code = code | q;
                code = code | c;
                memcpy(codebuff + bytes, &code, 1);
                bytes += 1;
                c = 0;
                code = 0;
                q = bin8qual(*(++strqual));
            }
            strqual++;
            continue;
        }
        //Encode
        code = code | q;
        code = code | c;
        memcpy(codebuff + bytes, &code, 1);
        strqual++;
        q = bin8qual(*strqual);
        c = 0;
        bytes += 1;
        code = 0;
    }
    return bytes;

}


int main(int argc, char *argv[]) {
    gzFile fp = gzopen(argv[1], "r");
    if (!fp) return -1;
    kseq_t *seq = kseq_init(fp);
    if (!seq) return -1;
    int l;
    uint8_t *code = malloc(16*1024*1024);
    char *uncode = malloc(16*1024*1024);
    size_t bytes;
    size_t tcbytes = 0;
    size_t tqbytes = 0;
    while ((l = kseq_read(seq)) >= 0) {
        bytes = qencode0bitlength(seq->qual.s, code);
        tcbytes += bytes;
        tqbytes += seq->qual.l;
        qdecode0bitlength(code, uncode, seq->qual.l);
        if (strlen(uncode) != seq->qual.l)
            fprintf(stderr, "ERROR: decode: %lu og %lu\n", strlen(uncode), seq->qual.l);
    }
    fprintf(stderr,
            "total q bytes: %lu, total c bytes: %lu, ratio: %f\n",
            tqbytes, tcbytes, (float)tcbytes/(float)tqbytes);
    kseq_destroy(seq);
    gzclose(fp);
    free(code);
    free(uncode);
}

/*
char bin2qual(char q) {
  if (q >= 53) return 63;
  return  10;
}

size_t qdecode1bitlength(uint8_t *codebuff, char *uncode)
{
    size_t byte = 0;
    size_t offset = 0;
    uint8_t code;
    unsigned char last;
    unsigned char count;
    unsigned char q;
    int len = 0;
    while (1) {
        code = *(codebuff + byte); //get byte value
        last = code & (1<<7);      //get 7th bit state
        count = code & 15;         //get current symbol count
        q = (code & 112) >> 4;     //get symbol index value
        for (int i = 0; i <= count; i++) {
            uncode[offset] = qual_val_table[q];
            offset++;
        }
        len += ++count;
        byte++;
        if (last) break;
    }
    uncode[offset] = '\0';
    return len;
}


size_t qencode1bitlength(const char *strqual, uint8_t *codebuff)
{
    size_t bytes = 0;
    unsigned char last = 0;
    char q = bin8qual1(*strqual);
    unsigned char  c = 0;
    unsigned char code = 0;
    while(*strqual) {
        if ( bin8qual1(*(strqual+1)) == q) {
            c++;
            if (c == 15) {
                //Encode
                if (!(*(strqual+2))) last = 128;
                code = code | last;
                code = code | q;
                code = code | c;
                memcpy(codebuff + bytes, &code, 1);
                bytes += 1;
                c = 0;
                code = 0;
                q = bin8qual1(*(++strqual));
            }
            strqual++;
            continue;
        }
        //Encode
        if (!(*(strqual+1))) last = 128;
        code = code | last;
        code = code | q;
        code = code | c;
        memcpy(codebuff + bytes, &code, 1);
        strqual++;
        q = bin8qual1(*strqual);
        c = 0;
        bytes += 1;
        code = 0;
    }
    return bytes;
}


char bin8qual1(char q) {
    if (q >= 73) return 7<<4;
    if ( (q < 73) && (q >= 68) ) return 6<<4;
    if ( (q < 68) && (q >= 63) ) return 5<<4;
    if ( (q < 63) && (q >= 58) ) return 4<<4;
    if ( (q < 58) && (q >= 53) ) return 3<<4;
    if ( (q < 53) && (q >= 43) ) return 2<<4;
    if ( (q < 43) && (q >= 35) ) return 1<<4;
    return 0;
}

*/
