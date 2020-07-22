#include <string.h>
#include "sqz_qualencode.h"


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


size_t sqz_qualencode(const char *strqual, uint8_t *codebuff)
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


size_t sqz_qualdecode(uint8_t *codebuff, char *uncode, size_t length)
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
