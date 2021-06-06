#include <stdio.h>

#define SQZLIB
#define KLIB
#include "sqz_data.h"

char zbytes[4] = {0, 0, 0, 0};
char magic[4] = {5, 8, 5, 9};
unsigned char cmpflag = 1;

char sqz_filetail(size_t numseqs, FILE *ofp)
{
    if ( 4 != fwrite(zbytes, 1, 4, ofp) ) {
        return 0;
    }
    if ( 1 != fwrite(&numseqs, sizeof(numseqs), 1, ofp) ) {
        return 0;
    }
    if ( 4 != fwrite(zbytes, 1, 4, ofp) ) {
        return 0;
    }
    return 1;
}


char sqz_filehead(unsigned char fmt, FILE *ofp)
{
    char wbytes = 0;
    if ( 4 != (wbytes += fwrite(magic, 1, 4, ofp)) ) return 0;
    if ( 5 != (wbytes += fwrite(&fmt,  1, 1, ofp)) ) return 0;
    //Compression library
    if ( 6 != (wbytes += fwrite(&cmpflag, 1, 1, ofp)) ) return 0;
    if ( 8 != (wbytes += fwrite(zbytes,   1, 2, ofp)) ) return 0;
    return wbytes;
}
