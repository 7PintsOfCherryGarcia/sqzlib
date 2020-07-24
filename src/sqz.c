#include <stdio.h>
#include "sqz.h"
#include "squeezmalib.h"



int main(int argc, char *argv[]) {
    int ret = -1;
    if (argc < 2) goto exit;

    if (!sqz_compress(argv[1])) goto exit;

    ret = 0;
    exit:
        return ret;
}


int sqz_compress(const char *filename)
{
    int ret = 0;
    //Initialize data main sqz data structure
    sqzfastx_t *sqz = sqz_fastxinit(filename, LOAD_SIZE);
    if (!sqz) {
        fprintf(stderr, "[sqz ERROR: INIT] Failed to start data structures.\n");
        goto exit;
    }
    //Check for format
    switch (sqz->fmt) {
    case 1:
        //if (!sqz_fasta(sqz, seq)) goto exit;
        break;
    case 2:
        if (!sqz_squeezefastq(sqz)) goto exit;
        break;
    }

    ret = 1;
    exit:
        sqz_kill(sqz);
        return ret;
}


