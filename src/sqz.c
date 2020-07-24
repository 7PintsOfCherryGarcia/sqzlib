#include <stdio.h>
#include "sqz.h"
#include "squeezmalib.h"



int main(int argc, char *argv[]) {
    int ret = -1;
    if (argc < 2) goto exit;

    if (!sqz_main(argv[1])) goto exit;

    ret = 0;
    exit:
        return ret;
}


int sqz_main(const char *filename)
{
    int ret = 0;
    //Initialize data structures
    //gzFile fp = gzopen(filename, "r");
    //if (!fp) return ret;
    //kseq_t *seq = kseq_init(fp);
    //if (!seq) {
    //    gzclose(fp);
    //    return ret;
    //}
    sqzfastx_t *sqz = sqz_fastxinit(filename, LOAD_SIZE);
    if (!sqz) {
        //kseq_destroy(seq);
        //gzclose(fp);
        sqz = NULL;
        goto exit;
    }
    //Check for format
    //switch (sqz->fmt) {
    //case 1:
    //    if (!sqz_loadfasta(sqz, seq)) goto exit;
    //    break;
    //case 2:
    //    if (!sqz_loadfastq(sqz, seq)) goto exit;
    //    break;
    //}

    ret = 1;
    exit:
        sqz_kill(sqz);
        //kseq_destroy(seq);
        //gzclose(fp);
        return ret;
}
