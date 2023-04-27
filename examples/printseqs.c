#include "sqzlib.h"
#include "klib/kseq.h"
KSEQ_INIT(sqzFile, sqzread)


int main(int argc, char *argv[])
{
    int ret = 1;
    if (argc == 1) {
        printf("Usage: countseqs <file>\n");
        return 0;
    }

    sqzFile sqzfp;
    kseq_t *seq;
    int l;
    sqzfp = sqzopen(argv[1], "r");
    if ( !sqzfp ) return -1;

    seq = kseq_init(sqzfp);
    while ( (l = kseq_read(seq)) >= 0)
        fprintf(stdout, ">%s\n%s\n", seq->name.s, seq->seq.s);
    if (l != -1) {
        fprintf(stderr, "ERROR: malformated FASTX %d\n", l);
        goto exit;
    }
    ret = 0;
    exit:
        kseq_destroy(seq);
        sqzclose(sqzfp);
    return ret;
}
