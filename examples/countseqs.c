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
    size_t n = 0;
    sqzfp = sqzopen(argv[1], "r");
    if ( !sqzfp ) return -1;

    seq = kseq_init(sqzfp);
    while ( (l = kseq_read(seq)) >= 0)
        ++n;
    if (l != -1) {
        fprintf(stderr, "ERROR: malformated FASTX\n");
        goto exit;
    }
    fprintf(stderr, "%lu sequences\n", n);
    ret = 0;
    exit:
        kseq_destroy(seq);
        sqzclose(sqzfp);
    return ret;
}
