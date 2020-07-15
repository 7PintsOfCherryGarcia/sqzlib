#include <stdio.h>

int sqz_main(const char *filename);

int main(int argc, char *argv[]) {
    int ret = -1;
    if (argc < 2) goto exit;

    if (!sqz_main(argv[1])) goto exit;

    ret = 0;
    exit:
        return ret;
}
