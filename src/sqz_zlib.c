#include "sqz_zlib.h"

size_t sqz_deflate(sqzblock_t *blk, int level)
//size_t sqz_deflate(void *seq,
//                   size_t seqlength,
//                   uint8_t *dest,
//                   size_t destlen,
//                   int level)
{
    int ret;
    size_t wbytes = 0;
    size_t have;
    z_stream strm;
    unsigned char out[CHUNK];
    // allocate deflate state
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = deflateInit(&strm, level);
    if (ret != Z_OK) {
        fprintf(stderr, "[sqzlib ZLIB ERROR]: Failed to init stream.\n");
        return (size_t)ret;
    }
    //Main compression loop
    strm.avail_in = blk->blksize;
    strm.next_in = blk->codebuff;
    do {
        strm.avail_out = CHUNK;
        strm.next_out = out;
        ret = deflate(&strm, Z_FINISH);    /* no bad return value */
        if ( ret == Z_STREAM_ERROR ) {
            fprintf(stderr, "Some error\n");
            return (size_t)ret;  /* state not clobbered */
        }
        //check there is enough space in buffer
        have = CHUNK - strm.avail_out;
        if ( (wbytes + have) > blk->cmpsize) {
            fprintf(stderr,
                    "Compression error %lu %lu %lu\n",
                    wbytes, have, blk->cmpsize);
            ret = Z_ERRNO;
            break;
        }
        //Write to output buffer, keep track of written bytes
        memcpy(blk->cmpbuff + wbytes, out, have);
        wbytes += have;
    } while (strm.avail_out == 0);
    //Check compression went well
    if (ret != Z_STREAM_END) wbytes = 0;
    deflateEnd(&strm);
    return wbytes;
}
