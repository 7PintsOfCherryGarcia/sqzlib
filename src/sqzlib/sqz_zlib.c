#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include "sqz_data.h"

size_t sqz_deflate(sqzblock_t *blk, int level)
{
    int ret;
    size_t wbytes = 0;
    size_t have;
    z_stream strm;
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
    strm.avail_in = blk->blkpos;
    strm.next_in  = blk->blkbuff;
    do {
        strm.avail_out = blk->cmpsize;
        strm.next_out = blk->cmpbuff;
        ret = deflate(&strm, Z_FINISH);    /* no bad return value */
        if ( ret == Z_STREAM_ERROR ) {
            fprintf(stderr, "[sqzlib ZLIB ERROR]: Deflate error.\n");
            return (size_t)ret;  /* state not clobbered */
        }
        //check there is enough space in buffer
        have = blk->cmpsize - strm.avail_out;
        if ( (wbytes + have) > blk->cmpsize) {
            fprintf(stderr, "[sqzlib ZLIB ERROR]: Compression error.\n");
            ret = Z_ERRNO;
            break;
        }
        wbytes += have;
    } while (strm.avail_out == 0);
    //Check compression went well
    if (ret != Z_STREAM_END) wbytes = 0;
    deflateEnd(&strm);
    return wbytes;
}


uint64_t sqz_inflate(sqzblock_t *blk)
{
    //TODO test with internal memory area to put decompressed data + memcpy
    int ret;
    size_t wbytes = 0;
    size_t have;
    z_stream strm;
    // allocate deflate state
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK) {
        fprintf(stderr, "[sqzlib ZLIB ERROR]: Failed to init stream.\n");
        return (size_t)ret;
    }
    //Main compression loop
    strm.avail_in = blk->cmpsize;
    strm.next_in = blk->cmpbuff;
    do {
        strm.avail_out = blk->blksize;
        strm.next_out  = blk->blkbuff;
        ret = inflate(&strm, Z_SYNC_FLUSH);    /* no bad return value */
        switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;
                break;
            case Z_DATA_ERROR:
                break;
            case Z_MEM_ERROR:
                break;
            case Z_STREAM_ERROR:
                fprintf(stderr, "[sqzlib ZLIB ERROR]: inflate error.\n");
                return (size_t)ret;  /* state not clobbered */
        }
        //check there is enough space in buffer
        have = blk->blksize - strm.avail_out;
        //fprintf(stderr, "**\t\tdcp: %lu\n", blk->blksize);
        if ( (wbytes + have) > blk->blksize) {
            fprintf(stderr,
                    "[sqzlib ZLIB ERROR]: Decompression error %lu.\n", have);
            ret = Z_ERRNO;
            break;
        }
        //Write to output buffer, keep track of written bytes
        wbytes += have;
    } while (strm.avail_out == 0);
    //Check compression went well
    if (ret != Z_STREAM_END) wbytes = 0;
    inflateEnd(&strm);
    return wbytes;
}
