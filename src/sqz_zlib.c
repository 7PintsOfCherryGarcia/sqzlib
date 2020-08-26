#include "sqz_zlib.h"


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
    strm.avail_in = blk->blksize;
    strm.next_in = blk->codebuff;
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


char sqz_zlibcmpdump(sqzblock_t *blk, size_t size, FILE *ofp)
{
    size_t bytes = sizeof(size_t);
    size_t wbytes = 0;
    fprintf(stderr, "[sqzlib INFO]: Dumping block.\n");
    //Write uncompressed number of bytes in block
    wbytes += fwrite(&(blk->blksize), bytes, 1, ofp);
    //Write compressed number of bytes in block
    wbytes += fwrite(&(size), bytes, 1, ofp);
    //Write block
    wbytes += fwrite(blk->cmpbuff, 1, size, ofp);
    if (wbytes != size+2) return 0;
    return 1;
}



size_t sqz_inflate(sqzblock_t *blk)
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
        strm.next_out = blk->codebuff;
        ret = inflate(&strm, Z_SYNC_FLUSH);    /* no bad return value */
        switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
            case Z_STREAM_ERROR:
                fprintf(stderr, "[sqzlib ZLIB ERROR]: inflate error.\n");
                return (size_t)ret;  /* state not clobbered */
        }
        //check there is enough space in buffer
        have = blk->blksize - strm.avail_out;
        if ( (wbytes + have) > blk->blksize) {
            fprintf(stderr,
                    "[sqzlib ZLIB ERROR]: Decompression error %lu.\n", have);
            ret = Z_ERRNO;
            break;
        }
        //Write to output buffer, keep track of written bytes
        //memcpy(blk->codebuff + wbytes, out, have);
        wbytes += have;
    } while (strm.avail_out == 0);
    //Check compression went well
    if (ret != Z_STREAM_END) wbytes = 0;
    inflateEnd(&strm);
    return wbytes;
}
