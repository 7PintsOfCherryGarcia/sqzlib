#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "sqz_data.h"

const char zbytes[4] = {0, 0, 0, 0};
const uint8_t magic1[4] = {5, 8, 5, 9};
const uint8_t magic2[4] = {9, 5, 8, 5};

unsigned char cmpflag = 1;


static void sqz_fastxreset(sqzfastx_t *sqz)
{
    sqz->endflag    = 0;
    sqz->cmpflag    = 0;
    sqz->offset     = 0;
    sqz->namepos    = 0;
    sqz->n          = 0;
    sqz->bases      = 0;
    sqz->rem        = 0;
    sqz->prevlen    = 0;
}

static void sqz_blkreset(sqzblock_t *blk)
{
    blk->blkpos  = 0;
    blk->namepos = 0;
    blk->newblk  = 1;
    blk->cmppos  = 0;
}

static void sqz_decode(sqzfastx_t *sqz,
                       sqzblock_t *blk,
                       uint8_t fmt,
                       uint64_t klibl)
{
    switch (fmt) {
    case 2:
        sqz->offset = sqz_fastXdecode(blk, sqz->readbuffer, klibl, 1);
        break;
    case 1:
        sqz->offset = sqz_fastXdecode(blk, sqz->readbuffer, klibl, 0);
        break;
    }
}

int64_t sqz_blkcompress(sqzblock_t *blk, int level, uint8_t libfmt)
{
    switch (libfmt) {
        case 1:
            return sqz_deflate(blk, level);
        case 2:
            return sqz_zstdcompress(blk, level);
    }
    return -1;
}

uint64_t sqzdecompress(sqzblock_t *blk, uint8_t libfmt)
{
    switch (libfmt) {
        case 1:
            return sqz_inflate(blk);
        case 2:
            return sqz_zstddecompress(blk);
    }
    return 0;
}

void sqzrewind(sqzFile sqzfp)
{
    sqz_fastxreset(sqzfp->sqz);
    sqz_blkreset(sqzfp->blk);
    sqz_gzseek(sqzfp, HEADLEN, SEEK_SET);
    sqzfp->filepos = sqz_gztell(sqzfp);
    sqzfp->ff = 0;
}

char sqz_filehead(uint8_t fmt, uint8_t libfmt, FILE *ofp)
{
    char wbytes = 0;
    if ( 4 != (wbytes += fwrite(magic1,  1, 4, ofp)) ) return 0;
    if ( 5 != (wbytes += fwrite(&fmt,    1, 1, ofp)) ) return 0;
    if ( 6 != (wbytes += fwrite(&libfmt, 1, 1, ofp)) ) return 0;
    if ( 8 != (wbytes += fwrite(zbytes,  1, 2, ofp)) ) return 0;
    return wbytes;
}

char sqz_filetail(uint64_t numseqs, uint32_t nblocks, FILE *ofp)
{
    if ( 4 != fwrite(zbytes,     1, 4, ofp) ) return 0;
    if ( 1 != fwrite(&numseqs, B64, 1, ofp) ) return 0;
    if ( 2 != fwrite(zbytes,     1, 2, ofp) ) return 0;
    if ( 1 != fwrite(&nblocks, B32, 1, ofp) ) return 0;
    if ( 2 != fwrite(zbytes,     1, 2, ofp) ) return 0;
    if ( 4 != fwrite(magic2,     1, 4, ofp) ) return 0;
    return 1;
}

uint8_t sqz_blkdump(sqzblock_t *blk, uint64_t cmpsize, FILE *ofp)
{
    size_t wbytes = 0;
    //Write uncompressed number of bytes in block
    wbytes += fwrite(&(blk->blkpos), B64, 1, ofp);
    //Write compressed number of bytes in block
    wbytes += fwrite(&cmpsize, B64, 1, ofp);
    //Write block
    wbytes += fwrite(blk->cmpbuff, cmpsize, 1, ofp);
    if (wbytes != 3) return 0;
    return 1;
}

void sqz_resetblk(sqzblock_t *blk)
{
    blk->blkpos  = 0;
}

uint8_t sqz_newblk(sqzblock_t *blk)
{
    return blk->newblk;
}

uint64_t sqz_filesize(const char *filename)
{
    if (!filename) return 0;
    FILE *fp = fopen(filename, "r");
    if ( fseek(fp, 0, SEEK_END) ) return 0;
    uint64_t s = ftell(fp);
    fclose(fp);
    return s - 22;
}

/*
int sqz_fseek(sqzFile sqzfp, long offset, int whence)
{
    return fseek(sqzfp->fp, offset, whence);
}
uint64_t sqz_fread(void *ptr, uint64_t size, uint64_t nmemb, sqzFile sqzfp)
{
    return fread(ptr, size, nmemb, sqzfp->fp);
}
*/

uint64_t sqz_getfilepos(sqzFile sqzfp)
{
    return ( sqzfp->filepos = sqz_gztell(sqzfp) );
}

uint8_t sqz_isfq(sqzFile sqzfp)
{
    return (sqzfp->fmt & 3) == 2 ? 1 : 0;
}

uint8_t sqz_sqzgetcmplib(sqzFile sqzfp)
{
    return sqzfp->fmt >> 3;
}

sqzblock_t *sqz_sqzgetblk(sqzFile sqzfp)
{
    return sqzfp->blk;
}

uint8_t sqz_format(sqzFile sqzfp)
{
    return sqzfp->fmt;
}

void sqzclose(sqzFile sqzfp)
{
    sqz_gzclose(sqzfp);
    if (sqzfp->sqz) sqz_fastxkill(sqzfp->sqz);
    if (sqzfp->blk) sqz_blkkill(sqzfp->blk);
    free(sqzfp);
}

uint8_t sqz_readblksize(sqzblock_t *blk, sqzFile sqzfp, uint8_t libfmt)
{
    uint8_t  ret = 1;
    uint64_t cmpsize;
    uint64_t dcpsize;
    uint64_t nelem;
    nelem = sqz_gzread(sqzfp, &dcpsize, B64);
    nelem += sqz_gzread(sqzfp, &cmpsize, B64);
    if ( nelem != 16 ) goto exit;
    if ( cmpsize > (blk->mcmpsize) ) {
        blk->cmpbuff = realloc(blk->cmpbuff, cmpsize);
        if ( !(blk->cmpbuff) ) goto exit;
        blk->mcmpsize = cmpsize;
    }
    if ( dcpsize > (blk->mblksize) ) {
        blk->blkbuff = realloc(blk->blkbuff, dcpsize);
        if ( !(blk->blkbuff) ) goto exit;
        blk->mblksize = dcpsize;
    }
    blk->cmpsize = cmpsize;
    blk->blksize = dcpsize;
    if ( (cmpsize != (uint64_t)sqz_gzread(sqzfp, blk->cmpbuff, cmpsize) ) )
        goto exit;
    if (dcpsize != sqzdecompress(blk, libfmt))
        goto exit;
    blk->newblk = 1;
    ret = 0;
    exit:
        return ret;
}

int64_t sqzread(sqzFile file, void *buff, uint64_t len)
{
    if (!file || !buff) return -1;
    if ( !(file->fmt & 4) )
        return sqz_gzread(file, buff, (uint32_t)len);
    sqzfastx_t *sqz  = file->sqz;
    sqzblock_t *blk  = file->blk;
    uint8_t     fmt  = file->fmt & 3;
    uint8_t *outbuff = (uint8_t *)buff;
    int64_t read     = 0;
    switch (file->ff & 127) {
    case 0:
        {
        if (!sqz_readblksize(file->blk, file, file->libfmt)) goto error;
        if ( sqz_gztell(file) == file->size)
            file->ff |= 128U;
        sqz_decode(sqz, blk, fmt, LOAD_SIZE);
        read = sqz->offset > len ? len : sqz->offset;
        memcpy(outbuff, sqz->readbuffer, read);
        sqz->rem += read;
        sqz->offset -= read;
        if (sqz->offset < len)
            file->ff = (file->ff & 128) | 2U;
        else
            file->ff = (file->ff & 128) | 1U;
        return read;
        }
    case 1: //Data just need to be copied
        {
        memcpy(outbuff, sqz->readbuffer + sqz->rem, len);
        sqz->rem += len;
        sqz->offset -= len;
        if (sqz->offset < len)
            file->ff = (file->ff & 128U) | 2U;
        return len;
        }
    case 2: //Data can be copied but entire buffer can't be filled
        {
        ;
        uint64_t outpos;
        uint64_t leftover = sqz->offset;
        memcpy(outbuff, sqz->readbuffer + sqz->rem, sqz->offset);
        outpos = leftover;
        if (blk->newblk)
            sqz_decode(sqz, blk, fmt, LOAD_SIZE);
        else {
            if (file->ff & 128U) {
                file->ff = 3U;
                return leftover;
            }
            else {
                if (!sqz_readblksize(file->blk, file, file->libfmt))
                    goto error;
                if ( sqz_gztell(file) == file->size )
                    file->ff |= 128U;
                sqz_decode(sqz, blk, fmt, LOAD_SIZE);
            }
        }
        read = sqz->offset > ( len - leftover) ? (len - leftover) : sqz->offset;
        memcpy(outbuff + outpos, sqz->readbuffer, read);
        sqz->offset -= read;
        outpos += read;
        sqz->rem = read;
        if (outpos < len) {
            if (file->ff & 128U) {
                file->ff = 3U;
                return outpos;
            }
            else {
                if (!sqz_readblksize(file->blk, file, file->libfmt))
                    goto error;
                if (sqz_gztell(file) == file->size)
                    file->ff |= 128U;
                sqz_decode(sqz, blk, fmt, LOAD_SIZE);
                read = sqz->offset > (len-outpos) ? (len-outpos) : sqz->offset;
                memcpy(outbuff + outpos, sqz->readbuffer, read);
                sqz->offset -= read;
                sqz->rem = read;
                outpos += read;
            }
        }
        if (sqz->offset < len)
            file->ff = (file->ff & 128U) | 2U;
        else
            file->ff = (file->ff & 128U) | 1U;
        return outpos;
        }
    case 3:
        return 0;
    }
    error:
        return -1;
}

uint32_t sqz_getblocks(sqzFile sqzfp)
{
    if (!sqzfp) return 0;
    if (!(sqzfp->fmt & 4)) return 0;
    FILE *fp = fopen(sqzfp->namestr, "r");
    if (!fp) return 0;
    if (fseek(fp, -10, SEEK_END)) return 0;
    uint32_t n;
    if (!fread(&n, B32, 1, fp)) return 0;
    return n;
}

uint8_t sqz_go2blockn(sqzFile sqzfp, uint64_t n)
{
    if ( (n > sqz_getblocks(sqzfp)) ) return 1;
    uint64_t blkn = 0;
    uint64_t blks = 0;
    uint64_t blkr = 0;
    sqz_gzseek(sqzfp, HEADLEN, SEEK_SET);
    while (blkn < n) {
        sqz_gzseek(sqzfp, 8, SEEK_CUR);
        blkr = sqz_gzread(sqzfp, &blks, B64);
        if (blkr) {
            sqz_gzseek(sqzfp, blks, SEEK_CUR);
            sqzfp->filepos = sqz_gztell(sqzfp);
            blkn++;
        }
        else return 1;
    }
    return 0;
}

uint64_t sqz_loadblk(sqzblock_t *blk, uint8_t **b, uint64_t *s, uint8_t fq)
{
    uint64_t p = 0;
    uint64_t d = *s;
    uint8_t *o = *b;
    uint64_t w = 0;
    while (blk->newblk) {
        p = sqz_fastXdecode(blk, o + w, d - w, fq);
        if ( (p == d) || (0 == p) ) {
            d <<= 2;
            *b = realloc(*b, d);
            o = *b;
            if (!b) return 0;
        }
        w += p;
    }
    *s = d;
    return w;
}

sqzFile sqzopen(const char *filename, const char *mode)
{
    sqzFile sqzfp = calloc(1, sizeof(struct sqzFile_s));
    if (!sqzfp) return NULL;
    sqzfp->namestr = filename;
    if ( sqz_gzopen(filename, sqzfp, mode) ) {
        free(sqzfp);
        return NULL;
    }
    sqz_getformat(sqzfp);
    //If not sqz format, no need for the rest of members
    if ( !(sqzfp->fmt & 4) ) return sqzfp;

    sqzfp->libfmt = sqzfp->fmt >> 3;
    sqzfp->ff = 0;
    sqzfp->size = sqz_filesize(filename);
    sqz_gzseek(sqzfp, HEADLEN, SEEK_SET);
    sqzfp->filepos = sqz_gztell(sqzfp);
    sqzfp->blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!sqzfp->blk) {
        sqz_gzclose(sqzfp);
        free(sqzfp);
        return NULL;
    }
    sqzfp->sqz = sqz_fastxinit(sqzfp->fmt, LOAD_SIZE);
    if (!sqzfp->sqz) {
        sqz_gzclose(sqzfp);
        free(sqzfp);
        return NULL;
    }
    return sqzfp;
}

uint8_t sqz_hasdata(sqzfastx_t *sqz)
{
    //bit 7 is on if there is still data
    return sqz->datflag & 128;
}

void sqz_setnodata(sqzfastx_t *sqz)
{
    //Unset bit 7 from rest of threads
    sqz->datflag &= 127;
}

void sqz_setlastread(sqzfastx_t *sqz)
{
    //Set bit 1 of threads that got some data
    sqz->datflag |= 1;
}

uint8_t sqz_readend(sqzfastx_t *sqz)
{
    //bit 1 is on if thread had data, but reader has finished
    return sqz->datflag & 1;
}

void sqz_resetsqz(sqzfastx_t *sqz)
{
    sqz->namepos = 0;
    sqz->endflag = 0;
    sqz->blks++;
}

uint64_t sqz_getn(sqzfastx_t *sqz)
{
    return sqz->n;
}
