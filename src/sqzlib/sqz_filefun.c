#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include <unistd.h>

#include "sqz_data.h"

const char zbytes[4] = {0, 0, 0, 0};
const uint8_t magic1[4] = {5, 8, 5, 9};
const uint8_t magic2[4] = {9, 5, 8, 5};

unsigned char cmpflag = 1;

static inline uint64_t sqz_getnamesize(uint8_t *buff, uint64_t size)
{
    return *(uint64_t *)( buff + ( size - B64 ) );
}

static void sqz_blkreset(sqzblock_t *blk)
{
    blk->blkbuff->pos  = 0;
    blk->newblk  = 1;
    blk->cmpbuff->pos  = 0;
}

static void sqz_fastxreset(sqzfastx_t *sqz)
{
    sqz_blkreset(sqz->blk);
    sqz->lseqflag = 0;
    sqz->cmpflag  = 0;
    sqz->offset   = 0;
    sqz->n        = 0;
    sqz->bases    = 0;
}

static uint64_t sqz_datasize(const char *filename)
{
    if (!filename) return 0;
    FILE *fp = fopen(filename, "r");
    if ( fseek(fp, 0, SEEK_END) ) return 0;
    uint64_t s = ftell(fp);
    fclose(fp);
    return s - 24;
}

uint64_t sqz_blkcompress(sqzfastx_t *sqz, int level, uint8_t libfmt)
{
    uint64_t cmpbytes = 0;
    switch (libfmt) {
        case 1:
            if ( !(cmpbytes = sqz_deflate(sqz->blk, level)) )
                return 0;
            return cmpbytes;
        case 2:
            return sqz_zstdcompress(sqz->blk, level);
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
    sqz_gzseek(sqzfp, HEADLEN, SEEK_SET);
    sqzfp->decoded = 0;
    sqzfp->bloaded = 0;
    sqzfp->rflag   = 0;
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

uint8_t sqz_filetail(uint64_t numseqs, uint32_t nblocks, FILE *ofp)
{
    if ( 4 != fwrite(zbytes,     1, 4, ofp) ) return 0;
    if ( 1 != fwrite(&numseqs, B64, 1, ofp) ) return 0;
    if ( 2 != fwrite(zbytes,     1, 2, ofp) ) return 0;
    if ( 1 != fwrite(&nblocks, B32, 1, ofp) ) return 0;
    if ( 2 != fwrite(zbytes,     1, 2, ofp) ) return 0;
    if ( 4 != fwrite(magic2,     1, 4, ofp) ) return 0;
    return 1;
}

uint8_t sqz_blkdump(sqzfastx_t *sqz, uint64_t cmpsize, FILE *ofp)
{
    size_t wbytes = 0;
    //Write uncompressed number of bytes in block
    wbytes += fwrite(&(sqz->blk->blkbuff->pos), B64, 1, ofp);
    //Write compressed number of bytes in block
    wbytes += fwrite(&cmpsize, B64, 1, ofp);
    //Write block
    wbytes += fwrite(sqz->blk->cmpbuff->data, cmpsize, 1, ofp);
    if (wbytes != 3) return 0;
    return 1;
}

void sqz_resetblk(sqzblock_t *blk)
{
    blk->blkbuff->pos = 0;
}

uint8_t sqz_newblk(sqzblock_t *blk)
{
    return blk->newblk;
}

//uint64_t sqz_getfilepos(sqzFile sqzfp)
//{
//    return ( sqzfp->filepos = sqz_gztell(sqzfp) );
//}

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
    return sqzfp->sqz->blk;
}

uint8_t sqz_format(sqzFile sqzfp)
{
    return sqzfp->fmt;
}

uint8_t sqz_readblksize(sqzFile sqzfp)
{
    uint8_t  ret = 1;
    sqzblock_t *blk = sqzfp->sqz->blk;
    sqzbuff_t *cmpbuff = blk->cmpbuff;
    sqzbuff_t *blkbuff = blk->blkbuff;
    uint64_t cmpsize;
    uint64_t dcpsize;
    uint64_t nelem;
    nelem = sqz_gzread(sqzfp, &dcpsize, B64);
    nelem += sqz_gzread(sqzfp, &cmpsize, B64);
    if ( nelem != 16 ) goto exit;
    if ( cmpsize > cmpbuff->size ) {
        cmpbuff = sqz_buffrealloc(cmpbuff, cmpsize);
        if ( !cmpbuff ) goto exit;
        blk->cmpbuff = cmpbuff;
    }
    if ( dcpsize > blkbuff->size ) {
        blkbuff = sqz_buffrealloc(blkbuff, dcpsize);
        if (!blkbuff) goto exit;
        blk->blkbuff = blkbuff;
    }

    cmpbuff->pos = cmpsize;
    blkbuff->pos = dcpsize;
    if ( (cmpsize != (uint64_t)sqz_gzread(sqzfp, cmpbuff->data, cmpsize) ) )
        goto exit;
    if (dcpsize != sqzdecompress(blk, sqz_fileformat(sqzfp)))
        goto exit;
    sqzbuff_t *names = sqzfp->sqz->namebuffer;
    uint64_t namesize = sqz_getnamesize(blkbuff->data, blkbuff->pos);
    uint64_t  datasize  = blkbuff->pos - B64 - namesize;
    if (names->size < namesize) {
        names->data = realloc(names->data, namesize + 1);
        if (names->data) goto exit;
        names->size = namesize;
    }
    names->data = memcpy(names->data, (uint8_t *)blkbuff->data + datasize, namesize);
    names->pos  = namesize;
    blk->datasize = blkbuff->pos - B64 - namesize;
    blk->newblk = 1;
    blkbuff->pos = 0;
    ret = 0;
    exit:
        return ret;
}

uint32_t sqz_getblocks(sqzFile sqzfp)
{
    //TODO change header format to include number of blocks
    //this is due to gzseek not supporting SEEK_END
    if (!sqzfp) return 0;
    if (!(sqzfp->fmt & 4)) return 0;
    FILE *fp = fopen(sqzfp->name, "r");
    if (!fp) return 0;
    if (fseek(fp, -10, SEEK_END)) return 0;
    uint32_t n;
    if (!fread(&n, B32, 1, fp)) return 0;
    return n;
}

uint8_t sqz_go2blockn(sqzFile sqzfp, uint64_t n)
{
    if ( (n > sqzfp->nblocks) ) return 1;
    uint64_t blkn = 0;
    uint64_t blks = 0;
    uint64_t blkr = 0;
    sqz_gzseek(sqzfp, HEADLEN, SEEK_SET);
    while (blkn < n) {
        sqz_gzseek(sqzfp, 8, SEEK_CUR);
        blkr = sqz_gzread(sqzfp, &blks, B64);
        if (blkr) {
            sqz_gzseek(sqzfp, blks, SEEK_CUR);
            //sqzfp->filepos = sqz_gztell(sqzfp);
            blkn++;
        }
        else return 1;
    }
    return 0;
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
    sqz_resetblk(sqz->blk);
    sqz->lseqbuff->pos = 0;
    sqz->namebuffer->pos = 0;
    sqz->lseqflag = 0;
    sqz->blks++;
}

uint64_t sqz_getn(sqzfastx_t *sqz)
{
    return sqz->n;
}

uint8_t sqz_loadblockn(sqzFile sqzfp, uint32_t n)
{
    uint8_t ret = 1;
    if ( sqz_go2blockn(sqzfp, n) ) goto exit;
    if ( sqz_readblksize(sqzfp) )  goto exit;
    ret = 0;
    exit:
        return ret;
}

uint8_t sqz_emptysqzfp(sqzFile sqzfp, uint8_t *buff)
{
    if (!buff) return 1;
    uint8_t *fpbuff = sqzfp->sqz->readbuffer->data;
    uint64_t size   = sqzfp->sqz->readbuffer->pos;
    memcpy(buff, fpbuff, size);
    sqzfp->sqz->readbuffer->pos = 0;
    return 0;
}

uint8_t sqz_sqzfp2buff(sqzFile sqzfp, uint8_t *buff, uint64_t size)
{
    if (!buff) return 1;
    uint8_t *fpbuff = sqzfp->sqz->readbuffer->data;
    uint32_t l = sqzfp->bloaded;
    memcpy(buff, fpbuff + l, size);
    sqzfp->bloaded += size;
    return 0;

}

sqzFile sqzopen(const char *filename, const char *mode)
{
    sqzFile sqzfp = calloc(1, sizeof(struct sqzFile_s));
    if (!sqzfp) return NULL;
    sqzfp->name = filename;
    if ( sqz_gzopen(filename, sqzfp, mode) ) {
        free(sqzfp);
        return NULL;
    }
    sqz_getformat(sqzfp);
    sqzfp->kseq = sqz_kseqinit(sqzfp);
    if ( !sqzfp->kseq ) {
        free(sqzfp);
        return NULL;
    }
    //If not sqz format, no need for the rest of members
    if ( !(sqzfp->fmt & 4) ) return sqzfp;

    sqzfp->sqz = sqz_fastxinit(sqzfp->fmt, LOAD_SIZE);
    if ( !sqzfp->sqz ) {
        sqz_gzclose(sqzfp);
        free(sqzfp);
        return NULL;
    }
    sqzfp->dataend = sqz_datasize(filename);
    sqzfp->nblocks = sqz_getblocks(sqzfp);
    sqzfp->rflag   = 0;
    sqzfp->decoded = 0;
    sqzfp->bloaded = 0;
    sqz_gzseek(sqzfp, HEADLEN, SEEK_SET);
    return sqzfp;
}

int64_t sqzread(sqzFile sqzfp, void *buff, uint64_t len)
{
    if (!sqzfp || !buff) return -1;
    if ( !(sqzfp->fmt & 4) )
        return sqz_gzread(sqzfp, buff, (uint32_t)len);
    uint8_t  *outbuff = buff;
    uint64_t read     = 0;
    uint64_t tocpy    = 0;
    switch (sqzfp->rflag) {
    case 0:
        {
        if ( sqz_readblksize(sqzfp) )
            goto error;
        read = sqz_decode(sqzfp);
        tocpy = read > len ? len : read;
        sqz_sqzfp2buff(sqzfp, outbuff, tocpy);
        if ( read < len )
            sqzfp->rflag = (sqzfp->rflag & 128) | 2U;
        else
            sqzfp->rflag = (sqzfp->rflag & 128) | 1U;
        return (int64_t)tocpy;
        }
    case 1: //Data just need to be copied
        {
        read = sqzfp->decoded - sqzfp->bloaded;
        tocpy = read > len ? len : read;
        sqz_sqzfp2buff(sqzfp, outbuff, tocpy);
        read = sqzfp->decoded - sqzfp->bloaded;
        if ( read < len )
            sqzfp->rflag = (sqzfp->rflag & 128) | 2U;
        return (int64_t)tocpy;
        }
    case 2: //Data can be copied but entire buffer can't be filled
        memset(outbuff, 0, len);
        tocpy = sqzfp->decoded - sqzfp->bloaded;
        sqz_sqzfp2buff(sqzfp, outbuff, tocpy);
        sqzfp->sqz->readbuffer->pos = 0;
        sqzfp->rflag = (sqzfp->rflag & 128) | 0U;
        if ( sqz_gztell(sqzfp) == sqzfp->dataend)
            sqzfp->rflag = 3U;
        return (int64_t)tocpy;
    case 3:
        return 0;
    }
    error:
        return -1;
}

void sqzclose(sqzFile sqzfp)
{
    sqz_gzclose(sqzfp);
    sqz_kseqdestroy(sqzfp);
    if (sqzfp->sqz) sqz_fastxkill(sqzfp->sqz);
    free(sqzfp);
}
