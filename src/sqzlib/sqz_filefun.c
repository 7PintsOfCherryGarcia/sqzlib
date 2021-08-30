#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "sqz_data.h"

char zbytes[4] = {0, 0, 0, 0};
char magic[4] = {5, 8, 5, 9};
unsigned char cmpflag = 1;

sqzfastx_t *sqz_fastxinit(uint8_t fmt, uint64_t size);
sqzblock_t *sqz_sqzblkinit(uint64_t size);
void sqz_fastxkill(sqzfastx_t *sqz);
void sqz_blkkill(sqzblock_t *blk);
uint64_t sqz_inflate(sqzblock_t *blk);
uint64_t sqz_fastXdecode(sqzblock_t *blk,
                         uint8_t *buff,
                         uint64_t buffsize,
                         uint8_t fqflag);
uint8_t sqz_getformat(const char *filename);
FILE *fdopen(int fd, const char *mode);
sqzFile sqz_gzopen(const char *filename, sqzFile sqzfp, const char *mode);
int64_t sqz_gzread(gzFile fp, void *buff, uint32_t len);
size_t sqz_deflate(sqzblock_t *blk, int level);
int64_t sqz_zstdcompress(sqzblock_t *blk, int level);
uint64_t sqz_zstddecompress(sqzblock_t *blk);

static void sqz_fastxreset(sqzfastx_t *sqz)
{
    sqz->endflag    = 0;
    sqz->cmpflag    = 0;
    sqz->offset     = 0;
    sqz->namepos    = 0;
    sqz->n          = 0;
    sqz->bases      = 0;
    sqz->rem        = 0;
    sqz->toread     = 0;
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


int64_t sqzcompress(sqzblock_t *blk, int level, uint8_t libfmt)
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
    fseek(sqzfp->fp, HEADLEN, SEEK_SET);
    sqzfp->filepos = ftell(sqzfp->fp);
    sqzfp->ff = 0;
}


char sqz_filehead(uint8_t fmt, uint8_t libfmt, FILE *ofp)
{
    char wbytes = 0;
    if ( 4 != (wbytes += fwrite(magic, 1, 4, ofp)) ) return 0;
    if ( 5 != (wbytes += fwrite(&fmt,  1, 1, ofp)) ) return 0;
    //Compression library
    if ( 6 != (wbytes += fwrite(&libfmt, 1, 1, ofp)) ) return 0;
    if ( 8 != (wbytes += fwrite(zbytes,   1, 2, ofp)) ) return 0;
    return wbytes;
}


char sqz_filetail(uint64_t numseqs, FILE *ofp)
{
    if ( 4 != fwrite(zbytes, 1, 4, ofp) ) {
        return 0;
    }
    if ( 1 != fwrite(&numseqs, sizeof(numseqs), 1, ofp) ) {
        return 0;
    }
    if ( 4 != fwrite(zbytes, 1, 4, ofp) ) {
        return 0;
    }
    return 1;
}


char sqz_blkdump(sqzblock_t *blk, uint64_t size, FILE *ofp)
{
    size_t wbytes = 0;
    //Write uncompressed number of bytes in block
    wbytes += fwrite(&(blk->blkpos), B64, 1, ofp);
    //Write compressed number of bytes in block
    wbytes += fwrite(&(size), B64, 1, ofp);
    //Write block
    wbytes += fwrite(blk->cmpbuff, 1, size, ofp);
    if (wbytes != size+2) return 0;
    return 1;
}


int64_t sqz_filesize(FILE *fp)
{
    fseek(fp, 0, SEEK_END);
    int64_t s = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    return s - 16;
}


sqzFile sqzopen(const char *filename, const char *mode)
{
    sqzFile sqzfp = calloc(1, sizeof(struct sqzFile_s));
    sqzfp->ff = 0;
    sqzfp->fmt = sqz_getformat(filename);
    fprintf(stderr, "From sqzopen: %u\n", sqzfp->fmt);
    //Check if not sqz format (bit 6)
    if ( !(sqzfp->fmt & 4) )
        return sqz_gzopen(filename, sqzfp, mode);
    //Check compression library
    sqzfp->libfmt = sqzfp->fmt >> 3;
    sqzfp->fp = fopen(filename, "rb");
    if (!sqzfp->fp) {
        free(sqzfp);
        return NULL;
    }
    sqzfp->size = sqz_filesize(sqzfp->fp);
    fseek(sqzfp->fp, HEADLEN, SEEK_SET);
    sqzfp->filepos = ftell(sqzfp->fp);
    sqzfp->blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!sqzfp->blk) {
        //TODO: Let killer function take care of this
        fclose(sqzfp->fp);
        free(sqzfp);
        return NULL;
    }
    sqzfp->sqz = sqz_fastxinit(sqzfp->fmt, LOAD_SIZE);
    if (!sqzfp->sqz) {
        //TODO: Let killer function take care of this
        fclose(sqzfp->fp);
        free(sqzfp);
        return NULL;
    }
    return sqzfp;
}


sqzFile sqzdopen(int fd, const char *mode)
{
    return NULL;
}


void sqzclose(sqzFile file)
{
    sqz_fastxkill(file->sqz);
    sqz_blkkill(file->blk);
    (file->fmt & 4) ? fclose(file->fp) : gzclose(file->gzfp);
    free(file);
}


char sqz_readblksize(sqzblock_t *blk, FILE *fp, uint8_t libfmt)
{
    char ret = 0;
    uint64_t cmpsize;
    uint64_t dcpsize;
    uint64_t nelem;
    nelem =  fread(&dcpsize, B64, 1, fp);
    nelem += fread(&cmpsize, B64, 1, fp);
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
    if ( (cmpsize != fread(blk->cmpbuff, 1, cmpsize, fp)) || (nelem != 2) )
        goto exit;
    if (dcpsize != sqzdecompress(blk, libfmt))
        goto exit;
    blk->newblk = 1;
    ret = 1;
    exit:
        return ret;
}


int64_t sqzread(sqzFile file, void *buff, uint64_t len)
{
    if (!file | !buff) return -1;
    if ( !(file->fmt & 4) ) return sqz_gzread(file->gzfp, buff, (uint32_t)len);
    sqzfastx_t *sqz  = file->sqz;
    sqzblock_t *blk  = file->blk;
    uint8_t     fmt  = file->fmt & 3;
    uint8_t *outbuff = (uint8_t *)buff;
    int64_t read     = 0;
    //Take lower 7 bits of flag
    switch (file->ff & 127) {
    case 0: //Buffers are empty
        // Decompress sqz block
        if (!sqz_readblksize(file->blk, file->fp, file->libfmt)) goto error;
        // Check if we have reached end of file
        if (ftell(file->fp) == (long)file->size) {
            //Set bit 7
            file->ff |= 128;
        }
        // Decode sqz block
        sqz_decode(sqz, blk, fmt, LOAD_SIZE);
        // Compute how much data can be loaded
        read = sqz->offset > len ? len : sqz->offset;
        // Load data
        memcpy(outbuff, sqz->readbuffer, read);
        // Update how much has been read
        sqz->rem += read;
        // Update how much data is left
        sqz->offset -= read;
        if (sqz->offset < len) {
            //All of leftover data will fit in buffer.
            //So more data will be needed
            //or reading has finished
            //Keep bit 7 status, switch bit 2
            file->ff = (file->ff & 128) | 2;
        }
        else {
            //Can keep loading data
            //Keep bit 7 status, switch flag to 1
            file->ff = (file->ff & 128) | 1;
        }
        return read;
    case 1: //Data just need to be copied
        read = sqz->offset > len ? len : sqz->offset;
        memcpy(outbuff, sqz->readbuffer + sqz->rem, read);
        sqz->rem += read;
        sqz->offset -= read;
        if (sqz->offset < len) {
            //All of leftover data will fit in buffer.
            //So more data will be needed
            //or reading has finished
            //Keep bit 7 status, switch bit 2
            file->ff = (file->ff & 128) | 2;
        }
        return read;
    case 2: //Data can be copied but entire buffer can't be filled
        ;
        uint64_t outpos;
        //Store how much data needs to be copied
        uint64_t leftover = sqz->offset;
        //Finish copying remaining data aka empty sqz->readbuffer
        memcpy(outbuff, sqz->readbuffer + sqz->rem, leftover);
        outpos = leftover;
        //Check if there is more data to decode
        if (blk->blkpos) {
            //There is more data to decode
            sqz_decode(sqz, blk, fmt, LOAD_SIZE);
        }
        else {
            //There is no more data to decode
            //We need to know if more data exists in file if there is more data,
            //data needs to be decompressed and decoded. Otherwise, reading data
            //has finished and we can exit
            if (file->ff & 128) {
                //We can terminate
                file->ff = 3;
                return leftover;
            }
            else {
                if (!sqz_readblksize(file->blk, file->fp, file->libfmt))
                    goto error;
                if (ftell(file->fp) == (long)file->size) {
                    //Set bit 7
                    file->ff |= 128;
                }
                sqz_decode(sqz, blk, fmt, LOAD_SIZE);
            }
        }
        //Determine how much new data can be copied
        read = sqz->offset > ( len - leftover) ? (len - leftover) : sqz->offset;
        //Copy data making sure to not overwrite previously copied data
        memcpy(outbuff + outpos, sqz->readbuffer, read);
        sqz->offset -= read;
        outpos += read;
        sqz->rem = read;
        if (outpos < len) {
            //Buffer still not full, end of file?
            //There is no more data to decode
            //We need to know if more data exists in file if there is more data,
            //data needs to be decompressed and decoded. Otherwise, reading data
            //has finished and we can exit
            if (file->ff & 128) {
                //We can terminate
                file->ff = 3;
                return outpos;
            }
            else {
                if (!sqz_readblksize(file->blk, file->fp, file->libfmt))
                    goto error;
                if (ftell(file->fp) == (long)file->size) {
                    //Set bit 7
                    file->ff |= 128;
                }
                sqz_decode(sqz, blk, fmt, LOAD_SIZE);
                read = sqz->offset > (len-outpos) ? (len-outpos) : sqz->offset;
                memcpy(outbuff + outpos, sqz->readbuffer, read);
                sqz->offset -= read;
                sqz->rem = read;
                outpos += read;
            }
        }
        if (sqz->offset < len) {
            //We remain in case 2
            file->ff = (file->ff & 128) | 2;
        }
        else
            file->ff = (file->ff & 128) | 1;
        return outpos;
    case 3:
        return 0;
    }
    error:
        return -1;
}
