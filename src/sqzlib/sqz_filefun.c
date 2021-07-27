#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "sqz_data.h"

char zbytes[4] = {0, 0, 0, 0};
char magic[4] = {5, 8, 5, 9};
unsigned char cmpflag = 1;

sqzblock_t *sqz_sqzblkinit(uint64_t size);
void sqz_fastxkill(sqzfastx_t *sqz);
void sqz_blkkill(sqzblock_t *blk);
uint64_t sqz_inflate(sqzblock_t *blk);
uint64_t sqz_fastXdecode(sqzblock_t *blk,
                         uint8_t *buff,
                         uint64_t buffsize,
                         uint8_t fqflag);
//uint64_t sqz_getblkbytes(uint8_t *blkbuff, uint64_t seqlength)
//{
//    uint64_t bases = 0;
//    uint64_t blkpos = 0;
//    while (bases != seqlength) {
//        bases += sqz_gotoblk(blkbuff + blkpos, &blkpos, 0, 0);
//    }
//   return blkpos;
//}


int64_t sqz_filesize(FILE *fp)
{
    fseek(fp, 0, SEEK_END);
    int64_t s = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    return s - 16;
}


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


void sqzrewind(sqz_File *sqzfp)
{
    sqz_fastxreset(sqzfp->sqz);
    sqz_blkreset(sqzfp->blk);
    fseek(sqzfp->fp, HEADLEN, SEEK_SET);
    sqzfp->filepos = ftell(sqzfp->fp);
    sqzfp->ff = 0;
}


char sqz_filehead(unsigned char fmt, FILE *ofp)
{
    char wbytes = 0;
    if ( 4 != (wbytes += fwrite(magic, 1, 4, ofp)) ) return 0;
    if ( 5 != (wbytes += fwrite(&fmt,  1, 1, ofp)) ) return 0;
    //Compression library
    if ( 6 != (wbytes += fwrite(&cmpflag, 1, 1, ofp)) ) return 0;
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


sqz_File sqz_sqzopen(char *filename)
{
    sqz_File sqzfile = {NULL, NULL, NULL, 0, 0, 0};
    sqzfile.fp = fopen(filename, "rb");
    if (!sqzfile.fp) {
        sqzfile.fp = NULL;
        goto exit;
    }
    sqzfile.size = sqz_filesize(sqzfile.fp);
    fseek(sqzfile.fp, HEADLEN, SEEK_SET);
    sqzfile.filepos = ftell(sqzfile.fp);

    //TODO Adjust comented code to reflect changes in sqzfastx_t struct
    //TODO Code was comented while reworking code to use multiple threads
    //sqzfile.sqz = sqz_fastxinit(filename, LOAD_SIZE);
    //if ( !sqzfile.sqz || !(sqzfile.sqz->fmt & 4) ) {
    //    sqz_kill(sqzfile.sqz);
    //    sqzfile.sqz = NULL;
    //    goto exit;
    //}

    sqzfile.blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!sqzfile.blk) {
        sqz_fastxkill(sqzfile.sqz);
        sqzfile.sqz = NULL;
        sqzfile.blk = NULL;
        goto exit;
    }
    exit:
        return sqzfile;
}


char sqz_readblksize(sqzblock_t *blk, FILE *fp)
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
    if (dcpsize != sqz_inflate(blk))
        goto exit;
    blk->newblk = 1;
    ret = 1;
    exit:
        return ret;
}


static void sqz_decode(sqzfastx_t *sqz, sqzblock_t *blk, uint64_t klibl)
{
    //TODO The changes below are non functionsl and only thereto permit building
    //while reworking code to use multiple threads
    int fmt = 1;
    switch (fmt) {
        //switch (sqz->fmt) {
    case 14:
        {
            sqz->offset = sqz_fastXdecode(blk, sqz->readbuffer, klibl, 1);
            break;
        }
    case 13:
        {
            sqz->offset = sqz_fastXdecode(blk, sqz->readbuffer, klibl, 0);
            break;
        }
    }
}


int64_t sqz_sqzread(sqz_File *file, void *buff, size_t len)
{
    if (!file | !buff) return 0;
    sqzfastx_t *sqz  = file->sqz;
    sqzblock_t *blk  = file->blk;
    uint8_t *outbuff = (uint8_t *)buff;
    int64_t read     = 0;
    //Take lower 7 bits of flag
    switch (file->ff & 127) {
    case 0: //Buffers are empty
        {
        // Decompress sqz block
        if (!sqz_readblksize(file->blk, file->fp)) goto error;
        // Check if we have reached end of file
        if (ftell(file->fp) == (long)file->size) {
            //Set bit 7
            file->ff |= 128;
        }
        // Decode sqz block
        sqz_decode(sqz, blk, LOAD_SIZE);
        // Compute how much data can be loaded
        read = sqz->offset > len ? len : sqz->offset;
        // Load data
        memcpy(outbuff, sqz->readbuffer, read);
        // Update how much has been read
        sqz->rem += read;
        // Update how much data is left
        sqz->offset -= read;
        if (sqz->offset < len) {
            //All of leftover data fits in buffer. So more data will be needed
            //or reading has finished
            //Keep bit 7 status, switch flag to 2
            file->ff = (file->ff & 128) | 2;
        }
        else {
            //Can keep loading data
            //Keep bit 7 status, switch flag to 1
            file->ff = (file->ff & 128) | 1;
        }
        //fprintf(stderr, "Returning case 0 %ld\n", read);
        return read;
        }
    case 1: //Data just need to be copied
        {
        read = sqz->offset > len ? len : sqz->offset;
        memcpy(outbuff, sqz->readbuffer + sqz->rem, read);
        sqz->rem += read;
        sqz->offset -= read;
        if (sqz->offset < len) {
            //Keep bit 7 status, switch flag to 2
            file->ff = (file->ff & 128) | 2;
        }
        return read;
        }
    case 2: //Data can be copied but entire buffer can't be filled
        {
        //Store how much data needs to be copied
        uint64_t leftover = sqz->offset;
        //Finish copying remaining data
        memcpy(outbuff, sqz->readbuffer + sqz->rem, leftover);
        //Check if there is more data to decode
        if (blk->blkpos) {
            //There is more data to decode
            sqz_decode(sqz, blk, LOAD_SIZE);
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
                if (!sqz_readblksize(file->blk, file->fp)) goto error;
                if (ftell(file->fp) == (long)file->size) {
                    //Set bit 7
                    file->ff |= 128;
                }
                sqz_decode(sqz, blk, LOAD_SIZE);
            }
        }
        //Determine how much new data can be copied
        read = sqz->offset > ( len - leftover) ? (len - leftover) : sqz->offset;
        //Copy data making sure to not overwrite previously copied data
        memcpy(outbuff + leftover, sqz->readbuffer, read);
        sqz->offset -= read;
        sqz->rem = read;
        if (sqz->offset < len)
            //Will need to exit or more decoding
            file->ff = (file->ff & 128) | 2;
        else
            //Can keep loading data
            file->ff = (file->ff & 128) | 1;
        return read + leftover;
        }
    case 3:
        return 0;
    }
    error:
        return -1;
}


void sqz_sqzclose(sqz_File file)
{
    sqz_fastxkill(file.sqz);
    sqz_blkkill(file.blk);
    fclose(file.fp);
}
