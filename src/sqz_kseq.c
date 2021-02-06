#include <stdio.h>
#include <zlib.h>
#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread)
#define SQZLIB
#include "sqz_kseq.h"


uint8_t  sqz_getformat(const char *filename)
{
    uint8_t ret = 0;
    if ( (ret = sqz_checksqz(filename)) ) return ret;

    gzFile fp = gzopen(filename, "r");
    if (!fp) return ret;

    kseq_t *seq = kseq_init(fp);
    if (!seq) {
        gzclose(fp);
        return ret;
    }
    int l = kseq_read(seq);
    //ERROR
    if (l < 0) {
        fprintf(stderr, "[sqzlib ERROR]: sequence file format not recognized\n");
        goto exit;
    }
    //FASTQ
    if (seq->qual.l > 0) {
        ret = 2;
        goto exit;
    }
    //FASTA
    ret = 1;
    exit:
        gzclose(fp);
        kseq_destroy(seq);
        return ret;
}


char     sqz_kseqinit(sqzfastx_t *sqz)
{
    char ret = 0;
    sqz->fp = gzopen(sqz->filename, "r");
    if (!sqz->fp) goto exit;
    sqz->seq = kseq_init(sqz->fp);
    if (!sqz->seq) {
        gzclose(sqz->fp);
        goto exit;
    }
    ret = 1;
    exit:
        return ret;
}


uint64_t sqz_loadfastq(sqzfastx_t *sqz)
{
    if (sqz->endflag) return sqz_fastqeblock(sqz);
    return sqz_fastqnblock(sqz);
}


uint64_t sqz_fastqnblock(sqzfastx_t *sqz)
{
    uint64_t offset = 0;
    uint64_t n = 0;
    uint64_t l;
    uint64_t maxlen = LOAD_SIZE - B64;
    //TODO remove pointer references inside loop
    while ( kseq_read(sqz->seq) >= 0 ) {
        l = sqz->seq->seq.l;
        sqz->bases += l;
        n++;
        if (!sqz_loadname(sqz, sqz->seq->name, n)) {
            offset = 0;
            goto exit;
        }
        //Determine if current sequence can be loaded completely
        if ( l + 1  > maxlen) {
            sqz->n = n;
            offset = sqz_fastqwrap(sqz, offset, maxlen);
            sqz->endflag = 1;
            goto exit;
        }
        memcpy(sqz->seqbuffer + offset, &l, B64);
        offset += B64;
        memcpy(sqz->seqbuffer + offset, sqz->seq->seq.s, l + 1);
        memcpy(sqz->qualbuffer + offset, sqz->seq->qual.s, l + 1);
        offset += l + 1;
        if ( maxlen <= l + 1 + B64 ) break;
        maxlen -= l + 1 + B64;
    }
    exit:
        sqz->n = n;
        sqz->offset = offset;
        return offset;
}


char     sqz_loadname(sqzfastx_t *sqz, kstring_t name, uint64_t n)
{
    char ret = 0;
    //fprintf(stderr, "\n\n%s\n\n", name.s);
    //sleep(2);
    memcpy(sqz->namebuffer + sqz->namepos, name.s, name.l + 1);
    sqz->namepos += name.l + 1;
    if (sqz->namepos + 100 >= sqz->namesize) {
        //fprintf(stderr, "At seq %lu\n", n);
        sqz->namebuffer = realloc(sqz->namebuffer, sqz->namesize*2);
        if (!(sqz->namebuffer)) {
            fprintf(stderr, "Name error\n");
            goto exit;
        }
        //fprintf(stderr, "Previous size: %lu\n", sqz->namesize);
        sqz->namesize *= 2;
        //fprintf(stderr, "Name memory %lu\n", sqz->namesize);
    }
    ret = 1;
    exit:
        return ret;
}


uint64_t sqz_fastqwrap(sqzfastx_t *sqz, uint64_t offset, uint64_t maxlen)
{
    //Compute how much buffer is available
    uint64_t l = sqz->seq->seq.l;
    //Copy sequence length data
    memcpy(sqz->seqbuffer + offset, &l, B64);
    offset += B64;
    //Copy as much seq data as we can fit in remaining buffer
    memcpy(sqz->seqbuffer + offset, sqz->seq->seq.s, maxlen);
    memcpy(sqz->qualbuffer + offset, sqz->seq->qual.s, maxlen);
    offset += maxlen;
    //Add null byte after loading data into buffers
    sqz->seqbuffer[offset] = '\0';
    sqz->qualbuffer[offset] = '\0';
    offset++;
    sqz->rem = l - maxlen;
    //Store length of sequence that could not complete loading
    sqz->prevlen = l;
    return offset;
}


uint64_t sqz_fastqeblock(sqzfastx_t *sqz)
{
    uint64_t l = sqz->seq->seq.l;
    //buffer can be completely filled with current sequence
    if (sqz->rem >= LOAD_SIZE) {
        memcpy(sqz->seqbuffer,
               sqz->seq->seq.s + (l + 1 - sqz->rem),
               LOAD_SIZE);
        memcpy(sqz->qualbuffer,
               sqz->seq->qual.s + (l + 1 - sqz->rem),
               LOAD_SIZE);
        sqz->rem -= LOAD_SIZE;
        sqz->offset = LOAD_SIZE;
        return LOAD_SIZE;
    }
    //Rest of sequence can go into buffer
    memcpy(sqz->seqbuffer,
           sqz->seq->seq.s + sqz->toread,
           sqz->rem + 1);
    memcpy(sqz->qualbuffer,
           sqz->seq->qual.s + sqz->toread,
           sqz->rem + 1);
    sqz->offset = sqz->rem + 1;
    sqz->endflag = 0;
    return sqz->offset;
}


uint64_t sqz_loadfasta(sqzfastx_t *sqz)
{
    if (sqz->endflag) return sqz_fastaeblock(sqz);
    return sqz_fastanblock(sqz);
}


uint64_t sqz_fastanblock(sqzfastx_t *sqz)
{
    uint64_t offset = 0;
    uint64_t n = 0;
    uint64_t l;
    uint64_t maxlen = LOAD_SIZE - B64;
    //TODO remove pointer references inside loop
    while ( (kseq_read(sqz->seq) >= 0) ) {
        l = sqz->seq->seq.l;
        sqz->bases += l;
        n++;
        if (!sqz_loadname(sqz, sqz->seq->name, n)) {
            offset = 0;
            goto exit;
        }
        //Determine if current sequence can be loaded completely
        if (l + 1 > maxlen) {
            sqz->n = n;
            offset = sqz_fastawrap(sqz, offset, maxlen);
            sqz->endflag = 1;
            goto exit;
        }
        memcpy(sqz->seqbuffer + offset, &l, B64);
        offset += B64;
        memcpy(sqz->seqbuffer+offset, sqz->seq->seq.s, l + 1);
        offset += l + 1;

        if ( maxlen <= l + 1 + B64 ) break;
        maxlen -= l + 1 + B64;
    }
    exit:
        sqz->n = n;
        sqz->offset = offset;
        return offset;
}


uint64_t sqz_fastawrap(sqzfastx_t *sqz, uint64_t offset, uint64_t maxlen)
{
    //Compute how much buffer is available
    uint64_t l = sqz->seq->seq.l;
    //Copy sequence length data
    memcpy(sqz->seqbuffer + offset, &l, B64);
    offset += B64;
    //Copy as much seq data as we can fit in remaining buffer
    memcpy(sqz->seqbuffer + offset, sqz->seq->seq.s, maxlen);

    offset += maxlen;
    //Add null byte after loading data into buffers
    sqz->seqbuffer[offset++] = '\0';


    sqz->rem = l - maxlen;
    //Store length of sequence that could not complete loading
    sqz->prevlen = l;
    return offset;
}


uint64_t sqz_fastaeblock(sqzfastx_t *sqz)
{
    uint64_t l = sqz->seq->seq.l;
    //buffer can be completely filled with current sequence
    if (sqz->rem >= LOAD_SIZE) {
        fprintf(stderr, ">>>case rem %lu\n", sqz->rem);
        memcpy(sqz->seqbuffer,
               sqz->seq->seq.s + (l + 1 - sqz->rem),
               LOAD_SIZE);



        sqz->rem -= LOAD_SIZE;
        sqz->offset = LOAD_SIZE;
        return LOAD_SIZE;
    }
    //Rest of sequence can go into buffer
    memcpy(sqz->seqbuffer,
           sqz->seq->seq.s + sqz->toread,
           sqz->rem + 1);



    sqz->offset = sqz->rem + 1;
    sqz->endflag = 0;
    return sqz->offset;
}


void     sqz_kill(sqzfastx_t *sqz)
{
    if (sqz) {
        if (sqz->fp)
            gzclose(sqz->fp);
        if (sqz->seq)
            kseq_destroy(sqz->seq);
        free(sqz->seqbuffer);
        free(sqz->namebuffer);
        free(sqz->readbuffer);
        if ( (sqz->fmt == 2) | (sqz->fmt == 14) ) free(sqz->qualbuffer);
        free(sqz);
    }
}


void     sqz_killblk(sqzblock_t *blk)
{
    free(blk->blkbuff);
    free(blk->cmpbuff);
    free(blk);
}


unsigned char sqz_checksqz(const char *filename)
{
    size_t tmp = 0;
    FILE *fp = fopen(filename, "rb");
    if (!fp) return 0;
    uint32_t magic;
    unsigned char fmt = 0;
    char sqz;
    //Read magic
    tmp += fread(&magic, 1, 4, fp);
    if (MAGIC ^ magic) {
        //Return 0 if not an sqz file
        fclose(fp);
        return 0;
    }
    //Set sqz flag
    fmt |= 4;
    //Read format (fastA or fastQ)
    tmp += fread(&sqz, 1, 1, fp);
    fmt |= sqz;
    //Read compression library (zlib or ...)
    tmp += fread(&sqz, 1, 1, fp);
    fmt |= sqz << 3;
    fclose(fp);
    return fmt;
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

    sqzfile.sqz = sqz_fastxinit(filename, LOAD_SIZE);
    if ( !sqzfile.sqz | !(sqzfile.sqz->fmt | 4) ) {
        sqz_kill(sqzfile.sqz);
        sqzfile.sqz = NULL;
        goto exit;
    }

    sqzfile.blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!sqzfile.blk) {
        sqz_kill(sqzfile.sqz);
        sqzfile.sqz = NULL;
        sqzfile.blk = NULL;
        goto exit;
    }
    exit:
        return sqzfile;
}


int64_t  sqz_sqzread(sqz_File *file, void *buff, size_t len)
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
        if (ftell(file->fp) == file->size) {
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
                if (ftell(file->fp) == file->size) {
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


void     sqz_sqzclose(sqz_File file)
{
    sqz_kill(file.sqz);
    sqz_killblk(file.blk);
    fclose(file.fp);
}


char     sqz_readblksize(sqzblock_t *blk, FILE *fp)
{
    //TODO Error handling
    char ret = 0;
    uint64_t cmpsize;
    uint64_t dcpsize;
    uint64_t cbytes;
    uint64_t nelem;
    nelem =  fread(&dcpsize, B64, 1, fp);
    nelem += fread(&cmpsize, B64, 1, fp);
    if ((cmpsize != fread(blk->cmpbuff, 1, cmpsize, fp)) | (nelem != 2)) {
        //TODO Error goes here
        goto exit;
    }
    blk->cmpsize = cmpsize;
    blk->blksize = LOAD_SIZE*2;
    cbytes = sqz_inflate(blk);
    if (cbytes != dcpsize) {
        //TODO Error goes here
        goto exit;
    }
    blk->blksize = dcpsize;
    ret = 1;
    exit:
        return ret;
}


void     sqz_decode(sqzfastx_t *sqz, sqzblock_t *blk, uint64_t klibl)
{
    switch (sqz->fmt) {
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
