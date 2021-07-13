#include <zlib.h>
#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread)
#include "sqz_kseq.h"


uint8_t  sqz_getformat(const char *filename)
{
    uint8_t ret = 0;
    if ( (ret = sqz_checksqz(filename)) ) return ret;
    ret = 0;
    //TODO: Adjust to other compression libraries
    gzFile fp = gzopen(filename, "r");
    if (!fp) return ret;

    kseq_t *seq = kseq_init(fp);
    if (!seq) {
        gzclose(fp);
        return ret;
    }
    int l = kseq_read(seq);

    //TODO: Add error signaling as this is a library so it should not print
    //ERROR
    if (l < 0)
        goto exit;

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


uint8_t  sqz_checksqz(const char *filename)
{
    size_t   tmp = 0;
    uint32_t magic;
    uint8_t  fmt = 0;
    char     sqz;

    FILE *fp = fopen(filename, "rb");
    if (!fp) return 0;
    //Read magic
    tmp += fread(&magic, 1, 4, fp);
    if (MAGIC ^ magic) {
        //Return 0 if not an sqz file
        fclose(fp);
        return 0;
    }
    //Set sqz flag (bit 6)
    fmt |= 4;

    //TODO: Check for valid entries. A file may have the magic number
    //but still not be an sqz file

    //Read format: 1 - fastA or 2 - fastQ (bit 5)
    tmp += fread(&sqz, 1, 1, fp);
    fmt |= sqz;
    //Read compression library: 1 - zlib (bits 4, 3, and 2)
    tmp += fread(&sqz, 1, 1, fp);
    fmt |= sqz << 3;
    fclose(fp);
    return fmt;
}


uint64_t sqz_loadfastX(sqzfastx_t *sqz, uint8_t fqflag, kseq_t *seq)
{
    if (sqz->endflag) {
        if (fqflag) return sqz_fastqeblock(sqz);
        return sqz_fastaeblock(sqz);
    }
    if (fqflag) return sqz_fastqnblock(sqz, seq);
    return sqz_fastanblock(sqz, seq);
}


static uint8_t sqz_loadname(sqzfastx_t *sqz, kseq_t *seq)
{
    uint8_t ret = 0;
    uint8_t *namebuffer = sqz->namebuffer;
    uint64_t pos = sqz->namepos;
    memcpy(namebuffer + pos, seq->name.s, seq->name.l + 1);
    pos += seq->name.l + 1;
    //If comment exists
    if (seq->comment.s) {
        //Substitute terminating null with space
        namebuffer[pos - 1] = ' ';
        //Append comment including terminating null
        memcpy(namebuffer + pos, seq->comment.s, seq->comment.l + 1);
        pos += seq->comment.l + 1;
    }
    if (pos + 100 >= sqz->namesize) {
        namebuffer = realloc(namebuffer, sqz->namesize*2);
        if (!(namebuffer))
            goto exit;
        sqz->namebuffer = namebuffer;
        sqz->namesize *= 2;
    }
    ret = 1;
    exit:
        sqz->namepos = pos;
        return ret;
}


static uint64_t sqz_fastqnblock(sqzfastx_t *sqz, kseq_t *seq)
{
    uint64_t offset = 0;
    uint64_t n      = 0;
    uint64_t l;
    uint64_t bases  = 0;
    uint64_t maxlen = LOAD_SIZE - B64;
    uint8_t  *seqbuffer = sqz->seqbuffer;
    uint8_t  *qltbuffer = sqz->qualbuffer;
    while ( kseq_read(seq) >= 0 ) {
        l = seq->seq.l;
        n++;
        if (!sqz_loadname(sqz, seq)) {
            offset = 0;
            goto exit;
        }
        //Determine if current sequence can be loaded completely
        if ( l + 1  > maxlen) {
            memcpy(seqbuffer + offset, &l, B64);
            offset += B64;
            memcpy(seqbuffer + offset, seq->seq.s, maxlen);
            memcpy(qltbuffer + offset, seq->qual.s, maxlen);
            offset += maxlen;
            //TODO this is not strictly necessary, but it brings me peaceof mind
            seqbuffer[offset] = 0;
            qltbuffer[offset] = 0;
            //From previous TODO
            //offset should be LOAD_SIZE + 1 after this increment
            offset++;
            bases += maxlen;
            sqz->endflag = 1;
            sqz->seqread = maxlen;
            sqz->prevlen = l;
            /*
              The sqzfastx_t struct must remember the sequence (and qualities)
              of this partially decoded record. Thus the contents of seq->seq.s
              (and seq->qual.s) must be transfered to the corresponding members
              of sqzfastx_t. Making shure those buffers are big enough.
            */
            sqz_rememberseq(sqz, seq, 1);
            //fprintf(stderr, "\tExiting from partial\n");
            goto exit;
        }
        bases += l;
        memcpy(seqbuffer + offset, &l, B64);
        offset += B64;
        memcpy(seqbuffer + offset, seq->seq.s, l + 1);
        memcpy(qltbuffer + offset, seq->qual.s, l + 1);
        offset += l + 1;
        if ( maxlen <= l + 1 + B64 ) break;
        maxlen -= l + 1 + B64;
    }
    //fprintf(stderr, "\tWe are donsies my friend\n");
    exit:
        sqz->n = n;
        sqz->bases = bases;
        sqz->offset = offset;
        return offset;
}


static uint64_t sqz_fastqeblock(sqzfastx_t *sqz)
{
    uint64_t l = sqz->prevlen;
    uint64_t seqleft = l - sqz->seqread;
    char *seq = sqz->pseq;
    char *qlt = sqz->pqlt;
    //buffer can be completely filled with current sequence
    if (seqleft >= LOAD_SIZE) {
        memcpy(sqz->seqbuffer, seq + sqz->seqread, LOAD_SIZE);
        memcpy(sqz->qualbuffer, qlt + sqz->seqread, LOAD_SIZE);
        sqz->seqread += LOAD_SIZE;
        sqz->offset   = LOAD_SIZE;
        sqz->bases   += LOAD_SIZE;
        return LOAD_SIZE;
    }
    //Rest of sequence can go into buffer
    memcpy(sqz->seqbuffer, seq + sqz->seqread, seqleft);
    memcpy(sqz->qualbuffer, qlt + sqz->seqread, seqleft);
    sqz->seqbuffer[seqleft] = 0;
    sqz->qualbuffer[seqleft] = 0;
    seqleft++;
    sqz->endflag = 0;
    sqz->offset = seqleft;
    sqz->bases += seqleft;
    return seqleft;
}


static uint64_t sqz_fastanblock(sqzfastx_t *sqz, kseq_t *seq)
{

    uint64_t offset = 0;
    uint64_t n      = 0;
    uint64_t l;
    uint64_t bases  = 0;
    uint64_t maxlen = LOAD_SIZE - B64;
    uint8_t *seqbuffer = sqz->seqbuffer;

    while ( (kseq_read(seq) >= 0) ) {
        l = seq->seq.l;
        n++;
        if (!sqz_loadname(sqz, seq)) {
            offset = 0;
            goto exit;
        }
        //Determine if current sequence can be loaded completely
        if (l + 1 > maxlen) {
            memcpy(seqbuffer + offset, &l, B64);
            offset += B64;
            memcpy(seqbuffer + offset, seq->seq.s, maxlen);

            offset += maxlen;

            seqbuffer[offset++] = 0;




            bases += maxlen;
            sqz->endflag = 1;
            sqz->seqread = maxlen;
            sqz->prevlen = l;
            /*
              The sqzfastx_t struct must remember the sequence (and qualities)
              of this partially decoded record. Thus the contents of seq->seq.s
              (and seq->qual.s) must be transfered to the corresponding members
              of sqzfastx_t. Making shure those buffers are big enough.
            */
            sqz_rememberseq(sqz, seq, 0);
            //fprintf(stderr, "\tExiting from partial\n");
            goto exit;
        }
        bases += l;
        memcpy(seqbuffer + offset, &l, B64); //Copy sequence length
        offset += B64;
        memcpy(seqbuffer + offset, seq->seq.s, l + 1); //Copy sequence + NULL
        offset += l + 1;
        if ( maxlen <= l + 1 + B64 ) break;
        maxlen -= l + 1 + B64;
    }
    exit:
        sqz->n = n;
        sqz->bases = bases;
        sqz->offset = offset;
        return offset;
}


static uint64_t sqz_fastaeblock(sqzfastx_t *sqz)
{
    uint64_t l = sqz->prevlen;
    uint64_t seqleft = l - sqz->seqread;
    char *seq = sqz->pseq;
    //buffer can be completely filled with current sequence
    if (seqleft >= LOAD_SIZE) {
        memcpy(sqz->seqbuffer, seq + sqz->seqread, LOAD_SIZE);
        sqz->seqread += LOAD_SIZE;
        sqz->offset   = LOAD_SIZE;
        sqz->bases   += LOAD_SIZE;
        return LOAD_SIZE;
    }
    //Rest of sequence can go into buffer
    memcpy(sqz->seqbuffer, seq + sqz->seqread, seqleft);
    sqz->seqbuffer[seqleft++] = '\0';
    sqz->endflag = 0;
    sqz->offset = seqleft;
    sqz->bases += seqleft;
    return seqleft;
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


void     sqz_sqzclose(sqz_File file)
{
    sqz_fastxkill(file.sqz);
    sqz_sqzblkkill(file.blk);
    fclose(file.fp);
}


char     sqz_readblksize(sqzblock_t *blk, FILE *fp)
{
    char ret = 0;
    uint64_t cmpsize;
    uint64_t dcpsize;
    uint64_t nelem;
    nelem =  fread(&dcpsize, B64, 1, fp);
    nelem += fread(&cmpsize, B64, 1, fp);
    if ( cmpsize > (blk->cmpsize) ) {
        blk->cmpbuff = realloc(blk->cmpbuff, cmpsize);
        if ( !(blk->cmpbuff) ) goto exit;
        blk->cmpsize = cmpsize;
    }
    if ( dcpsize > (blk->blksize) ) {
        blk->blkbuff = realloc(blk->blkbuff, dcpsize);
        if ( !(blk->blkbuff) ) goto exit;
        blk->blksize = dcpsize;
    }

    if ( (cmpsize != fread(blk->cmpbuff, 1, cmpsize, fp)) || (nelem != 2) )
        goto exit;
    blk->cmpsize = cmpsize;

    //cbytes = sqz_inflate(blk);
    if (dcpsize != sqz_inflate(blk))
        goto exit;
    blk->blksize = dcpsize;
    blk->newblk = 1;
    ret = 1;
    exit:
        return ret;
}


void     sqz_decode(sqzfastx_t *sqz, sqzblock_t *blk, uint64_t klibl)
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


static void sqz_rememberseq(sqzfastx_t *sqz, kseq_t *seq, char fqflag)
{
    if (sqz->plen < seq->seq.l) {
        sqz->plen = seq->seq.l + 1024;
        if (fqflag) sqz->pqlt = realloc(sqz->pqlt, sqz->plen);
        sqz->pseq = realloc(sqz->pseq, sqz->plen);
    }
    if (fqflag) strcpy(sqz->pqlt, seq->qual.s);
    strcpy(sqz->pseq, seq->seq.s);
}