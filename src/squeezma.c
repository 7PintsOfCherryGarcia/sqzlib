#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>

#include "squeezma.h"

uint64_t djb2(char *str)
{
    uint64_t hash = 5381;
    int c;
    while ( (c = *str++) )
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    return hash;
}


uint64_t bit2encode(const unsigned char *str, uint32_t strlen)
{
    uint64_t result = 0;
    for (uint32_t i = 0; i < strlen; i++) {
        result = (result << 2) | seq_nt4_table[str[i]];
    }
    return result;
}


unsigned char bit2decode(const uint64_t *mer, char *decoded, uint32_t len)
{
    unsigned char byte;
    unsigned char nbase = 0;
    uint64_t code = *mer;
    --len;
    do {
        byte = code & TWO_BIT_MASK;
        decoded[len] = seq_dec_table[byte];
        nbase += 1;
        code >>= 2;
    } while (len-- != 0);
    return nbase;
}


size_t sqz_zlibsqueez(void *seq, size_t seqlength, uint8_t *dest, size_t destlen, int level)
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
        fprintf(stderr, "Some other error\n");
        return (uint32_t)ret;
    }
    //Main compression loop
    strm.avail_in = seqlength;
    strm.next_in = seq;
    do {
        strm.avail_out = CHUNK;
        strm.next_out = out;
        ret = deflate(&strm, Z_FINISH);    /* no bad return value */
        if ( ret == Z_STREAM_ERROR ) {
            fprintf(stderr, "Some error\n");
            return (uint32_t)ret;  /* state not clobbered */
        }
        //check there is enough space in buffer
        have = CHUNK - strm.avail_out;
        if ( (wbytes + have) > destlen) {
            fprintf(stderr, "Compression error %lu %lu %lu\n",wbytes, have, destlen);
            ret = Z_ERRNO;
            break;
        }
        //Write to output buffer, keep track of written bytes
        memcpy(dest + wbytes, out, have);
        wbytes += have;
    } while (strm.avail_out == 0);
    //Check compression went well
    if (ret != Z_STREAM_END) wbytes = 0;
    deflateEnd(&strm);
    return wbytes;
}


size_t zlibpunch(void *cmpbuff, size_t cmplen, uint8_t *dest, size_t destlen)
{
    int ret;
    size_t wbytes = 0;
    unsigned have;
    z_stream strm;
    unsigned char out[CHUNK];
    //allocate inflate state
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK)
        return 0;
    strm.avail_in = cmplen;
    strm.next_in = cmpbuff;
    // run inflate() on input until output buffer not full
    do {
        strm.avail_out = CHUNK;
        strm.next_out = out;
        ret = inflate(&strm, Z_NO_FLUSH);
        if (ret == Z_STREAM_ERROR) {
            wbytes = 0;
            break;
        }  /* state not clobbered */
        switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                inflateEnd(&strm);
                return 0;
        }
        have = CHUNK - strm.avail_out;
        if ( (wbytes + have) > destlen) {
            fprintf(stderr, "Compression error %lu %u %lu\n",wbytes, have, destlen);
            ret = Z_ERRNO;
            break;
        }
        //Copy decompressed data into dest, keep tracked of amount of data decompressed
        memcpy(dest + wbytes, out, have);
        wbytes += have;
    } while (strm.avail_out == 0);
    /* clean up and return */
    inflateEnd(&strm);
    if ( ret != Z_STREAM_END ) {
        wbytes = 0;
    }
    return wbytes;
}


int sqz_cmpnflush(sqzfastx_t *sqz)
{
    fprintf(stderr, "Code buffer size: %lu\n", sqz->codeblk->offset);
    sqz->codeblk->offset = 0;
    return 1;
}


char sqz_getformat(const char *filename)
{
    char ret = 0;
    gzFile fp = gzopen(filename, "r");
    //ERROR
    if (!fp) return ret;
    kseq_t *seq = kseq_init(fp);
    //ERROR
    if (!seq) {
        gzclose(fp);
        return ret;
    }
    int l = kseq_read(seq);
    //ERROR
    if (l < 0) {
        fprintf(stderr, "[squeezma ERROR] sequence file format not recognized\n");
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


sqzfastx_t *sqz_fastxinit(const char *filename, size_t buffersize)
{
    sqzfastx_t *sqz = malloc(sizeof(sqzfastx_t));
    if (!sqz) {
        fprintf(stderr, "[squeezma ERROR] Memory error\n");
        return NULL;
    }
    sqz->n = 0;
    sqz->endflag = 0;
    sqz->toread = 0;
    sqz->prevlen = 0;
    sqz->codeblk = sqz_codeblkinit(LOAD_SIZE);
    if (!sqz->codeblk) {
        fprintf(stderr, "[squeezma ERROR] Memory error\n");
        free(sqz);
        return NULL;
    }
    //Get file format
    char fmt = sqz_getformat(filename);
    switch (fmt) {
        case 0:
            free(sqz);
            return NULL;
        //fasta
        case 1:
            fprintf(stderr, "[squeezma INFO] Detected FASTA format\n");
            sqz->fmt = fmt;
            //No need for quality buffer
            sqz->qualbuffer = NULL;
            sqz->seqbuffer = malloc(LOAD_SIZE + 1);
            if (!sqz->seqbuffer) {
                fprintf(stderr, "[squeezma ERROR] Failed to allocate sequence buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer[LOAD_SIZE] = 0;
            sqz->seqlen = 0;
            sqz->namebuffer = malloc(1*1024*1024);
            if (!sqz->namebuffer) {
                fprintf(stderr, "[squeezma ERROR] Failed to allocate name buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->namelen = 0;
            return sqz;
        case 2:
            fprintf(stderr, "[squeezma INFO] Detected FASTQ format\n");
            sqz->fmt = fmt;
            sqz->qualbuffer = malloc(LOAD_SIZE);
            if (!sqz->qualbuffer) {
                fprintf(stderr, "[squeezma ERROR] Failed to allocate quality buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer = malloc(LOAD_SIZE + 1);
            if (!sqz->seqbuffer) {
                fprintf(stderr, "[squeezma ERROR] Failed to allocate sequence buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->seqbuffer[LOAD_SIZE] = 0;
            sqz->seqlen = 0;
            sqz->namebuffer = malloc(1*1024*1024);
            if (!sqz->namebuffer) {
                fprintf(stderr, "[squeezma ERROR] Failed to allocate name buffer\n");
                free(sqz);
                return NULL;
            }
            sqz->namelen = 0;
            return sqz;
    }
    return NULL;
}


sqzcodeblock_t *sqz_codeblkinit(size_t size)
{
    sqzcodeblock_t *blk = malloc(sizeof(sqzcodeblock_t));
    if (!blk) return NULL;
    blk->codebuff = malloc(size);
    if (!blk->codebuff) {
        free(blk);
        return NULL;
    }
    blk->offset = 0;
    return blk;
}


void sqz_kill(sqzfastx_t *sqz)
{
    free(sqz->seqbuffer);
    free(sqz->namebuffer);
    free(sqz->codeblk->codebuff);
    free(sqz->codeblk);
    if (sqz->fmt == 2) free(sqz->qualbuffer);
    free(sqz);
}


int sqz_encodencompress(sqzfastx_t *sqz, size_t sqzsize)
{
    unsigned char *seq = NULL;
    unsigned char *qual = NULL;
    size_t *seqlen = NULL;
    size_t seqread;
    size_t lenbytes = sizeof(size_t);
    size_t k = 0;
    size_t n = 0;
    //Check if reading leftover sequence
    if (sqz->endflag) {
        seqlen = &(sqz->prevlen);
        seq = (unsigned char *)sqz->seqbuffer + k;
        qual = (unsigned char *)sqz->qualbuffer + k;
        //encode sequence
        size_t num = sqz_seqencode(seq, sqzsize, sqz->codeblk);
        //do something with the qualities
        //If this happens the buffer will only contain data from the sequence being read
        k = sqzsize;
        //Update how much sequence has been read
        sqz->toread += sqzsize;
        seqread = sqz->toread;
        sqz->endflag = 0;
    }
    //Process sequence data in buffer
    while ( k < sqzsize) {
        n++;
        //Extract length of sequences in buffer
        seqlen = (size_t *)(sqz->seqbuffer + k);
        k += lenbytes;
        seq = (unsigned char *)sqz->seqbuffer + k;
        qual = (unsigned char *)sqz->qualbuffer + k;
        //Determine how much sequence is contained in the buffer
        seqread = *seqlen < (sqzsize - k)?*seqlen + 1:sqzsize - k;
        //fprintf(stderr, "In sequence: %lu read %lu len %lu\n", n, seqread, *seqlen);
        k += seqread;
        //encode sequence
        size_t num = sqz_seqencode(seq, seqread, sqz->codeblk);
        //do something with the qualities
    }
    //Indicate if there is more sequence to read
    if (seqread < *seqlen) {
        //fprintf(stderr, "More seq please!! %lu %lu\n", seqread, *seqlen);
        sqz->endflag = 1;
        //Indicate how much sequence has been read
        sqz->toread = seqread;
        //Indicate length of sequence still needing loading
        sqz->prevlen = *seqlen;
    }
    if (k != sqzsize) {
        fprintf(stderr, "[ERROR] Unloading buffer\n");
        return 0;
    }
    return 1;
}


char sqz_loadfasta(sqzfastx_t *sqz, kseq_t *seq)
{
    //TODO
    /*
        In future pthread implementation, an sqzfastx array will be used. With number of
        elements equal to the number of threads being used.
    */
    size_t bleftover;
    size_t remaining;
    size_t offset = 0;
    size_t lenbytes = sizeof(seq->seq.l);
    size_t n = 0;
    int l;
    //Loop over sequences and load data into sqzfastx_t struct
    while ((l = kseq_read(seq)) >= 0) {
        n++;
        //process data until seqbuffer can't hold more data
        if (offset + seq->seq.l + 1 + lenbytes > LOAD_SIZE) {
            //Compute hom much buffer is available
            bleftover = LOAD_SIZE - offset;
            //Copy sequence length data
            memcpy(sqz->seqbuffer + offset, &(seq->seq.l), lenbytes);
            offset += lenbytes;
            bleftover -= lenbytes;
            //Copy as much seq data as we can fit in remaining buffer
            memcpy(sqz->seqbuffer + offset, seq->seq.s, bleftover);
            sqz->n = n;
            //Compress, buffer is guranteed to be full
            sqz_encodencompress(sqz, LOAD_SIZE);
            offset = 0;
            //Continue filling buffer with rest of sequence
            while ( bleftover != seq->seq.l + 1) {
                //Compute how much sequence remains
                remaining = seq->seq.l + 1 - bleftover;
                //buffer can be completely filled with current sequence
                if (remaining >= LOAD_SIZE) {
                    memcpy(sqz->seqbuffer + offset,
                           seq->seq.s + bleftover,
                           LOAD_SIZE);
                    sqz->n = 0;
                    //compress
                    sqz_encodencompress(sqz, LOAD_SIZE);
                    bleftover += LOAD_SIZE;
                    offset = 0;
                }
                //Rest of sequence can go into buffer
                else {
                    memcpy(sqz->seqbuffer + offset,
                           seq->seq.s + bleftover,
                           remaining);
                    bleftover += remaining;
                    offset += remaining;
                    sqz->n = 0;
                    sqz_encodencompress(sqz, remaining);
                    offset = 0;
                    //Indicate a sequence block must be closed
                }
            }
            //Current data block is finished. Buffers should be finilized and compressed
            fprintf(stderr, "Data ready for compression and flushing\n");
            sqz_cmpnflush(sqz);
            n = 0;
        }
        //copy sequence data into buffers
        else {
            //Copy sequence length data
            memcpy(sqz->seqbuffer + offset, &(seq->seq.l), lenbytes);
            offset += lenbytes;
            //Copy sequence content plus terminating null byte
            memcpy(sqz->seqbuffer + offset, seq->seq.s, seq->seq.l + 1);
            offset += seq->seq.l + 1;
        }
    }
    sqz->n = n;
    sqz_encodencompress(sqz, offset);
    return 1;
}


char sqz_loadfastq(sqzfastx_t *sqz, kseq_t *seq)
{
    //TODO
    /*
        In future pthread implementation, an sqzfastx array will be used. With number of
        elements equal to the number of threads being used.
    */
    size_t bleftover;
    size_t remaining;
    size_t offset = 0;
    size_t lenbytes = sizeof(seq->seq.l);
    size_t n = 0;
    int l;
    //Loop over sequences and load data into sqzfastx_t struct
    while ((l = kseq_read(seq)) >= 0) {
        n++;
        /*
          When a buffer is filled (Can't hold the entire kseq sequence + the lenth of
          the next sequence), the length of the current sequence is stored as well as
          any bases the buffer can accomodate.
        */
        if (offset + seq->seq.l + 1 + lenbytes + lenbytes > LOAD_SIZE) {
            //Compute how much buffer is available
            bleftover = LOAD_SIZE - offset;
            //Copy sequence length data
            fprintf(stderr, "Blockof size: %lu\n", n);
            memcpy(sqz->seqbuffer + offset, &(seq->seq.l), lenbytes);
            offset += lenbytes;
            bleftover -= lenbytes;
            //Compute how much sequence can be loaded to the buffer. There are two option:
            //as much sequence as leftover buffer space can be copied or the entire sequence
            remaining = bleftover < seq->seq.l?bleftover:seq->seq.l + 1;
            //Copy as much seq data as we can fit in remaining buffer
            memcpy(sqz->seqbuffer + offset, seq->seq.s, remaining);
            memcpy(sqz->qualbuffer + offset, seq->qual.s, remaining);
            offset += remaining;
            sqz->n = n;
            //Encode
            sqz_encodencompress(sqz, offset);
            offset = 0;
            //Keep track of how much sequence has been loaded
            bleftover = remaining;
            //Compute how much sequence is left to load
            remaining = seq->seq.l + 1 - remaining;
            //Continue filling buffer until there is no mre sequence to fill
            while ( remaining != 0) {
                //buffer can be completely filled with current sequence
                if (remaining >= LOAD_SIZE) {
                    memcpy(sqz->seqbuffer + offset,
                           seq->seq.s + bleftover,
                           LOAD_SIZE);
                    memcpy(sqz->qualbuffer + offset,
                           seq->qual.s + bleftover,
                           LOAD_SIZE);
                    sqz->n = 0;
                    //compress
                    sqz_encodencompress(sqz, LOAD_SIZE);
                    bleftover += LOAD_SIZE;
                    offset = 0;
                    remaining -= LOAD_SIZE;
                }
                //Rest of sequence can go into buffer
                else {
                    memcpy(sqz->seqbuffer + offset,
                           seq->seq.s + bleftover,
                           remaining);
                    memcpy(sqz->qualbuffer + offset,
                           seq->qual.s + bleftover,
                           remaining);
                    bleftover += remaining;
                    offset += remaining;
                    sqz->n = 0;
                    sqz_encodencompress(sqz, remaining);
                    offset = 0;
                    remaining = 0;
                }
            }
            //Current data block is finished. Buffers should be finilized and compressed
            fprintf(stderr, "Data ready for compression and flushing\n");
            if(!sqz_cmpnflush(sqz)) return 0;
            n = 0;
        }
        //copy sequence data into buffers
        else {
            //Copy sequence length data
            memcpy(sqz->seqbuffer + offset, &(seq->seq.l), lenbytes);
            offset += lenbytes;
            //Copy sequence string plus terminating null byte
            memcpy(sqz->seqbuffer + offset, seq->seq.s, seq->seq.l + 1);
            //Copy quality string plus terminating null byte
            memcpy(sqz->qualbuffer + offset, seq->qual.s, seq->seq.l + 1);
            offset += seq->seq.l + 1;
        }
    }
    if (n != 0) {
        sqz->n = n;
        sqz_encodencompress(sqz, offset);
        if(!sqz_cmpnflush(sqz)) return 0;
    }
    return 1;
}


int sqz_main(const char *filename)
{
    int ret = 0;
    //Initialize data structures
    gzFile fp = gzopen(filename, "r");
    if (!fp) return ret;
    kseq_t *seq = kseq_init(fp);
    if (!seq) {
        gzclose(fp);
        return ret;
    }
    sqzfastx_t *sqz = sqz_fastxinit(filename, LOAD_SIZE);
    if (!sqz) {
        kseq_destroy(seq);
        gzclose(fp);
        return ret;
    }

    switch (sqz->fmt) {
        case 1:
            if (!sqz_loadfasta(sqz, seq)) goto exit;
            break;
        case 2:
            if (!sqz_loadfastq(sqz, seq)) goto exit;
            break;
    }

    ret = 1;
    exit:
        sqz_kill(sqz);
        kseq_destroy(seq);
        gzclose(fp);
        return ret;
}


/*
  Returns length of block up to first null byte
*/
size_t sqz_seqlen(uint8_t *buffer, uint8_t *EOB)
{
    uint8_t *offset = buffer;
    while ( (offset != EOB) ) {
        if ( !(*offset) ) return (offset - buffer) + 1;
        offset++;
    }
    return offset - buffer;
}
