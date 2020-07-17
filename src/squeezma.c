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


const unsigned char *findn(const unsigned char *strseq)
{
    //TODO bugyMcbuggerson, this function reads beyond allocated space if no non ACGTacgt
    //base is ever found. Solution: last byte of *strseq must be set to 0
    do {
    }
    while(seq_nt4_table[*(++strseq)] < 4);
    return strseq;
}


int sqz_cmpnflush(sqzfastx_t *sqz)
{
    sqz->codeblk->offset = 0;
    return 1;
}


//TODO control for cmpbuff length
size_t sqz_loopencode(const unsigned char *str, uint32_t strlen, sqzcodeblock_t *codeblk)
{
    //Set coding buffer to the next unused position in the block
    uint8_t *cmpbuff = codeblk->codebuff + codeblk->offset;
    const unsigned char *lstop = str;         //Track position within sequence
	  const unsigned char *nptr = str + strlen; //End of string
	  uint32_t nn;                              //Number of Ns
    unsigned char wn;                         //127 N block. 1 bit flag 7 bit count
	  const unsigned char *npos;                //Track positions where N occurs
	  uint32_t blocklen = 0;                    //Length of segment before first N
    uint64_t code;                            //2 bit encoded sequence
    size_t wbytes = 0;                        //Number of bytes written
    //Write sequence length
    memcpy(cmpbuff + wbytes, &strlen, sizeof(uint32_t));
    wbytes += sizeof(uint32_t);
    //Search for stretches of Ns
    do {
        //Get position of first non ACGTacgt base
	      npos = findn(lstop);
	      if (*npos) {
            //Determine block length up to found N
            blocklen = npos - lstop;
            //Determine number of consecutive Ns until next base
            nn = 0;
	          while ( seq_nt4_table[*npos] == 4) {
	              nn++;
		            npos++;
	          }
            //Write block length [blen]
            memcpy(cmpbuff + wbytes, &blocklen, sizeof(uint32_t));
            wbytes += sizeof(uint32_t);
            //Loop over sequence block in chunks of 32 bases
            for (uint32_t i = 0; i < blocklen; i+=32) {
                //Encode 32mer or what is left of sequence block
                code = bit2encode(lstop + i, ((i+32) <= blocklen)?32:blocklen-i);
                memcpy(cmpbuff + wbytes, &code, sizeof(uint64_t));
                wbytes += sizeof(uint64_t);
            }
            //Loop over Ns in chunks of 127
            for (uint32_t i = 0; i < nn; i +=127) {
                //Encode 127 Ns or what is left
                wn = ((i+127) <= nn)?127:nn-i;
                //Set first bit if last N block before next base
                if (i+127 >= nn) wn = wn | NEND;
                memcpy(cmpbuff + wbytes, &wn, 1);
                wbytes++;
            }
            //Trace next base after sequence of Ns
            lstop = npos;
	      }
    } while (*npos);
    //Detect and encode trailing bases
    blocklen = nptr - lstop;
    if (blocklen) {
        memcpy(cmpbuff + wbytes, &blocklen, sizeof(uint32_t));
        wbytes += sizeof(uint32_t);
        for (uint32_t i = 0; i < blocklen; i+=32) {
            code = bit2encode(lstop + i, ((i+32) <= blocklen)?32:blocklen-i);
            memcpy(cmpbuff + wbytes, &code, sizeof(uint64_t));
            wbytes += sizeof(uint64_t);
        }
    }
    //Move offset by number of bytes written
    codeblk->offset += wbytes;
    return wbytes;
}


unsigned char writens(unsigned char numn, char *decoded)
{
    unsigned char nwritten = 0;
    while (numn-- > 0) {
        decoded[numn] = 'N';
        nwritten++;
    }
    return nwritten;
}


size_t loopdecode(const uint8_t *buff)
{
    uint32_t strlen;
    char *seqpos;
    uint32_t blocklen;
    const uint64_t *mer;
    unsigned char merlen;
    uint32_t blockidx;
    const unsigned char *nstr;
    size_t bytes = 0;
    uint32_t decoded = 0;
    char *seqstr;
    //Get sequence length of current block
    memcpy(&strlen, buff, sizeof(uint32_t));
    bytes += sizeof(uint32_t);
    //Allocate memory for current sequence
    seqstr = malloc(strlen + 1);
    if (!seqstr) {
        fprintf(stderr, "MEM ERROR\n");
    }
    uint32_t strlen2 = strlen;
    seqpos = seqstr;
    //Decode for as long as there are bases to decode
    while (strlen > 0) {
        //Get block length
        memcpy(&blocklen, buff + bytes, sizeof(uint32_t));
        bytes += sizeof(uint32_t);
        //Get starting porition of 64bit integers
        mer = (uint64_t*)(buff + bytes);
        //Track block position
        blockidx = 0;
        //Loop over blocks
        //fprintf(stderr, "DECODING block of length: %u\n", blocklen);
        for (uint32_t i = 0; i < blocklen; i+=32) {
            //Blocks code 32mers except last block which may code a shorter kmer
            merlen = ((i+32) <= blocklen)?32:blocklen - i;
            decoded += bit2decode(mer + blockidx, seqpos, merlen);
            //Keep traked of amount of decoded data: 64bit per mer
            bytes += sizeof(uint64_t);
            //Move sequence pointer by abount of bases decoded (32 except last mer)
            seqpos += merlen;
            //Track next block to decompress
            blockidx++;
        }
        //TODO make sure decoded bases equals block length
        //Update amount of data left to decompress
        strlen -= blocklen;
        //strlen is != 0 when blocks are followed by an Nblock
        if (strlen) {
            //Read nbyte
            nstr = buff + bytes;
            bytes++;
            while (1) {
                unsigned char numn = *nstr & ~(1<<7);
                //fprintf(stderr, "N value of: %u\n", numn);
                decoded += writens(*nstr & ~(1<<7), seqpos);
                seqpos += numn;
                strlen-= numn;
                if (*nstr & 128)
                    break;
                nstr++;
                bytes++;
            }
        }
    }
    fprintf(stderr, "djb2 of seq: %lu\n", djb2(seqstr));

    FILE *raw = fopen("decomseq", "wb");
    fwrite(seqstr, 1, strlen2+1, raw);
    fclose(raw);


    free(seqstr);
    return strlen;
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
            sqz->quallen = 0;
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
            sqz->quallen = 0;
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
    blk->codesize = size;
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
    size_t *seqlen = NULL;
    size_t seqread;
    size_t lenbytes = sizeof(size_t);
    size_t k = 0;
    size_t n = 0;
    if (sqz->endflag) {
        seqlen = &(sqz->prevlen);
        seq = (unsigned char *)sqz->seqbuffer + k;
        //encode sequence
        fprintf(stderr, "||Compressing size: %lu\n", sqzsize);
        size_t num = sqz_loopencode(seq, sqzsize, sqz->codeblk);
        fprintf(stderr, "%lu bytes written ratio: %lu\n", num, sqz->codeblk->offset);

        //If this happens the buffer will only contain data from the sequence being read
        k = sqzsize;
        //Update how much sequence has been read
        sqz->toread += sqzsize;
        seqread = sqz->toread;
        sqz->endflag = 0;
    }
    //Extract length of sequences in buffer
    while ( k < sqzsize) {
        n++;
        seqlen = (size_t *)(sqz->seqbuffer + k);
        k += lenbytes;
        seq = (unsigned char *)sqz->seqbuffer + k;
        seqread = *seqlen < sqzsize - k?*seqlen + 1:sqzsize - k;
        k += seqread;
        //encode sequence
        fprintf(stderr, "Compressing size: %lu\n", *seqlen);
        size_t num = sqz_loopencode(seq, *seqlen, sqz->codeblk);
        fprintf(stderr, "%lu bytes written ratio: %lu\n", num, sqz->codeblk->offset);
    }

    if (seqread < *seqlen) {
        //Indicate there is more sequence to read
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
            memcpy(sqz->qualbuffer + offset, seq->qual.s, bleftover);
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
                    memcpy(sqz->qualbuffer + offset,
                           seq->qual.s + bleftover,
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
                    memcpy(sqz->qualbuffer + offset,
                           seq->qual.s + bleftover,
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
            if(!sqz_cmpnflush(sqz)) return 0;
            n = 0;
        }
        //copy sequence data into buffers
        else {
            //Copy sequence length data
            memcpy(sqz->seqbuffer + offset, &(seq->seq.l), lenbytes);
            offset += lenbytes;
            //Copy sequence content plus terminating null byte
            memcpy(sqz->seqbuffer + offset, seq->seq.s, seq->seq.l + 1);
            memcpy(sqz->qualbuffer + offset, seq->qual.s, seq->seq.l + 1);
            offset += seq->seq.l + 1;
        }
    }
    fprintf(stderr, "On exit\n");
    sqz->n = n;
    sqz_encodencompress(sqz, offset);
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
