uint64_t djb2(char *str)
{
    uint64_t hash = 5381;
    int c;
    while ( (c = *str++) )
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    return hash;
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
        sqz_seqencode(seq, sqzsize, sqz->codeblk);
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
        sqz_seqencode(seq, seqread, sqz->codeblk);
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









#define LOAD_SIZE 16*1024
#define CHUNK 16384
#define TWO_BIT_MASK (3)



//djb2 http://www.cse.yorku.ca/~oz/hash.html
uint64_t djb2(char *str);


uint64_t bit2encode(const unsigned char *str, uint32_t strlen);



/*
Returns number of bases decoded
*/
unsigned char bit2decode(const uint64_t *mer, char *decoded, uint32_t len);


//Comprsion functions shamelessly copypasata from https://zlib.net/zlib_how.html
/*
Compresses an array of data with zlibs deflate into the provided dest array
returns the number of compressed bytes
*/
size_t zlibsqueezze(void *nuts, size_t nutlength, uint8_t *dest, size_t destlen,int level);


/*
Decompress a zlib data stream stored in cmpbuff into dest
*/
size_t zlibpunch(void *cmpbuff, size_t cmplen, uint8_t *dest, size_t destlen);


/*
Returns pointer to first non ACGTacgt base in *strseq, If no nonACGTacgt bases are present
in the stirng, findn returns a pointer to the teminating null byte of strseq like in
strchrnul() in GNU
*/
const unsigned char *findn(const unsigned char *strseq);

/*
  //TODO update function documentation to reflect change of input arguments
    Sequence compression loop
    This function takes a sequence and encodes it via 2 bit enconding for bases
    AaCcGgTt and runlength encoding for N. The following format is used to store the
    encoded sequence:
        [slen|32bit|][blen|32bit]*[bblock|n*64bit](n*[nblock|n*8bit])
    slen   - sequence length
    blen   - block length: Sequence length up to an N bases
    bblock - 64 bit integers encoding up to 32 base kmers each for a total of blen bases
    nblock - 7 bit integers encoding up to 127 consecutive Ns. First bit indicates if
             current nblock is the last one before next blen,bblock or end of sequence

    Encoding algorithm:
        1.Read DNA sequence
        2.Write sequence length to [slen]
        While there is sequence to compress:
            3.Find first occurance of non ACGTacgt base
            If a non ACGTacgt base is found
                4.Compute sequence length up to found non ACGTacgt base
                5.Compute number of consecutive Ns
                6.Write sequence length up to found N to [blen]
                Loop over sequence (up to found N in chunks of 32 bases)
                    7.Encode 32mer or whatever is left of sequence into [bblock]
                Loop over number of Ns in chunks of 127
                    If no more Ns left
                        8.Encode 127 Ns or number of Ns left in [nblock]
                        9.Set bit 7 of currnt [nblock]
                    If more Ns need to be coded
                        8.Encode 127Ns as byte=127 in [nblock]
            10.Compute how much sequence there is left
            If there is sequence to be encoded (trailing bases after last occurance of N)
                11.Encode sequence as in 6
*/
size_t sqz_seqencode(const unsigned char *str, uint32_t strlen, sqzcodeblock_t *codeblk);


unsigned char writens(unsigned char numn, char *decoded);


size_t loopdecode(const uint8_t *buff);


int sqz_encodencompress(sqzfastx_t *sqz, size_t size);


char sqz_getformat(const char *filename);


sqzfastx_t *sqz_fastxinit(const char *filename, size_t buffersize);


void sqz_kill(sqzfastx_t *sqz);


char sqz_loadfasta(sqzfastx_t *sqz, kseq_t *seq);


char sqz_loadfastq(sqzfastx_t *sqz, kseq_t *seq);


sqzcodeblock_t *sqz_codeblkinit(size_t size);


int sqz_cmpnflush(sqzfastx_t *sqz);
