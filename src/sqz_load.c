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
