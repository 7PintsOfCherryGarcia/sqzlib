#include <stdio.h>
#include <zlib.h>
#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread)
#define KLIB
#include "sqz_kseq.h"


unsigned char sqz_getformat(const char *filename)
{
    unsigned char ret = 0;
    if ( (ret = sqz_checksqz(filename)) ) {
        fprintf(stderr, "DETECTED sqz\n");
        return ret;
    }
    gzFile fp = gzopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "[sqzlib ERROR]: Failed to open file: %s\n", filename);
        return ret;
    }

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


char sqz_kseqinit(sqzfastx_t *sqz)
{
    char ret = 0;
    sqz->fp = gzopen(sqz->filename, "r");
    //ERROR
    if (!sqz->fp) {
        fprintf(stderr,
                "[libsqz ERROR: zlib] Failed to open file: %s\n", sqz->filename);
        goto exit;
    }
    sqz->seq = kseq_init(sqz->fp);
    //ERROR
    if (!sqz->seq) {
        gzclose(sqz->fp);
        fprintf(stderr,
                "[libsqz ERROR: kseq] Failed to init kseq: %s\n", sqz->filename);
        goto exit;
    }
    ret = 1;
    exit:
        return ret;
}

/*
    Loads sqz's seqbuffer and qualbuffer with at most LOAD_SIZE bases/quality
    data. Returns number of bytes used to store the data
*/
size_t sqz_loadfastq(sqzfastx_t *sqz)
{
    if (!sqz->endflag) return sqz_newblock(sqz);
    return sqz_endblock(sqz);
}


size_t sqz_newblock(sqzfastx_t *sqz)
{
    size_t bleftover;
    size_t offset = 0;
    size_t lenbytes = sizeof(sqz->seq->seq.l);
    size_t n = 0;
    //Loop over sequences and load data into sqzfastx_t struct
    while ( kseq_read(sqz->seq) >= 0) {
        //TODO store name
        //fprintf(stderr, "NAME: %s LEN: %lu LEN: %lu\n",
        //        sqz->seq->name.s,
        //        sqz->seq->name.l,
        //        strlen(sqz->seq->name.s));
        sqz->bases += sqz->seq->seq.l;
        n++;
        /*
          When a buffer is filled (Can't hold the entire kseq sequence + the
          lenth of the next sequence), the length of the current sequence is
          stored as well as any bases the buffer can accomodate.
        */
        //TODO move to function
        if (offset + sqz->seq->seq.l + 1 + lenbytes + lenbytes > LOAD_SIZE) {
            //Compute how much buffer is available
            bleftover = LOAD_SIZE - offset;
            //Copy sequence length data
            memcpy(sqz->seqbuffer + offset, &(sqz->seq->seq.l), lenbytes);
            offset += lenbytes;
            bleftover -= lenbytes;
            /*
            Compute how much sequence can be loaded to the buffer.
            There are two option: as much sequence as leftover buffer, or the
            entire sequence
            TODO - Problem: If only a fraction of a sequence is loaded, no null
                   bytes is copied which messes up the buffer unloading (I think)
            */
            sqz->rem = bleftover < sqz->seq->seq.l ? bleftover: sqz->seq->seq.l;
            //Copy as much seq data as we can fit in remaining buffer
            memcpy(sqz->seqbuffer + offset, sqz->seq->seq.s, sqz->rem);
            memcpy(sqz->qualbuffer + offset, sqz->seq->qual.s, sqz->rem);
            offset += sqz->rem;
            //Add null byte after loading data into buffers
            sqz->seqbuffer[offset] = '\0';
            sqz->qualbuffer[offset] = '\0';
            offset++;
            sqz->n = n;
            sqz->rem = sqz->seq->seq.l - sqz->rem;

            if (sqz->rem != 0) {
                //Store length of sequence that could not complete loading
                sqz->prevlen = sqz->seq->seq.l;
                //Set flag to indicate there is more sequence to load
                sqz->endflag = 1;
            }
            sqz->offset = offset;
            return offset;
        }
        //copy sequence data into buffers
        else {
            //Copy sequence length data
            memcpy(sqz->seqbuffer + offset, &(sqz->seq->seq.l), lenbytes);
            offset += lenbytes;
            //Copy sequence string plus terminating null byte
            memcpy(sqz->seqbuffer+offset, sqz->seq->seq.s, sqz->seq->seq.l+1);
            //Copy quality string plus terminating null byte
            memcpy(sqz->qualbuffer+offset, sqz->seq->qual.s, sqz->seq->seq.l+1);
            offset += sqz->seq->seq.l + 1;
        }
    }

    if (n != 0)
        sqz->n = n;
    sqz->offset = offset;
    return offset;
}


size_t sqz_endblock(sqzfastx_t *sqz)
{
    size_t lsize;
    //buffer can be completely filled with current sequence
    if (sqz->rem >= LOAD_SIZE) {
        memcpy(sqz->seqbuffer,
               sqz->seq->seq.s + (sqz->seq->seq.l + 1 - sqz->rem),
               LOAD_SIZE);
        memcpy(sqz->qualbuffer,
               sqz->seq->qual.s + (sqz->seq->qual.l + 1 - sqz->rem),
               LOAD_SIZE);
        sqz->rem -= LOAD_SIZE;
        lsize = LOAD_SIZE;
        sqz->offset = LOAD_SIZE;
        return lsize;
    }

    //Rest of sequence can go into buffer
    memcpy(sqz->seqbuffer,
           sqz->seq->seq.s + sqz->toread,
           sqz->rem + 1);
    memcpy(sqz->qualbuffer,
           sqz->seq->qual.s + sqz->toread,
           sqz->rem + 1);
    lsize = sqz->rem + 1;
    sqz->endflag = 0;
    sqz->offset = lsize;
    return lsize;
}


size_t sqz_loadfasta(sqzfastx_t *sqz)
{
    if (!sqz->endflag) return sqz_fastanblock(sqz);
    return sqz_fastaeblock(sqz);
}


size_t sqz_fastanblock(sqzfastx_t *sqz)
{
    //TODO handle errors
    size_t bleftover;
    size_t offset = 0;
    //TODO change to constant
    size_t lenbytes = sizeof(sqz->seq->seq.l);
    size_t n = 0;
    sqz->namelen = 0;
    //Loop over sequences and load data into sqzfastx_t struct
    while ( kseq_read(sqz->seq) >= 0) {
        memcpy(sqz->namebuffer + sqz->namelen,
               sqz->seq->name.s,
               sqz->seq->name.l + 1);
        sqz->namelen += sqz->seq->name.l + 1;
        if (sqz->namelen + 100 >= sqz->maxname) {
            sqz->namebuffer = realloc(sqz->namebuffer, sqz->maxname*2);
            if (!(sqz->namebuffer)) {
                fprintf(stderr, "[sqzlib ERROR]: Memory error\n");
            }
            sqz->maxname *= 2;
            fprintf(stderr, "More mem for names: %lu\n", sqz->maxname);
        }
        sqz->bases += sqz->seq->seq.l;
        n++;
        /*
          When a buffer is filled (Can't hold the entire kseq sequence + the
          lenth of the next sequence), the length of the current sequence is
          stored as well as any bases the buffer can accomodate.
        */
        //TODO move to function
        if (offset + sqz->seq->seq.l + 1 + lenbytes + lenbytes > LOAD_SIZE) {
            //Compute how much buffer is available
            bleftover = LOAD_SIZE - offset;
            //Copy sequence length data
            memcpy(sqz->seqbuffer + offset, &(sqz->seq->seq.l), lenbytes);
            offset += lenbytes;
            bleftover -= lenbytes;
            /*
            Compute how much sequence can be loaded to the buffer.
            There are two option: as much sequence as leftover buffer, or the
            entire sequence
            TODO - Problem: If only a fraction of a sequence is loaded, no null
                   bytes is copied which messes up the buffer unloading (I think)
            */
            sqz->rem = bleftover < sqz->seq->seq.l ? bleftover: sqz->seq->seq.l;
            //Copy as much seq data as we can fit in remaining buffer
            memcpy(sqz->seqbuffer + offset, sqz->seq->seq.s, sqz->rem);
            offset += sqz->rem;
            //Add null byte after loading data into buffers
            sqz->seqbuffer[offset] = '\0';
            offset++;
            sqz->n = n;
            sqz->rem = sqz->seq->seq.l - sqz->rem;

            if (sqz->rem != 0) {
                //Store length of sequence that could not complete loading
                sqz->prevlen = sqz->seq->seq.l;
                //Set flag to indicate there is more sequence to load
                sqz->endflag = 1;
            }
            sqz->offset = offset;
            return offset;
        }
        //copy sequence data into buffers
        else {
            //Copy sequence length data
            memcpy(sqz->seqbuffer + offset, &(sqz->seq->seq.l), lenbytes);
            offset += lenbytes;
            //Copy sequence string plus terminating null byte
            memcpy(sqz->seqbuffer+offset, sqz->seq->seq.s, sqz->seq->seq.l+1);
            offset += sqz->seq->seq.l + 1;
        }
    }

    if (n != 0)
        sqz->n = n;
    sqz->offset = offset;
    return offset;
}


size_t sqz_fastaeblock(sqzfastx_t *sqz)
{
    size_t lsize;
    //buffer can be completely filled with current sequence
    if (sqz->rem >= LOAD_SIZE) {
        memcpy(sqz->seqbuffer,
               sqz->seq->seq.s + (sqz->seq->seq.l + 1 - sqz->rem),
               LOAD_SIZE);
        sqz->rem -= LOAD_SIZE;
        lsize = LOAD_SIZE;
        sqz->offset = LOAD_SIZE;
        return lsize;
    }

    //Rest of sequence can go into buffer
    memcpy(sqz->seqbuffer,
           sqz->seq->seq.s + sqz->toread,
           sqz->rem + 1);
    lsize = sqz->rem + 1;
    sqz->endflag = 0;
    sqz->offset = lsize;
    return lsize;
}


void sqz_kill(sqzfastx_t *sqz)
{
    if (sqz) {
        gzclose(sqz->fp);
        kseq_destroy(sqz->seq);
        free(sqz->seqbuffer);
        free(sqz->namebuffer);
        if ( (sqz->fmt == 2) | (sqz->fmt == 14) ) free(sqz->qualbuffer);
        free(sqz);
    }
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
    //Read format
    tmp += fread(&sqz, 1, 1, fp);
    fmt |= sqz;
    //Read compression library
    tmp += fread(&sqz, 1, 1, fp);
    fmt |= sqz << 3;
    fclose(fp);
    return fmt;
}
