#include <stdio.h>
#include <zlib.h>
#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread)
#define KLIB
#include "sqz_kseq.h"


unsigned char sqz_getformat(const char *filename)
{
    unsigned char ret = 0;
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


char sqz_kseqinit(sqzfastx_t *sqz)
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
    if (!sqz->endflag) return sqz_fastqnblock(sqz);
    return sqz_fastqeblock(sqz);
}


uint64_t sqz_fastqnblock(sqzfastx_t *sqz)
{
    uint64_t offset = 0;
    uint64_t n = 0;
    uint64_t l;
    while ( kseq_read(sqz->seq) >= 0 ) {
        l = sqz->seq->seq.l;
        sqz->bases += l;
        n++;
        if (!sqz_loadname(sqz, sqz->seq->name)) {
            offset = 0;
            goto exit;
        }
        if (offset + l + 1 + B64 + B64 > LOAD_SIZE) {
            sqz->n = n;
            offset = sqz_fastqwrap(sqz, offset);
            goto exit;
        }
        memcpy(sqz->seqbuffer + offset, &l, B64);
        offset += B64;
        memcpy(sqz->seqbuffer + offset, sqz->seq->seq.s, l + 1);
        memcpy(sqz->qualbuffer + offset, sqz->seq->qual.s, l + 1);
        offset += l + 1;
    }
    sqz->n = n;
    sqz->offset = offset;
    exit:
        return offset;
}


char sqz_loadname(sqzfastx_t *sqz, kstring_t name)
{
    char ret = 0;
    memcpy(sqz->namebuffer + sqz->namelen, name.s, name.l + 1);
    sqz->namelen += name.l + 1;
    if (sqz->namelen + 100 >= sqz->maxname) {
        sqz->namebuffer = realloc(sqz->namebuffer, sqz->maxname*2);
        if (!(sqz->namebuffer)) {
            fprintf(stderr, "[sqzlib ERROR]: Memory error\n");
            goto exit;
        }
        sqz->maxname *= 2;
        fprintf(stderr,
                "[sqzlib INFO]:Allocating memory %lu\n", sqz->maxname);
    }
    ret = 1;
    exit:
        return ret;
}


uint64_t sqz_fastqwrap(sqzfastx_t *sqz, uint64_t offset)
{
    //Compute how much buffer is available
    uint64_t bleftover = LOAD_SIZE - offset;
    uint64_t l = sqz->seq->seq.l;
    //Copy sequence length data
    memcpy(sqz->seqbuffer + offset, &l, B64);
    offset += B64;
    bleftover -= B64;
    /*
    Compute how much sequence can be loaded to the buffer.
    There are two option: as much sequence as leftover buffer, or the
    entire sequence
    TODO - Problem: If only a fraction of a sequence is loaded, no null
           bytes is copied which messes up the buffer unloading (I think)
    */
    sqz->rem = bleftover < l ? bleftover: l;
    //Copy as much seq data as we can fit in remaining buffer
    memcpy(sqz->seqbuffer + offset, sqz->seq->seq.s, sqz->rem);
    memcpy(sqz->qualbuffer + offset, sqz->seq->qual.s, sqz->rem);
    offset += sqz->rem;
    //Add null byte after loading data into buffers
    sqz->seqbuffer[offset] = '\0';
    sqz->qualbuffer[offset] = '\0';
    offset++;
    //sqz->n = n;
    sqz->rem = l - sqz->rem;
    if (sqz->rem != 0) {
        //Store length of sequence that could not complete loading
        sqz->prevlen = l;
        //Set flag to indicate there is more sequence to load
        sqz->endflag = 1;
    }
    sqz->offset = offset;
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


size_t sqz_loadfasta(sqzfastx_t *sqz)
{
    if (!sqz->endflag) return sqz_fastanblock(sqz);
    return sqz_fastaeblock(sqz);
}


uint64_t sqz_fastanblock(sqzfastx_t *sqz)
{
    uint64_t offset = 0;
    uint64_t n = 0;
    uint64_t l;
    while ( kseq_read(sqz->seq) >= 0) {
        l = sqz->seq->seq.l;
        sqz->bases += l;
        n++;
        if (!sqz_loadname(sqz, sqz->seq->name)) {
            offset = 0;
            goto exit;
        }
        if (offset + l + 1 + B64 + B64 > LOAD_SIZE) {
            sqz->n = n;
            offset = sqz_fastawrap(sqz, offset);
            goto exit;
        }
        memcpy(sqz->seqbuffer + offset, &(sqz->seq->seq.l), B64);
        offset += B64;
        memcpy(sqz->seqbuffer+offset, sqz->seq->seq.s, l + 1);
        offset += l + 1;
    }
    sqz->n = n;
    sqz->offset = offset;
    exit:
        return offset;
}


uint64_t sqz_fastawrap(sqzfastx_t *sqz, uint64_t offset)
{
    //Compute how much buffer is available
    uint64_t bleftover = LOAD_SIZE - offset;
    uint64_t l = sqz->seq->seq.l;
    //Copy sequence length data
    memcpy(sqz->seqbuffer + offset, &l, B64);
    offset += B64;
    bleftover -= B64;
    /*
    Compute how much sequence can be loaded to the buffer.
    There are two option: as much sequence as leftover buffer, or the
    entire sequence
    TODO - Problem: If only a fraction of a sequence is loaded, no null
           bytes is copied which messes up the buffer unloading (I think)
    */
    sqz->rem = bleftover < l ? bleftover: l;
    //Copy as much seq data as we can fit in remaining buffer
    memcpy(sqz->seqbuffer + offset, sqz->seq->seq.s, sqz->rem);
    offset += sqz->rem;
    //Add null byte after loading data into buffers
    sqz->seqbuffer[offset] = '\0';
    offset++;
    sqz->rem = l - sqz->rem;
    if (sqz->rem != 0) {
        //Store length of sequence that could not complete loading
        sqz->prevlen = l;
        //Set flag to indicate there is more sequence to load
        sqz->endflag = 1;
    }
    sqz->offset = offset;
    return offset;
}


uint64_t sqz_fastaeblock(sqzfastx_t *sqz)
{
    //buffer can be completely filled with current sequence
    if (sqz->rem >= LOAD_SIZE) {
        memcpy(sqz->seqbuffer,
               sqz->seq->seq.s + (sqz->seq->seq.l + 1 - sqz->rem),
               LOAD_SIZE);
        sqz->rem -= LOAD_SIZE;
        sqz->offset = LOAD_SIZE;
        return sqz->offset;
    }

    //Rest of sequence can go into buffer
    memcpy(sqz->seqbuffer,
           sqz->seq->seq.s + sqz->toread,
           sqz->rem + 1);
    sqz->offset = sqz->rem + 1;
    sqz->endflag = 0;
    return sqz->offset;
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
    fprintf(stderr, "fmt: %u\n", fmt);
    return fmt;
}
