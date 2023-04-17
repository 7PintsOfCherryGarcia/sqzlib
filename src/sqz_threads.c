#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <pthread.h>

#include "../sqzlib/sqzlib.h"
#include "klib/kseq.h"
KSEQ_INIT(sqzFile, sqzread)


typedef struct {
    pthread_mutex_t mtx;
    pthread_cond_t  conscond;
    pthread_cond_t  readcond;
    pthread_cond_t  intraconscond;
    pthread_attr_t  thatt;
    int goread;
    int gocons;
    int wakethreadn;
    uint8_t doneq;
    uint8_t threadid;
    uint8_t nthread;
    const char *ifile;
    const char *ofile;
    sqzfastx_t **sqzqueue;
    FILE *ofp;
    uint8_t fqflag;
    uint8_t libfmt;
} sqzthread_t;


static void sqz_threadwait(int *flag,
                           pthread_cond_t *cond,
                           pthread_mutex_t *mtx)
{
    while (!(*flag)) {
        pthread_cond_wait(cond, mtx);
    }
}

static void sqz_syncwakeup(int *flag,
                           pthread_cond_t *cond,
                           pthread_mutex_t *mtx,
                           int n)
{
    if (n == *flag) {
        pthread_cond_broadcast(cond);
        return;
    }
    while (*flag != n)
        pthread_cond_wait(cond, mtx);
}

static int sqz_getthreadid(sqzthread_t *sqzthread)
{
    pthread_mutex_lock(&(sqzthread->mtx));
    int id = sqzthread->threadid++;
    sqz_threadwait(&(sqzthread->gocons),
                   &(sqzthread->conscond),
                   &(sqzthread->mtx));
    //Thread wakes up, now it has to wait for all other threads to wake up
    sqzthread->wakethreadn++;
    //Go to sleep again until all threads have woken up
    sqz_syncwakeup(&(sqzthread->wakethreadn),
                   &(sqzthread->intraconscond),
                   &(sqzthread->mtx),
                   sqzthread->nthread);
    pthread_mutex_unlock(&(sqzthread->mtx));
    return id;
}

static void sqz_wakereader(sqzthread_t *sqzthread)
{
    pthread_mutex_lock(&(sqzthread->mtx));
    sqzthread->doneq++;
    if (sqzthread->nthread == sqzthread->doneq) {
        sqzthread->doneq = 0;
        sqzthread->goread = 1;
        pthread_cond_signal(&(sqzthread->readcond));
    }
    sqzthread->gocons = 0;
    sqz_threadwait(&(sqzthread->gocons),
                   &(sqzthread->conscond),
                   &(sqzthread->mtx));
    //Thread wakes up and increases flag
    sqzthread->wakethreadn++;
    //Thread goes to sleep again until all threads have woken up
    sqz_syncwakeup(&(sqzthread->wakethreadn),
                   &(sqzthread->intraconscond),
                   &(sqzthread->mtx),
                   sqzthread->nthread);
    pthread_mutex_unlock(&(sqzthread->mtx));
}

static sqzfastx_t **sqz_sqzqueueinit(uint8_t n, uint8_t fmt)
{
    sqzfastx_t *sqz;
    sqzfastx_t **sqzqueue = calloc(n, sizeof(sqzfastx_t*));
    if (!sqzqueue) return NULL;
    uint8_t i, j;
    for (i = 0; i < n; i++) {
        sqz = sqzqueue[i];
        if ( !( sqz = sqz_fastxinit(fmt, 8LU*1024LU*1024LU)) ) {
            for(j = 0; j < i; j++)
                sqz_fastxkill(sqzqueue[j]);
            return NULL;
        }
        sqzqueue[i] = sqz;
    }
    return sqzqueue;
}

static void sqz_wakeconsumers(sqzthread_t *sqzthread)
{
    pthread_mutex_lock(&(sqzthread->mtx));
    sqzthread->goread = 0;
    sqzthread->gocons = 1;
    sqzthread->wakethreadn = 0;
    pthread_cond_broadcast(&(sqzthread->conscond));
    sqz_threadwait(&(sqzthread->goread),
                   &(sqzthread->readcond),
                   &(sqzthread->mtx));
    pthread_mutex_unlock(&(sqzthread->mtx));
}

static void sqz_threadkill(sqzthread_t *sqzthread)
{
    if (sqzthread) {
        if( sqzthread->sqzqueue )
            for (int i = 0; i < sqzthread->nthread; i++)
                if ( sqzthread->sqzqueue[i] )
                    sqz_fastxkill(sqzthread->sqzqueue[i]);
        free(sqzthread->sqzqueue);
        pthread_attr_destroy(&(sqzthread->thatt));
        pthread_mutex_destroy(&(sqzthread->mtx));
        pthread_cond_destroy(&(sqzthread->conscond));
        pthread_cond_destroy(&(sqzthread->readcond));
        pthread_cond_destroy(&(sqzthread->intraconscond));
        free(sqzthread);
    }
}

static void *sqz_dcmpthread(void *thread_data)
{
    sqzthread_t *sqzthread = (sqzthread_t *)thread_data;
    uint8_t n = sqzthread->nthread;

    uint8_t *outbuff = NULL;
    uint8_t *filebuff = NULL;
    uint64_t dsize = 0;

    pthread_mutex_lock(&(sqzthread->mtx));
    uint8_t id = sqzthread->threadid++;
    pthread_mutex_unlock(&(sqzthread->mtx));

    sqzFile sqzfp = sqzopen(sqzthread->ifile, "rb");
    if (!sqzfp) goto exit;

    sqzblock_t *blk =  sqz_sqzgetblk(sqzfp);

    outbuff = malloc(LOAD_SIZE);
    if (!outbuff) goto exit;

    filebuff = malloc(LOAD_SIZE);
    if (!filebuff) goto exit;
    uint64_t fbsize = LOAD_SIZE;
    uint64_t fbpos = 0;

    uint32_t nblk = sqz_getblocks(sqzfp);
    if (!nblk) goto exit;
    uint32_t currentblk = (uint32_t)id - 1;
    while (currentblk < nblk) {
        if ( sqz_loadblockn(sqzfp, currentblk) ) {
            fprintf(stderr, "[sqz ERROR]: Failed to load block %u\n", currentblk);
            goto exit;
        }
        do {
            dsize = sqz_decode(sqzfp, outbuff, LOAD_SIZE);
            fprintf(stderr, "\tdecoded size: %lu\n", dsize);
            if ( (fbpos + dsize) >= fbsize) {
                fbsize <<= 2;
                filebuff = realloc(filebuff, fbsize);
                if (!filebuff) goto exit;
            }
            memcpy(filebuff + fbpos, outbuff, dsize);
            fbpos += dsize;
        } while ( sqz_newblk(blk) );
        fprintf(stderr, "Done!!! about to dump!!\n");
        sleep(100);
        pthread_mutex_lock(&(sqzthread->mtx));
        fwrite(filebuff, fbpos, 1, sqzthread->ofp);
        fflush(sqzthread->ofp);
        pthread_mutex_unlock(&(sqzthread->mtx));
        fbpos = 0;
        currentblk += n;
    }
    exit:
        if(sqzfp) sqzclose(sqzfp);
        if(outbuff) free(outbuff);
        if (filebuff) free(filebuff);
        pthread_exit(NULL);
}

static void *sqz_decompressor(void *thrdata)
{
    sqzthread_t *sqzthread = (sqzthread_t *)thrdata;
    sqzthread->threadid++;
    sqzFile sqzfp = sqzopen(sqzthread->ifile, "rb");
    if (!sqzfp) goto exit;
    uint8_t fmt = sqz_format(sqzfp);
    if ( !(fmt & 4) ) {
        fprintf(stderr, "[sqz WARNING]: Not an sqz file\n");
        sqz_gzdump(sqzfp, sqzthread->ofile);
        sqzclose(sqzfp);
    }
    else {
        sqzclose(sqzfp);
        sqzthread->ofp = fopen(sqzthread->ofile, "wb");
        if (!sqzthread->ofp) goto exit;
        pthread_t *consumer_pool = malloc(sqzthread->nthread * sizeof(pthread_t));
        for (int i = 0; i < sqzthread->nthread; i++)
            if (pthread_create(consumer_pool + i,
                               &(sqzthread->thatt),
                               sqz_dcmpthread,
                               (void *)thrdata))
                goto exit;

        for (int i = 0; i < sqzthread->nthread; i++)
            if (pthread_join(consumer_pool[i], NULL))
                fprintf(stderr, "[sqz ERROR]: Thread error join\n");
        free(consumer_pool);
        fclose(sqzthread->ofp);
    }
    exit:
        return NULL;
}

static void *sqz_cmprthread(void *thread_data)
{
    sqzthread_t *sqzthread = thread_data;
    uint64_t cbytes = 0;
    int id = sqz_getthreadid(sqzthread);
    sqzfastx_t *sqz = sqzthread->sqzqueue[id - 1];
    uint8_t fqflag  = sqzthread->fqflag;
    int libfmt      = sqzthread->libfmt;
    //Do some work if there is available data
    while ( sqz_hasdata(sqz) ) {
        //Encode data
        sqz_fastXencode(sqz, fqflag);
        //Compress block
        cbytes = sqz_blkcompress(sqz, 9, libfmt);
        //Write compressed block to output file
        pthread_mutex_lock(&(sqzthread->mtx));
        sqz_blkdump(sqz, cbytes, sqzthread->ofp);
        fflush(sqzthread->ofp);
        sqz_resetsqz(sqz);
        pthread_mutex_unlock(&(sqzthread->mtx));
        if (sqz_readend(sqz))
            break;
        //Thread done, signal reader and go to sleep until new data arrives
        sqz_wakereader(sqzthread);
    }
    pthread_exit(NULL);
}

static uint32_t sqz_cmpreadloop(sqzthread_t *sqzthread, sqzFile sqzfp)
{
    uint32_t ret = 0, b = 0;
    sqzfastx_t **sqzqueue = sqzthread->sqzqueue;
    kseq_t *seq = kseq_init(sqzfp);
    uint8_t n = sqzthread->nthread, t = 0, i = 0, j = 0;
    uint8_t fqflag = sqz_isfq(sqzfp);
    //Start consumer pool
    pthread_t *consumer_pool = malloc(n * sizeof(pthread_t));
    if (!consumer_pool) goto exit;
    for (i = 0; i < n; i++)
        //BUG - If one thread fails to be created, the rest will hang
        if (pthread_create(consumer_pool + i,
                           &(sqzthread->thatt),
                           sqz_cmprthread,
                           (void *)sqzthread))
            goto exit;
    while (sqz_loadfastX(sqzqueue[t], fqflag, seq)) {
        b++;
        t++;
        if (n == t) {
            t = 0;
            sqz_wakeconsumers(sqzthread);
        }
    }
    if (t % n)
        for (i = 0; i < t; i++)
            sqz_setlastread(sqzqueue[i]);
    for (i = t; i < n; i++)
        sqz_setnodata(sqzqueue[i]);
    //Wake consumers one last time
    sqzthread->goread = 0;
    sqzthread->gocons = 1;
    sqzthread->wakethreadn = 0;
    pthread_cond_broadcast(&(sqzthread->conscond));
    ret = b;
    exit:
        //Wait until done
        if (seq) kseq_destroy(seq);
        for (j = 0; j < i; j++)
            if (pthread_join(consumer_pool[j], NULL))
                fprintf(stderr, "[sqz ERROR]: Thread error join\n");
        if (consumer_pool) free(consumer_pool);
        return ret;
}

static void *sqz_compressor(void *thrdata)
{
    sqzthread_t *sqzthread = (sqzthread_t *)thrdata;
    sqzthread->threadid++;
    sqzFile sqzfp = sqzopen(sqzthread->ifile, "r");
    if (!sqzfp) goto exit;
    uint8_t fmt = sqz_format(sqzfp);
    if (fmt & 4) {
        fprintf(stderr, "[sqz]: File already sqz encoded and compressed\n");
        goto exit;
    }

    uint8_t n = sqzthread->nthread;
    sqzthread->sqzqueue = sqz_sqzqueueinit(n, fmt);
    if ( !(sqzthread->sqzqueue) ) goto exit;
    sqzthread->fqflag = sqz_isfq(sqzfp);

    sqzthread->ofp = fopen(sqzthread->ofile, "wb");
    if (!sqzthread->ofp) goto exit;
    if ( !sqz_filehead(sqzthread->fqflag ? 2 : 1,
                       sqzthread->libfmt,
                       sqzthread->ofp) ) {
        fclose(sqzthread->ofp);
        return NULL;
    }
    fflush(sqzthread->ofp);

    uint32_t b = sqz_cmpreadloop(sqzthread, sqzfp);

    sqzfastx_t **sqzqueue = sqzthread->sqzqueue;
    //Log number of sequences and blocks
    uint64_t nseq = 0;
    for (int i = 0; i < n; i++)
        nseq += sqz_getn(sqzqueue[i]);
    sqz_filetail(nseq, b, sqzthread->ofp);
    exit:
        if (sqzfp) sqzclose(sqzfp);
        return NULL;
}

static sqzthread_t *sqz_threadinit(const char *i,
                                   const char *o,
                                   uint8_t lib,
                                   uint8_t n)
{
    sqzthread_t *sqzthread = calloc(1, sizeof(sqzthread_t));
    if (!sqzthread) return NULL;
    sqzthread->ifile   = i;
    sqzthread->ofile   = o;
    sqzthread->libfmt  = lib;
    sqzthread->nthread = n;
    pthread_mutex_init(&(sqzthread->mtx), NULL);
    pthread_cond_init(&(sqzthread->conscond), NULL);
    pthread_cond_init(&(sqzthread->readcond), NULL);
    pthread_cond_init(&(sqzthread->intraconscond), NULL);
    pthread_attr_init(&(sqzthread->thatt));
    return sqzthread;
}

uint8_t sqz_decompress(const char *i, const char *o, uint8_t n)
{
    uint8_t ret = 1;
    sqzthread_t *sqzthread = sqz_threadinit(i, o, 0, n);
    if (!sqzthread) goto exit;
    pthread_t r;
    if (pthread_create(&r,
                       &(sqzthread->thatt),
                       sqz_decompressor,
                       (void *)sqzthread)) {
        fprintf(stderr, "[sqz ERROR]: thread failed to launch.\n");
        goto exit;
    }
    if (pthread_join(r, NULL)) {
        fprintf(stderr, "\t[sqz ERROR]: thread join failed.\n");
        goto exit;
    }
    ret = 0;
    exit:
        sqz_threadkill(sqzthread);
        return ret;
}

uint8_t sqz_compress(const char *i, const char *o, uint8_t lib, uint8_t n)
{
    uint8_t ret = 1;
    sqzthread_t *sqzthread = sqz_threadinit(i, o, lib, n);
    if (!sqzthread) goto exit;
    pthread_t r;
    if (pthread_create(&r,
                       &(sqzthread->thatt),
                       sqz_compressor,
                       (void *)sqzthread)) {
        fprintf(stderr, "[sqz ERROR]: thread failed to launch.\n");
        goto exit;
    }
    if (pthread_join(r, NULL)) {
        fprintf(stderr, "\t[sqz ERROR]: thread join failed.\n");
        goto exit;
    }
    ret = 0;
    exit:
        sqz_threadkill(sqzthread);
        return ret;
}
