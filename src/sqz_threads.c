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
    int doneq;
    int wakethreadn;
    int threadid;
    int nthread;
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

static sqzfastx_t **sqz_sqzqueueinit(uint8_t nthread, uint8_t fmt)
{
    sqzfastx_t **sqzqueue = calloc(nthread, sizeof(sqzfastx_t*));
    if (!sqzqueue) return NULL;
    for (int i = 0; i < nthread; i++)
        if ( !(sqzqueue[i] = sqz_fastxinit(fmt, 64LU*1024LU*1024LU)) ) {
            for(int j = 0; j < i; j++)
                sqz_fastxkill(sqzqueue[j]);
            return NULL;
        }
    return sqzqueue;
}

static void *sqz_consumerthread(void *thread_data)
{
    sqzthread_t *sqzthread = thread_data;
    uint64_t cbytes = 0;
    sqzblock_t *blk = sqz_sqzblkinit(100UL*1024UL*1024UL);
    if (!blk) goto exit;
    int id = sqz_getthreadid(sqzthread);
    sqzfastx_t *sqz = sqzthread->sqzqueue[id - 1];
    uint8_t fqflag  = sqzthread->fqflag;
    int libfmt      = sqzthread->libfmt;
    //Do some work if there is available data
    while ( sqz_hasdata(sqz) ) {
        //Encode data
        sqz_fastXencode(sqz, blk, fqflag);
        //Compress block
        cbytes = sqz_blkcompress(blk, 9, libfmt);
        //Write compressed block to output file
        pthread_mutex_lock(&(sqzthread->mtx));
        sqz_blkdump(blk->cmpbuff, &(blk->blkpos), cbytes, sqzthread->ofp);
        fflush(sqzthread->ofp);
        blk->blkpos  = 0;
        sqz_newblk(sqz);
        pthread_mutex_unlock(&(sqzthread->mtx));
        if (sqz_endread(sqz))
            break;
        //Thread done, signal reader and go to sleep until new data arrives
        sqz_wakereader(sqzthread);
    }
    exit:
        sqz_blkkill(blk);
        pthread_exit(NULL);
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

static void *sqz_dcpthread(void *thread_data)
{
    sqzthread_t *sqzthread = (sqzthread_t *)thread_data;
    uint8_t nthread = sqzthread->nthread;

    uint8_t *outbuff = NULL;
    uint8_t *filebuff = NULL;
    uint64_t dsize = 0;

    pthread_mutex_lock(&(sqzthread->mtx));
    int id = sqzthread->threadid++;
    pthread_mutex_unlock(&(sqzthread->mtx));

    sqzFile sqzfp = sqzopen(sqzthread->ifile, "rb");
    if (!sqzfp) goto exit;
    uint8_t fqflag = sqz_sqzisfq(sqzfp);
    uint8_t libfmt = sqz_sqzgetcmplib(sqzfp);
    sqzblock_t *blk =  sqz_sqzgetblk(sqzfp);

    outbuff = malloc(LOAD_SIZE);
    if (!outbuff) goto exit;

    filebuff = malloc(LOAD_SIZE);
    if (!filebuff) goto exit;
    uint64_t fbsize = LOAD_SIZE;
    uint64_t fbpos = 0;

    uint64_t nblk = sqz_getblocks(sqzfp);
    if (!nblk) goto exit;

    uint64_t currentblk = (uint64_t)id;
    while (currentblk < nblk) {
        if (sqz_go2blockn(sqzfp, currentblk)) {
            fprintf(stderr, "[sqz ERROR]: Failed to read number of blocks.\n");
            goto exit;
        }
        if (!sqz_readblksize(blk, sqzfp, libfmt)) {
            fprintf(stderr, "[sqz ERROR]: Failed to read block %lu.\n",
                            currentblk);
            goto exit;
        }
        do {
            dsize = sqz_fastXdecode(blk, outbuff, LOAD_SIZE, fqflag);
            if ( (fbpos + dsize) >= fbsize) {
                fbsize *= 2;
                filebuff = realloc(filebuff, fbsize);
                if (!filebuff) goto exit;
            }
            memcpy(filebuff + fbpos, outbuff, dsize);
            fbpos += dsize;
        } while (blk->newblk);
        pthread_mutex_lock(&(sqzthread->mtx));
        fwrite(filebuff, fbpos, 1, sqzthread->ofp);
        fflush(sqzthread->ofp);
        pthread_mutex_unlock(&(sqzthread->mtx));
        fbpos = 0;
        currentblk += nthread;
    }
    exit:
        if(sqzfp) sqzclose(sqzfp);
        if(outbuff) free(outbuff);
        if (filebuff) free(filebuff);
        pthread_exit(NULL);
}

static void *sqz_decompressor(void *thread_data)
{
    sqzthread_t *sqzthread = (sqzthread_t *)thread_data;
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
                               sqz_dcpthread,
                               (void *)thread_data))
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

static uint8_t sqz_inflatefastX(sqzFile sqzfp,
                                FILE *ofp,
                                char fqflag,
                                uint8_t libfmt)
{
    uint8_t ret      = 1;
    uint8_t *outbuff = NULL;
    uint64_t dsize   = 0;
    uint64_t size    = 0;
    sqzblock_t *blk  = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;
    size = sqz_filesize(sqzfp);
    outbuff = malloc(LOAD_SIZE);
    if (!outbuff) goto exit;
    sqz_fseek(sqzfp, HEADLEN, SEEK_SET);
    while ( sqz_getfilepos(sqzfp) < size )
        {
            if (!sqz_readblksize(blk, sqzfp, libfmt)) goto exit;
            do {
                dsize = sqz_fastXdecode(blk, outbuff, LOAD_SIZE, fqflag);
                fwrite(outbuff, 1, dsize, ofp);
                fflush(ofp);
            } while (blk->newblk);
        }
    ret = 0;
    exit:
        sqz_blkkill(blk);
        free(outbuff);
        return ret;
}

static void *sqz_compressor(void *thrdata)
{
    sqzthread_t *sqzthread = (sqzthread_t *)thrdata;
    sqzthread->threadid++;

    sqzfastx_t **sqzqueue = NULL;
    kseq_t *seq = NULL;
    pthread_t *consumer_pool = NULL;

    sqzFile sqzfp = sqzopen(sqzthread->ifile, "r");
    if (!sqzfp) goto exit;
    uint8_t fmt = sqz_format(sqzfp);
    fprintf(stderr, "FMT: %u\n", fmt);
    uint8_t nthread = sqzthread->nthread;
    fprintf(stderr, "SLEEP TIME\n");
    sleep(100);
    sqzqueue = sqz_sqzqueueinit(nthread, fmt & 3);
    if (!sqzqueue) goto exit;
    sqzthread->sqzqueue = sqzqueue;

    uint8_t fqflag  = sqz_sqzisfq(sqzfp);
    sqzthread->fqflag = fqflag;
    //Initialize kseq object
    seq = kseq_init(sqzfp);
    if (!seq) goto exit;

    //Open output file and write sqz header
    sqzthread->ofp = fopen(sqzthread->ofile, "wb");
    if (!sqzthread->ofp) goto exit;
    if ( !sqz_filehead(sqzthread->fqflag ? 2 : 1,
                       sqzthread->libfmt,
                       sqzthread->ofp) ) {
        fclose(sqzthread->ofp);
        return NULL;
    }
    fflush(sqzthread->ofp);

    //Start consumer pool
    consumer_pool = malloc(nthread * sizeof(pthread_t));
    if (!consumer_pool) goto exit;
    for (int i = 0; i < nthread; i++)
        //TODO indicate error
        if (pthread_create(consumer_pool + i,
                           &(sqzthread->thatt),
                           sqz_consumerthread,
                           (void *)thrdata))
            goto exit;
    uint8_t  t = 0; //Count number of threads
    uint64_t b = 0; //Count number of blocks
    while (sqz_loadfastX(sqzqueue[t], fqflag, seq)) {
        b++;
        t++;
        if (nthread == t) {
            t = 0;
            sqz_wakeconsumers(sqzthread);
        }
    }
    if (t % nthread)
        for (int i = 0; i < t; i++)
            sqz_setlastread(sqzqueue[i]);
    for (int i = t; i < nthread; i++)
        sqz_endthread(sqzqueue[i]);
    //Wake consumers one last time
    sqzthread->goread = 0;
    sqzthread->gocons = 1;
    sqzthread->wakethreadn = 0;
    pthread_cond_broadcast(&(sqzthread->conscond));
    //Wait until done
    for (int i = 0; i < nthread; i++)
        if (pthread_join(consumer_pool[i], NULL))
            fprintf(stderr, "[sqz ERROR]: Thread error join\n");
    //Log number of sequences and blocks
    uint64_t n = 0;
    for (int i = 0; i < nthread; i++)
        n += sqz_getn(sqzqueue[i]);
    sqz_filetail(n, b, sqzthread->ofp);
    exit:
        if (consumer_pool) free(consumer_pool);
        if (seq) kseq_destroy(seq);
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

uint8_t sqz_decompress(const char *ifile,
                             const char *ofile,
                             uint8_t nthread)
{
    uint8_t ret = 1;
    sqzthread_t *sqzthread = sqz_threadinit(ifile,
                                            ofile,
                                            0,
                                            nthread);
    if (!sqzthread) goto exit;
    pthread_t rthread;
    if (pthread_create(&rthread,
                       &(sqzthread->thatt),
                       sqz_decompressor,
                       (void *)sqzthread)) {
        fprintf(stderr, "[sqz ERROR]: thread launch error.\n");
        goto exit;
    }
    if (pthread_join(rthread, NULL)) {
        fprintf(stderr, "\t[sqz ERROR]: thread join failed.\n");
        goto exit;
    }
    ret = 0;
    exit:
        sqz_threadkill(sqzthread);
        return ret;
}

uint8_t sqz_compress(const char *i, const char *o, uint8_t lib, uint8_t nthread)
{
    uint8_t ret = 1;
    sqzthread_t *sqzthread = sqz_threadinit(i, o, lib, nthread);
    if (!sqzthread) goto exit;
    pthread_t r;
    if (pthread_create(&r,
                       &(sqzthread->thatt),
                       sqz_compressor,
                       (void *)sqzthread)) {
        fprintf(stderr, "[sqz ERROR]: thread launch error.\n");
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
