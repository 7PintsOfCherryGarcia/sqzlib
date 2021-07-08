#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <pthread.h>

#include "../sqz_data.h"
#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread)

//TODO Move to header file
sqzblock_t *sqz_sqzblkinit(uint64_t size);
void sqz_sqzblkkill(sqzblock_t *blk);
sqzfastx_t *sqz_fastxinit(uint8_t fmt, uint64_t bsize);
void sqz_fastxkill(sqzfastx_t *sqz);
uint64_t sqz_loadfastX(sqzfastx_t *sqz, uint8_t fqflag, kseq_t *seq);
char sqz_fastXencode(sqzfastx_t *sqz, sqzblock_t *blk, uint8_t fqflag);
size_t sqz_deflate(sqzblock_t *blk, int level);
char sqz_zlibcmpdump(sqzblock_t *blk, uint64_t size, FILE *ofp);
void sqz_blkdestroy(sqzblock_t *blk);


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
    const char *filename;
    sqzfastx_t **sqzqueue;
    char fqflag;
    FILE *ofp;
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
        //fprintf(stderr, "broadcasting!!!\n");
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
    //fprintf(stderr, "Thread %d wokeup\n", id);
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


static void *sqz_consumerthread(void *thread_data)
{
    sqzthread_t *sqzthread = thread_data;
    uint64_t cbytes = 0;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;
    int id = sqz_getthreadid(sqzthread);
    sqzfastx_t *sqz = sqzthread->sqzqueue[id - 1];
    char fqflag = sqzthread->fqflag;
    //Do some work if there is available data
    while (sqz->endthread & 128) { //bit 7 is on if there is still data
        //Encode data
        sqz_fastXencode(sqz, blk, fqflag);
        //Check if there is leftover sequence that needs loading
        while (sqz->endflag) {
            sqz_loadfastX(sqz, fqflag, NULL);
            sqz_fastXencode(sqz, blk, fqflag);
        }
        //Compress block
        cbytes = sqz_deflate(blk, 9);

        //Write compressed block to output file
        pthread_mutex_lock(&(sqzthread->mtx));
        sqz_zlibcmpdump(blk, cbytes, sqzthread->ofp);
        fflush(sqzthread->ofp);
        blk->blkpos  = 0;
        sqz->namepos = 0;
        sqz->endflag = 0;
        blk->newblk  = 1;
        pthread_mutex_unlock(&(sqzthread->mtx));
        //bit 1 is on if thread had data, but reader has finished
        if (sqz->endthread & 1) break;
        //Thread done, signal reader and go to sleep until new data arrives
        sqz_wakereader(sqzthread);
    }
    exit:
        sqz_sqzblkkill(blk);
        return NULL;
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


static void *sqz_readerthread(void *thread_data)
{
    //TODO indicate errors
    sqzthread_t *sqzthread = thread_data;
    sqzfastx_t **sqzqueue   = sqzthread->sqzqueue;
    sqzthread->threadid++;
    int nthread = sqzthread->nthread;
    char fqflag = sqzthread->fqflag;
    //Initialize kseq object
    gzFile fp = gzopen(sqzthread->filename, "r");
    if (!fp) return NULL;
    kseq_t *seq = kseq_init(fp);
    if (!seq) return NULL;
    //Start consumer pool
    pthread_t *consumer_pool = malloc(nthread * sizeof(pthread_t));
    if (!consumer_pool) goto exit;
    for (int i = 0; i < nthread; i++)
        //TODO indicate error
        if (pthread_create(consumer_pool + i,
                           &(sqzthread->thatt),
                           sqz_consumerthread,
                           thread_data))
            goto exit;
    int thcounter = 0;
    /*
    While there is data to read, read data in blocks accorind to number
    of threads
    */
    while (sqz_loadfastX(sqzqueue[thcounter], fqflag, seq)) {
        thcounter++;
        if (nthread == thcounter) {
            thcounter = 0;
            sqz_wakeconsumers(sqzthread);
        }
    }

    if (thcounter % nthread) {
        //Set bit 1 of threads that got some data
        for (int i = 0; i < thcounter; i++)
            sqzqueue[i]->endthread |= 1;
    }

    //Unset bit 7 from rest of threads
    for (int i = thcounter; i < nthread; i++)
        sqzqueue[i]->endthread &= 127;
    //Wake consumers one last time
    sqzthread->goread = 0;
    sqzthread->gocons = 1;
    sqzthread->wakethreadn = 0;
    pthread_cond_broadcast(&(sqzthread->conscond));
    //Wait until done
    for (int i = 0; i < nthread; i++)
        if (pthread_join(consumer_pool[i], NULL))
            fprintf(stderr, "Thread error join\n");
    exit:
        free(consumer_pool);
        kseq_destroy(seq);
        gzclose(fp);
        return NULL;
}


static sqzthread_t *sqz_threadinit(FILE *ofp,
                                   const char *filename,
                                   char fqflag,
                                   int nthread,
                                   uint8_t fmt)
{
    sqzthread_t *sqzthread = calloc(1, sizeof(sqzthread_t));
    if (!sqzthread) return NULL;
    sqzthread->nthread = nthread;
    sqzthread->goread = 0;
    sqzthread->gocons = 0;
    sqzthread->doneq  = 0;
    sqzthread->wakethreadn = 0;
    sqzthread->threadid    = 0;
    sqzthread->fqflag      = fqflag;
    sqzthread->filename    = filename;
    sqzthread->ofp         = ofp;
    sqzthread->sqzqueue    = calloc(nthread, sizeof(sqzfastx_t*));
    if (!(sqzthread->sqzqueue)) {
        free(sqzthread);
        return NULL;
    }
    for (int i = 0; i < nthread; i++)
        if (!(sqzthread->sqzqueue[i] = sqz_fastxinit(fmt, LOAD_SIZE))) {
            //TODO kill succesfully initialized sqzfastx_t objects
            free(sqzthread);
            return NULL;
        }
    pthread_mutex_init(&(sqzthread->mtx), NULL);
    pthread_cond_init(&(sqzthread->conscond), NULL);
    pthread_cond_init(&(sqzthread->readcond), NULL);
    pthread_cond_init(&(sqzthread->intraconscond), NULL);

    pthread_attr_init(&(sqzthread->thatt));
    pthread_attr_setdetachstate(&(sqzthread->thatt), PTHREAD_CREATE_JOINABLE);
    return sqzthread;
}


static void sqz_threadkill(sqzthread_t *sqzthread)
{
    if (sqzthread) {
        for (int i = 0; i < sqzthread->nthread; i++)
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


uint8_t sqz_threadlauncher(FILE *ofp,
                           const char *filename,
                           char fqflag,
                           int nthread,
                           uint8_t fmt)
{
    uint8_t ret = 1;
    //Start thread object
    sqzthread_t *sqzthread = sqz_threadinit(ofp, filename, fqflag, nthread, fmt);
    //Launch reader thread
    pthread_t rthread;
    if (pthread_create(&rthread,
                       &(sqzthread->thatt),
                       sqz_readerthread,
                       (void *)sqzthread)) {
        fprintf(stderr, "[sqzlib]: Error Thread launch error.\n");
        goto exit;
    }
    //Wait for reader thred to finish
    if (pthread_join(rthread, NULL)) {
        fprintf(stderr, "\t[sqzlib]: Error, thread join failed.\n");
        goto exit;
    }
    ret = 0;
    exit:
        //Clean up
        sqz_threadkill(sqzthread);
        return ret;
}


static void sqz_singleread(kseq_t *seq)
{
    kseq_read(seq);
}

