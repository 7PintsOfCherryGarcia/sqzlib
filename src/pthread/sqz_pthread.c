#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
//#define SQZLIB
//#define KLIB
#include "../sqz_data.h"

sqzblock_t *sqz_sqzblkinit(uint64_t size);
char sqz_fastxinit(sqzfastx_t *sqz, unsigned char fmt, uint64_t bsize);


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
    sqzfastx_t *sqzqueue;
    unsigned char fmt;
} sqzthread_t;


void sqz_threadwait(int *flag, pthread_cond_t *cond, pthread_mutex_t *mtx)
{
    while (!(*flag)) {
        pthread_cond_wait(cond, mtx);
    }
}


void sqz_syncwakeup(int *flag, pthread_cond_t *cond, pthread_mutex_t *mtx, int n)
{
    if (n == *flag) {
        fprintf(stderr, "broadcasting!!!\n");
        pthread_cond_broadcast(cond);
        return;
    }
    while (*flag != n)
        pthread_cond_wait(cond, mtx);
}


int sqz_getthreadid(sqzthread_t *sqzthread)
{
    pthread_mutex_lock(&(sqzthread->mtx));
    int id = sqzthread->threadid++;
    sqz_threadwait(&(sqzthread->gocons),
                   &(sqzthread->conscond),
                   &(sqzthread->mtx));
    //Thread wakes up, now it has to wait for all other threadsto wake up
    sqzthread->wakethreadn++;
    fprintf(stderr, "Thread %d wokeup\n", id);
    //Go to sleep again until all threads have woken up
    sqz_syncwakeup(&(sqzthread->wakethreadn),
                   &(sqzthread->intraconscond),
                   &(sqzthread->mtx),
                   sqzthread->nthread);
    pthread_mutex_unlock(&(sqzthread->mtx));
    return id;
}


void *sqz_consumerthread(void *thread_data)
{
    sqzthread_t *sqzthread = thread_data;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;
    int id = sqz_getthreadid(sqzthread);
    fprintf(stderr, "Thread with ID %d\n", id);
    sqzfastx_t sqz = sqzthread->sqzqueue[id-1];
    fprintf(stderr, "Initializing sqz structs\n");
    fprintf(stderr, "\ttest: %d\n", sqz.endflag);
    //Do some work
    while (1)
        sleep(10);
    exit:
        return NULL;
}


void sqz_wakeconsumers(sqzthread_t *sqzthread)
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


void *sqz_readerthread(void *thread_data)
{
    sqzthread_t *sqzthread = thread_data;
    int id = sqzthread->threadid++;
    int nthread = sqzthread->nthread;
    sqzfastx_t *sqzqueue = sqzthread->sqzqueue;
    fprintf(stderr, "This is printed from a thread. %d\n", id);
    fprintf(stderr, "\tWe have %d threads\n", nthread);
    //Start consumer pool
    //TODO indicate error
    pthread_t *consumer_pool = malloc(nthread * sizeof(pthread_t));
    if (!consumer_pool) goto exit;
    for (int i = 0; i < nthread; i++)
        //TODO indicate error
        if (pthread_create(consumer_pool + i,
                           &(sqzthread->thatt),
                           sqz_consumerthread,
                           thread_data))
            goto exit;

    //Initialize kseq object
    gzFile fp = gzopen(sqzthread->filename, "r");
    kseq_t *seq = kseq_init();
    int thcounter = 0;
    //While there is data to read
    while (1) {
        thcounter++;
        if (nthread == thcounter) {
            thcounter = 0;
            sqz_wakeconsumers(sqzthread);
        }
        //Read data in blocks accorind to number of threads
    }
    //int counter = 0;
    //while (1) {
    //    counter++;
    //    fprintf(stderr, "counter: %d\n", counter);
    //    sleep(1);
    //    if (10 == counter) {
    //        counter = 0;
    //        fprintf(stderr, "Wake up consumers!!!\n");
    //        sqz_wakeconsumers(sqzthread);
    //        sleep(100);
    //    }
    //}
    exit:
        for (int i = 0; i < nthread; i++)
            if (pthread_join(consumer_pool[i], NULL))
                fprintf(stderr, "Thread error join\n");
        fprintf(stderr, "Cleaning up\n");
        sleep(10);
        return NULL;
}


char sqz_threadlauncher(FILE *ofp,
                        const char *filename,
                        char fqflag,
                        int nthread,
                        unsigned char fmt)
{
    //TODO Move to initialization routine
    sqzthread_t sqzthread;
    sqzthread.nthread = nthread;
    sqzthread.goread = 0;
    sqzthread.gocons = 0;
    sqzthread.doneq  = 0;
    sqzthread.wakethreadn = 0;
    sqzthread.threadid    = 0;
    sqzthread.fmt         = fmt;
    sqzthread.filename    = filename;
    sqzthread.sqzqueue    = calloc(nthread, sizeof(sqzfastx_t));
    for (int i = 0; i < nthread; i++)
        if (!sqz_fastxinit(sqzthread.sqzqueue + i, fmt, LOAD_SIZE))
            return 1;
    if (!sqzthread.sqzqueue) return 1;

    pthread_mutex_init(&(sqzthread.mtx), NULL);
    pthread_cond_init(&(sqzthread.conscond), NULL);
    pthread_cond_init(&(sqzthread.readcond), NULL);
    pthread_cond_init(&(sqzthread.intraconscond), NULL);

    pthread_attr_init(&(sqzthread.thatt));
    pthread_attr_setdetachstate(&(sqzthread.thatt), PTHREAD_CREATE_JOINABLE);

    fprintf(stderr, "This is where threads are launched\n");
    //Launch launch reader thread
    pthread_t rthread;
    if (pthread_create(&rthread,
                       &sqzthread.thatt,
                       sqz_readerthread,
                       (void *)&sqzthread))
        fprintf(stderr, "[sqz]: Thread error\n");
    //Wait for reader thred to finish
    if (pthread_join(rthread, NULL))
        fprintf(stderr, "\t[ERROR]: Thread error join.\n");


    //Clean up
    pthread_attr_destroy(&(sqzthread.thatt));
    pthread_mutex_destroy(&(sqzthread.mtx));
    pthread_cond_destroy(&(sqzthread.conscond));
    pthread_cond_destroy(&(sqzthread.readcond));
    pthread_cond_destroy(&(sqzthread.intraconscond));
    return 0;
    /*
    //Loop will be moved to multithread
    while ( (lbytes += sqz_loadfastX(sqz, fqflag)) > 0 ) {
        if (!sqz_fastXencode(sqz, blk, fqflag)) {
            fprintf(stderr, "[sqz ERROR]: Encoding error.\n");
            goto exit;
        }
        if (sqz->cmpflag) {
            numseqs += sqz->n;
            cbytes = sqz_deflate(blk, 9);
            if ( !sqz_zlibcmpdump(blk, cbytes, ofp) ) {
                fprintf(stderr, "[sqz ERROR]: Failed to write to output.\n");
                goto exit;
            }
            blk->blkpos  = 0;
            sqz->bases   = 0;
            sqz->namepos = 0;
            sqz->endflag = 0;
            sqz->cmpflag = 0;
            blk->newblk  = 1;
            sqz->n = 0;
            lbytes = 0;
        }
    }
    */
}
