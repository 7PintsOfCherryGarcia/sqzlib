#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "../sqz_data.h"
#include "klib/kseq.h"
KSEQ_INIT(gzFile, gzread)

sqzblock_t *sqz_sqzblkinit(uint64_t size);
char sqz_fastxinit(sqzfastx_t *sqz, unsigned char fmt, uint64_t bsize);
uint64_t sqz_loadfastX(sqzfastx_t *sqz, uint8_t fqflag, kseq_t *seq);
char sqz_fastXencode(sqzfastx_t *sqz, sqzblock_t *blk, uint8_t fqflag);
size_t sqz_deflate(sqzblock_t *blk, int level);
char sqz_zlibcmpdump(sqzblock_t *blk, uint64_t size, FILE *ofp);


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
    char fqflag;
    FILE *ofp;
    //unsigned char endflag;
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
    //Thread wakes up, now it has to wait for all other threads to wake up
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


static void sqz_wakereader(sqzthread_t *sqzthread, int id)
{
    pthread_mutex_lock(&(sqzthread->mtx));
    sqzthread->doneq++;
    if (sqzthread->nthread == sqzthread->doneq) {
        fprintf(stderr, "\tthread: %d waking reader\n", id);
        sqzthread->doneq = 0;
        sqzthread->goread = 1;
        pthread_cond_signal(&(sqzthread->readcond));
    }
    sqzthread->gocons = 0;
    fprintf(stderr, "\tthread %d done\n", id);
    sqz_threadwait(&(sqzthread->gocons),
                   &(sqzthread->conscond),
                   &(sqzthread->mtx));
    //Thread wakes up and increases flag
    sqzthread->wakethreadn++;
    fprintf(stderr, "%d threads awake %d\n", sqzthread->wakethreadn, id);
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
    sqzfastx_t *sqz = sqzthread->sqzqueue + (id - 1);
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
        fprintf(stderr, "Compressed %d to %lu bytes\n", id, cbytes);
        //Write compressed block to output file
        pthread_mutex_lock(&(sqzthread->mtx));
        sqz_zlibcmpdump(blk, cbytes, sqzthread->ofp);
        fflush(sqzthread->ofp);
        fprintf(stderr, "Thread %d done writing\n", id);
        blk->blkpos  = 0;
        sqz->namepos = 0;
        sqz->endflag = 0;
        blk->newblk  = 1;
        pthread_mutex_unlock(&(sqzthread->mtx));
        if (sqz->endthread & 1) break; //bit 1 is on if thread had data, but reader has finished
        //We are done, we can signal reader and go to sleep while new data arrives
        sqz_wakereader(sqzthread, id);
        //fprintf(stderr, "Now what %d %u?!?\n", id, sqz->endthread);
    }
    fprintf(stderr, "Thread %d is out!!!\n", id);
    exit:
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
    sqzthread_t *sqzthread = thread_data;
    sqzfastx_t *sqzqueue   = sqzthread->sqzqueue;
    sqzthread->threadid++;
    int nthread = sqzthread->nthread;
    char fqflag = sqzthread->fqflag;
    //Initialize kseq object
    gzFile fp = gzopen(sqzthread->filename, "r");
    if (!fp) return NULL;
    kseq_t *seq = kseq_init(fp);
    if (!seq) return NULL;
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
    int thcounter = 0;
    //While there is data to read
    //Read data in blocks accorind to number of threads
    while (sqz_loadfastX(&(sqzqueue[thcounter]), fqflag, seq)) {
        fprintf(stderr, "Loaded in sqz #%d\n", thcounter);
        thcounter++;
        if (nthread == thcounter) {
            thcounter = 0;
            fprintf(stderr, "This is when threads will be woken up\n");
            sqz_wakeconsumers(sqzthread);
            fprintf(stderr, "The BEAST shall wakeup!!!\n");
        }
    }

    //TODO
    /*
      Now that reader is done, some signal must be sent to consumer threads
      indicating that they should exit once there is no more data. Problem is
      that there might be a mixture of consumers with data and others without.
      Those without should just terminate,while those with, should finish their
      consuming and then terminate.
    */
    fprintf(stderr, "Done with the reader loop\n");
    if (thcounter % nthread) {
        fprintf(stderr, "%d threads were not woken up\n", thcounter % nthread);
        //Set bit 1 of threads that got some data
        for (int i = 0; i < thcounter; i++)
            sqzqueue[i].endthread |= 1;
    }

    //Unset bit 7 from rest of threads
    for (int i = thcounter; i < nthread; i++) {
        //fprintf(stderr, "1 %d - %u\n", i, (unsigned char)sqzqueue[i].endthread);
        sqzqueue[i].endthread &= 127;
        //fprintf(stderr, "2 %d - %u\n\n", i, (unsigned char)sqzqueue[i].endthread);
    }
    fprintf(stderr, "Final waking should happen here\n");
    sqzthread->goread = 0;
    sqzthread->gocons = 1;
    sqzthread->wakethreadn = 0;
    pthread_cond_broadcast(&(sqzthread->conscond));

    exit:
        for (int i = 0; i < nthread; i++)
            if (pthread_join(consumer_pool[i], NULL))
                fprintf(stderr, "Thread error join\n");
        kseq_destroy(seq);
        gzclose(fp);
        fprintf(stderr, "Cleaning up, GTFO!!!\n");
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
    sqzthread.fqflag      = fqflag;
    sqzthread.filename    = filename;
    sqzthread.ofp         = ofp;
    sqzthread.sqzqueue    = calloc(nthread, sizeof(sqzfastx_t));
    if (!sqzthread.sqzqueue) return 1;
    for (int i = 0; i < nthread; i++)
        if (!sqz_fastxinit(sqzthread.sqzqueue + i, fmt, LOAD_SIZE))
            return 1;

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


static void sqz_singleread(kseq_t *seq)
{
    kseq_read(seq);
}
