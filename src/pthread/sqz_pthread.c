#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <pthread.h>

#include "../sqzlib/sqzlib.h"

#include "klib/kseq.h"
int64_t sqzread(sqzFile file, void *buff, uint64_t len);
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
        if ( !(sqzqueue[i] = sqz_fastxinit(fmt, LOAD_SIZE)) ) {
            fprintf(stderr, "ERROR\n");
            for(int j = 0; j < i; j++)
                sqz_fastxkill(sqzqueue[j]);
            return NULL;
        }
    fprintf(stderr, "INIT endflag: %u\n", sqzqueue[0]->endflag);
    sleep(10);
    return sqzqueue;
}


static void *sqz_consumerthread(void *thread_data)
{
    sqzthread_t *sqzthread = thread_data;
    uint64_t cbytes = 0;
    sqzblock_t *blk = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;
    int id = sqz_getthreadid(sqzthread);
    sqzfastx_t *sqz = sqzthread->sqzqueue[id - 1];
    uint8_t fqflag  = sqzthread->fqflag;
    int libfmt      = sqzthread->libfmt;
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
        cbytes = sqz_blkcompress(blk, 9, libfmt);
        //TODO Deal with error
        //Write compressed block to output file
        pthread_mutex_lock(&(sqzthread->mtx));
        sqz_blkdump(blk->cmpbuff, &(blk->blkpos), cbytes, sqzthread->ofp);
        fflush(sqzthread->ofp);
        blk->blkpos  = 0;
        blk->newblk  = 1;
        sqz->namepos = 0;
        sqz->endflag = 0;
        sqz->blks++;
        pthread_mutex_unlock(&(sqzthread->mtx));
        //bit 1 is on if thread had data, but reader has finished
        if (sqz->endthread & 1)
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


static void *sqz_readerthread(void *thread_data)
{
    sqzthread_t *sqzthread = (sqzthread_t *)thread_data;
    sqzthread->threadid++;

    int tmp = 5;
    sqzfastx_t **sqzqueue = NULL;
    kseq_t *seq = NULL;
    pthread_t *consumer_pool;

    sqzFile sqzfp = sqzopen(sqzthread->ifile, "r");
    if (!sqzfp) goto exit;

    uint8_t nthread = sqzthread->nthread;
    fprintf(stderr, "$$$%u\n", sqzfp->fmt);
    sqzqueue = sqz_sqzqueueinit(nthread, sqzfp->fmt & 3);
    if (!sqzqueue)goto exit;
    fprintf(stderr, "endflag: %u\n", sqzqueue[0]->endflag);
    sleep(10);

    sqzthread->sqzqueue = sqzqueue;
    uint8_t fqflag  = (sqzfp->fmt == 2) ? 1 : 0;;

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
                           (void *)thread_data))
            goto exit;
    int thcounter = 0;
    //TODO: Error checking on fastX parsing
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
            fprintf(stderr, "[sqz ERROR]: Thread error join\n");
    //Log number of sequences and blocks
    uint64_t n = 0;
    uint64_t b = 0;
    for (int i = 0; i < nthread; i++) {
        b += sqzqueue[i]->blks;
        n += sqzqueue[i]->n;
    }
    sqz_filetail(n, b, sqzthread->ofp);
    tmp = 0;
    exit:
    fprintf(stderr, "tmp: %d\n", tmp);
        if (consumer_pool) free(consumer_pool);
        if (seq) kseq_destroy(seq);
        if (sqzfp) sqzclose(sqzfp);
        return NULL;
}


static sqzthread_t *sqz_threadinit(const char *ifile,
                                   const char *ofile,
                                   uint8_t libfmt,
                                   uint8_t nthread)
{
    sqzthread_t *sqzthread = calloc(1, sizeof(sqzthread_t));
    if (!sqzthread) return NULL;
    sqzthread->nthread = nthread;
    sqzthread->ifile   = ifile;
    sqzthread->ofile   = ofile;
    sqzthread->libfmt  = libfmt;
    pthread_mutex_init(&(sqzthread->mtx), NULL);
    pthread_cond_init(&(sqzthread->conscond), NULL);
    pthread_cond_init(&(sqzthread->readcond), NULL);
    pthread_cond_init(&(sqzthread->intraconscond), NULL);
    pthread_attr_init(&(sqzthread->thatt));
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


uint8_t sqz_threadcompress(const char *ifile,
                           const char *ofile,
                           uint8_t libfmt,
                           uint8_t nthread)
{
    uint8_t ret = 1;
    //uint8_t fmt = sqz_getformat(ifile);
    //if (!fmt) {
    //    fprintf(stderr, "[sqz ERROR]: Not a fastX file.\n");
    //    return ret;
    //}
    //Start thread object
    sqzthread_t *sqzthread = sqz_threadinit(ifile,
                                            ofile,
                                            libfmt,
                                            nthread);
    if (!sqzthread) goto exit;
    //Launch reader thread
    pthread_t rthread;
    if (pthread_create(&rthread,
                       &(sqzthread->thatt),
                       sqz_readerthread,
                       (void *)sqzthread)) {
        fprintf(stderr, "[sqz]: Error Thread launch error.\n");
        goto exit;
    }
    //Wait for reader thred to finish
    if (pthread_join(rthread, NULL)) {
        fprintf(stderr, "\t[sqz]: Error, thread join failed.\n");
        goto exit;
    }
    ret = 0;
    exit:
        //Clean up
        sqz_threadkill(sqzthread);
        return ret;
}


uint8_t sqz_inflatefastX(FILE *ifp, FILE *ofp, char fqflag, uint8_t libfmt)
{
    uint8_t ret      = 1;
    uint8_t *outbuff = NULL;
    uint64_t dsize   = 0;
    int64_t size    = 0;
    sqzblock_t *blk  = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;
    size = sqz_filesize(ifp);
    outbuff = malloc(LOAD_SIZE);
    if (!outbuff) goto exit;
    fseek(ifp, HEADLEN, SEEK_SET);
    while ( ftell(ifp) < size )
        {
            if (!sqz_readblksize(blk, ifp, libfmt)) goto exit;
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


static uint8_t sqz_gzread(const char *ifile, const char *ofile)
{
    uint8_t ret = 1;
    uint8_t *buff = NULL;
    FILE *ofp = NULL;
    gzFile ifp = gzopen(ifile, "r");
    if (!ifp) goto exit;
    ofp = fopen(ofile, "w");
    buff = malloc(LOAD_SIZE);
    if (!buff) goto exit;
    uint32_t read;
    while ( (read = gzread(ifp, buff, LOAD_SIZE)) ) {
        fwrite(buff, read, 1, ofp);
    }
    ret = 0;
    exit:
         if (ifp) gzclose(ifp);
         if (ofp) fclose(ofp);
         if (buff) free(buff);
         return ret;
}


uint8_t sqz_threaddecompress(const char *ifile,
                             const char *ofile,
                             uint8_t nthread)
{
    uint8_t ret = 1;
    sqzFile sqzfp = sqzopen(ifile, "rb");

    if ( !(sqzfp->fmt & 4) ) {
        fprintf(stderr, "[sqz WARNING]: Not an sqz file\n");
        return sqz_gzread(ifile, ofile);
    }
    uint64_t nblk;
    if (sqz_getblocks(sqzfp, &nblk)) return ret;
    fprintf(stderr, "%lu blocks in file\n", nblk);
    sleep(10);
    FILE *ifp = fopen(ifile, "rb");
    FILE *ofp = fopen(ofile, "wb");
    if (sqz_inflatefastX(ifp, ofp, 1, 2)) fprintf(stderr, "^ERROR\n");
    fclose(ifp);
    fclose(ofp);
    ret = 0;
    return ret;
}
/*
char ret = 0;
switch (fmt & 7) {
case 1:
if (!sqz_deflatefastX(opts.ifile,
opts.ofile,
0,
opts.libfmt,
opts.nthread)) goto exit;
break;
case 2:
if (!sqz_deflatefastX(opts.ifile,
opts.ofile,
1,
opts.libfmt,
opts.nthread)) goto exit;
break;
case 5:
fprintf(stderr, "File %s alredy sqz encoded.\n", opts.ifile);
goto exit;
case 6:
fprintf(stderr, "File %s alredy sqz encoded.\n", opts.ifile);
goto exit;
}
ret = 1;
exit:
return ret;
*/


/*
char sqz_decompress(sqzopts_t opts)
{
char ret = 0;
FILE *ifp = NULL;
FILE *ofp = NULL;
ofp = fopen(opts.ofile, "wb");
if (!ofp)
goto exit;
ifp = fopen(opts.ifile, "rb");
if (!ifp)
goto exit;
uint8_t fmt = sqz_getformat(opts.ifile);
uint8_t libfmt = fmt >> 3;
//Check for format
switch (fmt & 7) {
case 1:
fprintf(stderr, "File %s already decoded.\n", opts.ifile);
goto exit;
case 2:
fprintf(stderr, "File %s already decoded.\n", opts.ifile);
goto exit;
case 5:
if (!sqz_inflatefastX(ifp, ofp, 0, libfmt)) {
fprintf(stderr, "[sqz ERROR]: Failed to decode data.\n");
goto exit;
}
break;
case 6:
if (!sqz_inflatefastX(ifp, ofp, 1, libfmt)) {
fprintf(stderr, "[sqz ERROR]: Failed to decode data.\n");
goto exit;
}
break;
}
ret = 1;
exit:
if(ifp) fclose(ifp);
if(ofp) fclose(ofp);
return ret;
}
*/


/*

*/
