#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <pthread.h>

#include "../sqzlib/sqzlib.h"

#include "klib/kseq.h"
//int64_t sqzread(sqzFile file, void *buff, uint64_t len);
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


static void sqz_gzread(sqzFile sqzfp, const char *ofile)
{
    uint8_t *buff = NULL;
    FILE *ofp = NULL;
    ofp = fopen(ofile, "w");
    if (!ofp) goto exit;
    buff = malloc(LOAD_SIZE);
    if (!buff) goto exit;
    uint32_t read;
    while ( (read = gzread(sqzfp->gzfp, buff, LOAD_SIZE)) ) {
        fwrite(buff, read, 1, ofp);
    }
    exit:
        if (ofp) fclose(ofp);
        if (buff) free(buff);
}


uint8_t sqz_go2blockn(sqzFile sqzfp, uint64_t n)
{
    //TODO Error out on request of block number greater than blocks available
    //fprintf(stderr, "\tRequested blk: %lu\n", n);
    uint64_t blkn = 0;
    uint64_t blks = 0;
    uint64_t blkr = 0;
    if (fseek(sqzfp->fp, HEADLEN, SEEK_SET)) return 1;
    while (blkn < n) {
        if (fseek(sqzfp->fp, 8, SEEK_CUR)) return 1;
        blkr = fread(&blks, B64, 1, sqzfp->fp);
        if (blkr) {
            //fprintf(stderr, "\t\tblk: %lu blksize: %lu\n",blkn ,blks);
            if(fseek(sqzfp->fp, blks, SEEK_CUR)) return 1;
            sqzfp->filepos = ftell(sqzfp->fp);
            blkn++;
        }
        else {
            return 1;
        }
    }
    return 0;
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
    uint8_t fqflag = (sqzfp->fmt & 3) == 2 ? 1 : 0;
    uint8_t libfmt = sqzfp->fmt >> 3;
    sqzblock_t *blk =  sqzfp->blk;

    outbuff = malloc(LOAD_SIZE);
    if (!outbuff) goto exit;

    filebuff = malloc(LOAD_SIZE);
    uint64_t fbsize = LOAD_SIZE;
    uint64_t fbpos = 0;

    uint64_t nblk;
    if (sqz_getblocks(sqzfp, &nblk)) goto exit;;

    uint64_t currentblk = (uint64_t)id;
    //currentblk = 395;
    while (currentblk < nblk) {
        //fprintf(stderr, "currentblk: %lu\n", currentblk);
        if (sqz_go2blockn(sqzfp, currentblk)) {
            fprintf(stderr, "ERROR\n");
            sleep(100);
        }
        if (!sqz_readblksize(blk, sqzfp->fp, libfmt))
            goto exit;
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

    uint8_t fmt = sqzfp->fmt;
    if ( !(fmt & 4) ) {
        fprintf(stderr, "[sqz WARNING]: Not an sqz file\n");
        sqz_gzread(sqzfp, sqzthread->ofile);
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


static void *sqz_compressor(void *thread_data)
{
    sqzthread_t *sqzthread = (sqzthread_t *)thread_data;
    sqzthread->threadid++;

    sqzfastx_t **sqzqueue = NULL;
    kseq_t *seq = NULL;
    pthread_t *consumer_pool = NULL;

    sqzFile sqzfp = sqzopen(sqzthread->ifile, "r");
    if (!sqzfp) goto exit;

    uint8_t nthread = sqzthread->nthread;
    sqzqueue = sqz_sqzqueueinit(nthread, sqzfp->fmt & 3);
    if (!sqzqueue) goto exit;
    sqzthread->sqzqueue = sqzqueue;

    uint8_t fqflag  = (sqzfp->fmt == 2) ? 1 : 0;;
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
    exit:
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
    sqzthread->ifile   = ifile;
    sqzthread->ofile   = ofile;
    sqzthread->libfmt  = libfmt;
    sqzthread->nthread = nthread;
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


static uint8_t sqz_inflatefastX(sqzFile sqzfp, FILE *ofp, char fqflag, uint8_t libfmt)
{
    uint8_t ret      = 1;
    uint8_t *outbuff = NULL;
    uint64_t dsize   = 0;
    int64_t size    = 0;
    sqzblock_t *blk  = sqz_sqzblkinit(LOAD_SIZE);
    if (!blk) goto exit;
    size = sqz_filesize(sqzfp);
    outbuff = malloc(LOAD_SIZE);
    if (!outbuff) goto exit;
    fseek(sqzfp->fp, HEADLEN, SEEK_SET);
    while ( ftell(sqzfp->fp) < size )
        {
            if (!sqz_readblksize(blk, sqzfp->fp, libfmt)) goto exit;
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


uint8_t sqz_threaddecompress(const char *ifile,
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


uint8_t sqz_threadcompress(const char *ifile,
                           const char *ofile,
                           uint8_t libfmt,
                           uint8_t nthread)
{
    uint8_t ret = 1;
    sqzthread_t *sqzthread = sqz_threadinit(ifile,
                                            ofile,
                                            libfmt,
                                            nthread);
    if (!sqzthread) goto exit;
    pthread_t rthread;
    if (pthread_create(&rthread,
                       &(sqzthread->thatt),
                       sqz_compressor,
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


