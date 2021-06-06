#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>

typedef struct {
    pthread_mutex_t mtx;
    pthread_cond_t conscond;
    pthread_cond_t conscond2;
    pthread_cond_t readcond;
    int goread;
    int gocons;
    int doneq;
    int wthread;
    int threadid;
} threads_t;

void wait(int *flag, pthread_cond_t *cond, pthread_mutex_t *mtx)
{
    while (!(*flag)) {
        pthread_cond_wait(cond, mtx);
    }
}


/*
Same as wait, but indicates how many threads have been woken up
*/
void wait2(int *flag, pthread_cond_t *cond, pthread_mutex_t *mtx)
{
    if (10 == *flag) {
        fprintf(stderr, "broadcasting!!!\n");
        pthread_cond_broadcast(cond);
        return;
    }
    while (*flag != 10)
        pthread_cond_wait(cond, mtx);

}


int getqueue(threads_t *thread_data)
{
    pthread_mutex_lock(&(thread_data->mtx));
    int id = thread_data->threadid;
    thread_data->threadid += 1;
    wait(&(thread_data->gocons),
         &(thread_data->conscond),
         &(thread_data->mtx));
    //Thread wakes up and increases flag
    thread_data->wthread++;
    fprintf(stderr, "%d threads awake %d\n", thread_data->wthread, id);
    //Thread goes to sleep again until all threads have woken up
    wait2(&(thread_data->wthread),
          &(thread_data->conscond2),
          &(thread_data->mtx));

    pthread_mutex_unlock(&(thread_data->mtx));
    return id;
}


void coordinate_consumers(threads_t *thread_data)
{
    pthread_mutex_lock(&(thread_data->mtx));
    thread_data->goread = 0;
    thread_data->gocons = 1;
    thread_data->wthread = 0;
    pthread_cond_broadcast(&(thread_data->conscond));
    wait(&(thread_data->goread),
         &(thread_data->readcond),
         &(thread_data->mtx));
    pthread_mutex_unlock(&(thread_data->mtx));
}


void coordinate_reader(threads_t *thread_data, int id)
{
    pthread_mutex_lock(&(thread_data->mtx));
    thread_data->doneq++;
    if (10 == thread_data->doneq) {
        fprintf(stderr, "\tthread: %d waking reader\n", id);
        thread_data->doneq = 0;
        thread_data->goread = 1;
        pthread_cond_signal(&(thread_data->readcond));
    }
    thread_data->gocons = 0;
    fprintf(stderr, "\tthread %d done\n", id);
    wait(&(thread_data->gocons),
         &(thread_data->conscond),
         &(thread_data->mtx));
    //Thread wakes up and increases flag
    thread_data->wthread++;
    fprintf(stderr, "%d threads awake %d\n", thread_data->wthread, id);
    //Thread goes to sleep again until all threads have woken up
    wait2(&(thread_data->wthread),
          &(thread_data->conscond2),
          &(thread_data->mtx));
    pthread_mutex_unlock(&(thread_data->mtx));
}


void *consumer(void *thread_data)
{
    threads_t *data = thread_data;
    int id = getqueue(data);
    fprintf(stderr, "Got id: %d\n", id);
    while(1) {
        //Does some work
        coordinate_reader(data, id);
        //Wait for all other threads to wake up
    }
}


void *reader(void *thread_data)
{
    threads_t *data = thread_data;
    fprintf(stderr, "goread: %d\ngocond: %d\n", data->goread, data->gocons);
    int id = data->threadid;
    data->threadid += 1;
    fprintf(stderr, "threadid: %d\n", id);

    pthread_attr_t thattr;
    pthread_attr_init(&thattr);
    pthread_attr_setdetachstate(&thattr, PTHREAD_CREATE_JOINABLE);

    pthread_t *consumer_pool = NULL;
    consumer_pool = malloc(10 * sizeof(pthread_t));
    for (int i = 0; i < 10; i++)
        pthread_create(consumer_pool+i, &thattr, consumer, thread_data);
    fprintf(stderr, "Whileing\n");
    int counter = 0;
    while (1) {
        counter++;
        if (counter == 10) {
            counter = 0;
            fprintf(stderr, "Coordinating consumers\n");
            coordinate_consumers(data);
            fprintf(stderr, "Reader woke up\n");
        }
    }
}


int main()
{
    srand(time(NULL));
    threads_t tdata;
    tdata.gocons = 0;
    tdata.goread = 0;
    tdata.doneq  = 0;
    tdata.wthread  = 0;
    tdata.threadid = 0;
    pthread_mutex_init(&(tdata.mtx), NULL);
    pthread_cond_init(&(tdata.conscond), NULL);
    pthread_cond_init(&(tdata.conscond2), NULL);
    pthread_cond_init(&(tdata.readcond), NULL);

    pthread_attr_t thattr;
    pthread_attr_init(&thattr);
    pthread_attr_setdetachstate(&thattr, PTHREAD_CREATE_JOINABLE);

    pthread_t read_thread;
    if (pthread_create(&read_thread, &thattr, reader, (void *)&tdata))
        fprintf(stderr, "[ERROR]: Thread creation.\n");
    // Wait for threads to finish
    if (pthread_join(read_thread, NULL))
        fprintf(stderr, "\t[ERROR]: Thread error join.\n");

    pthread_mutex_destroy(&(tdata.mtx));
}



