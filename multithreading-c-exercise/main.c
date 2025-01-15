#include "includes.h"

sem_t sem;
pthread_mutex_t m;
int main
{
    int i;
    pthread_t thr[4];
    sem_init(&sem, 0, 0);
    pthread_mutex_init(&m, NULL);
    remove(PIPE_NAME);
    if (mkfifo(PIPE_NAME, 0777) == -1)
    {
        perror("mkfifo");
        return 1;
    }
    pthread_create(&thr[0], NULL, Generate_Signal, 0);
    pthread_create(&thr[1], NULL, CreateArray, 0);
    pthread_create(&thr[2], NULL, Transpose, 0);
    pthread_create(&thr[3], NULL, Search, 0);
    for (i = 0; i < 4; i++)
    {
        pthread_join(thr[i], NULL);
    }

    pthread_mutex_destroy(&m);
    sem_destroy(&sem) :
}
// gcc main. generatesignal. createarray. transpose. search. -o main -m -pthread -Irt