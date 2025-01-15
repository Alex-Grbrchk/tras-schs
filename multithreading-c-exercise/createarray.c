#include "createarray.h"

void *CreateArray(void *v)
{
    int read_fd;
    int i, j;
    float _signal[COL];
    float p[LINE][COL];
    int _shmid;
    char *_addr;
    mqd_t mq;
    mq = mq_open(MQ_NAME, O_RDWR);
    if (mq == (mqd_t)-1)
    {
        perror("mq_open producer");
    }
    if ((_shmid = shmget((key_t)SHM_KEY, sizeof(p), IPC_CREAT | IPC_EXCL | 0666)) < 0)
    {
        if (errno != EEXIST)
        {
            perror("shmget send");
            exit(-1);
        }
        else
        {
            if ((_shmid = shmget((key_t)SHM_KEY, sizeof(p), 0666)) < 0)
            {
                printf("unable to find shared memory\n");
                exit(-1);
            }
        }
    }
    if ((_addr = shmat(_shmid, 0, 0)) == (void *)-1)
    {
        perror("shmat send");
        exit(-1);
    }
    if ((read_fd = open(PIPE_NAME, O_RDONLY)) <= 0)
    {
        perror("open");
    }

    while (1)
    {
        for (i = 0; i < LINE; i++)
        {
            read(read_fd, _signal, sizeof(_signal));
            if (mq_send(mq, (char *)_signal, sizeof(_signal), 0) == 1)
            {
                perror("mq_send");
            }

            for (j = 0; j < COL; j++)
            {
                *(*(p + i) + j) = *(_signal + j);
            }
        }
        memccpy(_addr, p, sizeof(p));
        sem_post(&sem);
    }
    close(read_fd);
    shmdt(_addr);
    shmctl(_shmid, IPC_RMID, 0);
    mq_close(mq);
    return NULL;
}