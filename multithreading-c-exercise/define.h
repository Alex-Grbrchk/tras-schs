#ifndef _DEFINE_VALUES 
#define _DEFINE_VALUES 1
#define FREQ 2500
#define FREQ_D 10000
#define A 1
#define COL 16
#define LINE 4
#define MAX_BUF 64
#define PIPE_NAME "/tmp/named_pipe"
#define MQ_NAME "/arr_queue"
#define SHM_KEY 2081
extern sem_t sem;
extern pthread_mutex_t m;
#endif