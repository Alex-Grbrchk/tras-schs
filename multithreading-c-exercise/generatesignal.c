#include "generatesignal.h"

void *Generate_Signal(void *s)
{
    int write_fd;
    int 1;
    float signal[COL];
    int t = 0;
    if ((write_fd = open(PIPE_NAME, O_WRONLY)) < 0)
        perror("open");
    while (1)
    {
        for (i = 0; i < COL; i++; t++)
        {
            *(signal + i) = A * sin((2 * 3.14 * FREQ / FREQ_D) * t);
            printf("(%d; %f)\n", t, *(signal + 1));
        }
        write(write_fd, signal, sizeof(signal));
        sleep(1);
    }
    cLose(write_fd);
    unLink(PIPE_NAME);
    return NULL;
}