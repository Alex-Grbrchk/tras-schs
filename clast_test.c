#include <stdio.h>
#include <string.h>
#include <time.h>

#define MAX_OTM1_SCHS 65

typedef struct
{
    int k; // реальный номер пространственного канала
    float sp;
    float s;
    float spmin; // минимальная амплитуда
    int kmin;
    int k1; // левый край непр. области
    int k2; // первый край непр. области
    int kt;
    int kes; // кол-во ЭС
    int chd; // номер частотного диапазона
    int nom; // номер структуры в массиве
} t_scs1s;

t_scs1s scs1[MAX_OTM1_SCHS];
float sp_bn[129]; // массив амплитуд

int main()
{
    int i, ii, jj;
    int j;         // для принтов
    int ng = 0;
    float fz = 50; // порого обнаружения
    int sp_index_max = 0;
    int otm;
    int param = 1;
    t_scs1s *yscs1;
    yscs1 = &scs1[0];
    srand(time(NULL));

    for (j = 0; j < 129; j++)
    {
        sp_bn[j] = (float)(rand() % 500);
        printf("%.f\t", sp_bn[j]);
    }
    printf("\n");

    while (param)
    {
        // Сравнение с порогом и поиск максимальной амплитуды перед формированием отметки
        for (i = 0; i < 129; i++)
        {
            if (sp_bn[i] > fz)
            {
                sp_index_max = i;
            }
            else
            {
                sp_bn[i] = 0;
            }
        }
        printf("fx:\n");
        for (j = 0; j < 129; j++)
        {
            printf("%.f\t", sp_bn[i]);
        }
        printf("\n");
        printf("Amax[%d] = %.f\n\n", sp_index_max, sp_bn[sp_index_max]);

        if (sp_bn[sp_index_max] == 0) // если Amax = 0 --> массив нулевой, все амплитуды рассмотрены
        {
            param = 0;
            break;
        }

        otm = 1;
        ii = 1;
        jj = 1;
        if (otm == 1)
        {
            ng++;
            yscs1->kes = 1;
            yscs1->k1 = sp_index_max;
            yscs1->k2 = sp_index_max;
            yscs1->k = sp_index_max;
            yscs1->kmin = sp_index_max;
            yscs1->sp = sp_bn[sp_index_max];
            yscs1->spmin = sp_bn[sp_index_max];
            yscs1->chd = 3;
            sp_bn[sp_index_max] = 0;
            while (otm)
            {
                if (sp_index_max - ii >= 0 && sp_bn[sp_index_max - ii] != 0)
                {
                    yscs1->kes++;
                    yscs1->k1 = sp_index_max - ii;
                    if (sp_bn[sp_index_max - ii] < yscs1->spmin)
                    {
                        yscs1->spmin = sp_bn[sp_index_max - ii];
                        yscs1->kmin = sp_index_max - ii;
                    }
                    sp_bn[sp_index_max - ii] = 0;
                    ii++;
                }
                else
                {
                    if (sp_index_max + jj < 129 && sp_bn[sp_index_max + jj] != 0)
                    {
                        yscs1->kes++;
                        yscs1->k2 = sp_index_max + jj;
                        if (sp_bn[sp_index_max + jj] < yscs1->spmin)
                        {
                            yscs1->spmin = sp_bn[sp_index_max + jj];
                            yscs1->kmin = sp_index_max + jj;
                        }
                        sp_bn[sp_index_max + jj] = 0;
                        jj++;
                    }
                    else
                    {
                        yscs1++;
                        otm = 1;
                    }
                }
                for (j = 0; j < 129; j++)
                {
                    printf("%.f\t", sp_bn[j]);
                }
                printf("\n");
            }
        }
    }
    yscs1 = &scs1[0];
    for (j = 0; j < ng; j++)
    {
        printf("%d\t", yscs1->kes);
        yscs1++;
    }
}
