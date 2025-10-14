#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h> // -lm
#define MAX_OTM1_SCHS 65
#define CNST_NOS_T 14.
#define CNST_KORM_T 167.

typedef struct
{
    int k;          // реальный номер пространственного канала;
    float kci;      // номер канала для вычисления tr_ct
    float ky_tr;    // курсовой угол
    float tr_ct;    // траверзный угол по центру тяжести
    float tr_amax;  // траверзный угол - по Amax;
    int pr_ct_amax; // признак тр угла по максимуму
    int kmin;       // минимальное значение сигнал/помеха
    float sp;       // амплитуда сигнала нормированная и центрированная;
    float s;        // амплитуда сигнала центрированная;
    float spmin;    // минимальная А (для поиска акустической протяж)
    float prg_prvl; // порог провала
    int k1;         // левый край непр. обл
    int k2;         // правый край непр. обл.
    int kp;         // карйний правый канал, при котором А >= 0.5 Amax
    int kl;         // крайний левый канал, при котором А >= 0.5 Amax
    int kp5;        // kp (для дс)
    int kl5;        // kl (для дс)
    int kt;         //
    int kes;        // кол-во эс
    int chd;        // номер частотного диапазона
    int nom;        // номер структуры в массиве
} t_scs1s;

typedef struct
{

    int chd;      // номер ЧД (1/2/3 )
    float sp;     // максимальная нормированная величина(ОСП)
    float s;      // максимальная ненормированная(амплитуда)
    int k;        // номер пространственного канала(максимум sp)
    int prot;     // протяженность
    float a_prot; // акустическая протяженность
    float p_prot;
    int q;         // номер частотного отсчета(максимум sp для отметок ПЧС)
    float sh_ds5;  // ширина ДС на уровне 0.5
    float ky;      // курсовой угол
    float cko_pel; // СКО пеленга
    float f;       // частота
    float kci;
    float ky_tr; // траверзный угол(рад) Qq=arcsin(2q/(Q-1));
                 //   q=-64,...,0,...,64;
                 //   Q=129 - количество ХН;

    int pr_pchs_schs; // признак ДС/СЧС (0/1);
    int pr_brt;       // признак определенности борта (0/1/2-Н/ЛБ/ПБ);
    float cko_f;      // СКО частоты;
    float KA;         // коэффициент асимметрии;
    int KOL_PRVL;     // число провалов;
    int var_schs;     // вводится оператором (для СЧС по МН);
    float d_dfrnt;    // поправка на дифферент;
    int nomscol;
    int pr_korma;
    float pel;  // пеленг;
    float zpel; // зеркальный пеленг;
    int notm;   // номер отметки;

} t_str_pchs_schs;

typedef struct
{
    // t_IO_NAV_PAST mass_nav;
    int num_form_ds;
    t_str_pchs_schs formul_ds[100];
    int num_form_schs;
    t_str_pchs_schs formul_schs[30];
} t_IO_MASS_BN3;

static float pel_ky(float ky, float kurs);
static t_scs1s scs1[MAX_OTM1_SCHS];
t_IO_MASS_BN3 mass_bn_in[MAX_OTM1_SCHS]; // выходной массив сигналов ДС,СЧС BN(ЧД3) для ТА;
float sp_bn[129];                        // массив амплитуд (заполнить рандом)
float cnst_rdn;

int main()
{
    int i, ii, jj;
    int j;      // для принтов
    int ng = 0; // счётчик отметок
    float SA;   // Сумма амплитуд элементарных сигналов одной отметки;
    float SAl;  // Сумма произведений амплитуд элементарных сигналов на номер пространственного канала одной отметки;
    float kci;  // Координата центра тяжести(пространственный канал) отметки СЧС;
    float fz = 50;
    float arg = 0.;
    float kfct_prvl = .1; // коэффициент для вычисления числа провалов отметки СЧС;
    int sp_index_max = 0;
    int otm;
    int param = 1;
    // float d_dfrnt[129];     // поправки на дифферент (по 129 направлениям), поступающие из системы позиционирования;
    t_scs1s *yscs1;         // Указатель на массив первичных параметров отметок СЧС;
    t_str_pchs_schs *ytscs; // указатель на массив параметров отметки СЧС;
    yscs1 = &scs1[0];
    ytscs = &mass_bn_in->formul_schs[0];
    cnst_rdn = M_PI / 180.;

    srand(time(NULL));
    for (j = 0; j < 129; j++)
    {
        sp_bn[j] = (float)(rand() % 500);
        printf("%.0f\t", sp_bn[j]);
    }
    printf("\n");

    clock_t t = clock();
    while (param)
    {
        // сравнение с порогом и поиск максимальной амплитуды перед формированием отметки
        for (i = 0; i < 129; ++i)
        {
            if (sp_bn[i] > fz)
            {
                if (sp_bn[i] > sp_bn[sp_index_max])
                {
                    sp_index_max = i;
                }
            }
            else
            {
                sp_bn[i] = 0;
            }
        }
        /*printf("fz:\n");
        for (j = 0; j < 129; j++)
        {
            printf("%.0f\t", sp_bn[j]);
        }
        printf("\n");
        printf("Amax[%d] = %.f\n\n", sp_index_max, sp_bn[sp_index_max]);*/

        if (sp_bn[sp_index_max] == 0) // если Аmax = 0 --> массив нулевой, все амплитуды рассмотрены
        {
            param = 0;
            break;
        }

        otm = 1;
        ii = 1;
        jj = 1;
        if (otm == 1) // если в отметке
        {
            // первая амплитуда в отметке
            SA = 0.;
            SAl = 0.;
            kci = 0.;
            yscs1->nom = ng + 1;
            yscs1->kes = 1;
            yscs1->k1 = sp_index_max;
            yscs1->k2 = sp_index_max;
            yscs1->k = sp_index_max;
            yscs1->kmin = sp_index_max;
            yscs1->sp = sp_bn[sp_index_max];
            // yscs1->s = s1_bn[sp_index_max - ii];
            yscs1->spmin = sp_bn[sp_index_max];
            yscs1->tr_amax = asin(2. * (yscs1->k - 65.) / 128.); // траверзный угол - по Amax;
            yscs1->chd = 3;

            SA = SA + sp_bn[sp_index_max];
            SAl = SAl + sp_bn[sp_index_max] * sp_index_max;
            sp_bn[sp_index_max] = 0; // обнуляем максимальную амплитуду, вокруг которой строится отметка
            ng++;

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
                    SA = SA + sp_bn[sp_index_max - ii];
                    SAl = SAl + sp_bn[sp_index_max - ii] * (sp_index_max - ii);

                    sp_bn[sp_index_max - ii] = 0; // обнуляем рассмотренную амплитуду
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
                        SA = SA + sp_bn[sp_index_max + jj];
                        SAl = SAl + sp_bn[sp_index_max + jj] * (sp_index_max + jj);

                        sp_bn[sp_index_max + jj] = 0; // обнуляем рассмотренную амплитуду
                        jj++;
                    }
                    else // построение отметки закончилось
                    {
                        /*kci = SAl / SA;
                        yscs1->kci = kci;
                        arg = 2. * (kci - 65.) / 128.;
                        if (fabs(arg) > 1.)
                        {
                            if (arg < 0.)
                                arg = -1.;
                            else
                                arg = 1.;
                        }
                        yscs1->tr_ct = asin(arg); // траверзный угол по центру тяжести

                        // если траверзный угол вошел в зону продольных углов
                        
                        if (yscs1->tr_amax <= cnst_rdn * (CNST_NOS_T - 90.) || yscs1->tr_amax >= cnst_rdn * (CNST_KORM_T - 90.))
                        {
                            yscs1->ky_tr = yscs1->tr_amax;
                            yscs1->pr_ct_amax = 1;
                        }
                        else
                        {
                            yscs1->ky_tr = yscs1->tr_ct;
                            yscs1->pr_ct_amax = 0;
                        }

                        if (yscs1->ky_tr <= cnst_rdn * (CNST_NOS_T - 90.) || yscs1->ky_tr >= cnst_rdn * (CNST_KORM_T - 90.))
                        {
                            yscs1->nom = 0;
                        }
                        yscs1->prg_prvl = kfct_prvl * (yscs1->sp - yscs1->spmin);

                        if (yscs1->nom != 0)
                        {
                            ytscs->chd = yscs1->chd;
                            ytscs->ky_tr = yscs1->ky_tr;
                            ytscs->k = yscs1->k;
                            ytscs->prot = yscs1->k2 - yscs1->k1 + 1;
                            ytscs->sp = yscs1->sp;
                            // ytscs->s = yscs1->s;
                            ytscs->pr_pchs_schs = 1;
                            ytscs->ky = ytscs->ky_tr + M_PI / 2.;
                        }
                        */

                        yscs1++; // следующая отметка
                        otm = 0;
                    }
                }
                for (j = 0; j < 129; j++)
                {
                    printf("%.0f\t", sp_bn[j]);
                }
                printf("\n");
            }
        }
    }
    t = clock() - t;
    printf("%lf seconds\n", (double)(t) / CLOCKS_PER_SEC);

    yscs1 = &scs1[0];
    for (j = 0; j < ng; j++)
    {
        printf("%d\t", yscs1->kes);
        yscs1++;
    }
}

float pel_ky(float ky, float kurs)
{

    float pel;

    pel = kurs + ky;
    if (pel < 0.)
        pel = pel + 2. * M_PI;
    else if (pel > 2. * M_PI)
        pel = pel - 2. * M_PI;
    return (pel);
}