#ifndef STRUCT_H
#define STRUCT_H

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

#endif