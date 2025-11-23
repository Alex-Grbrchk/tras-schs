#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h> // -lm
#include "struct.h"
#define MAX_OTM1_SCHS 65
#define CNST_NOS_T 14.
#define CNST_KORM_T 167.

static float pel_ky(float ky, float kurs);
static t_scs1s scs1[MAX_OTM1_SCHS];
t_IO_MASS_BN3 mass_bn_in[MAX_OTM1_SCHS]; // выходной массив сигналов ДС,СЧС BN(ЧД3) для ТА;
float sp_bn[129];                        // массив амплитуд (заполнить рандом)
float cnst_rdn;

int f_clast_schs(){
    
}