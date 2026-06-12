#include "struct.h"  // тут arxs
#include "io_func.h" // тут t_str_pchs_schs
#include <math.h>

#define DLA 100
#define DLS 100
#define NM 100
#define CHSLDC 1
#define lz 3.5
#define PISR 3.14159265
#define PPI 6.28318530
#define STRB1 2.
#define MNJ PISR / 180.
static float cnst[3];     // частота дискретизации / 1024
static int iden[DLA];     // номер отметки для данного архивного сигнала
static int znomer[DLA];   // массив-счётчик кол-ва отметок для каждого сигнала (znomer[i]+1 - номер сигнала, i - кол-во отметок)
static int zmas[DLA][NM]; // массив отметок, проидентифицированных с сигналами (№сигнала = №строки+1, zmas[i][j] - №отметки)
static int mass_bn_a[DLS][NM];
static int progbrt[DLA][NM], prbrt[DLA][NM];
static int hhel[NM];      // номер отметки без трассы, если noms[i]==0, записывается как новая в банк трасс на свободное местось
static int svn[DLS];      // признак идентификации отметок (svn[i] = 0 или 1), i + 1 - номер отметки
static int prizn_kon_man; // из f_manev.c, переписывается по началу вызова функции
static str_svs *ykta;     // указатель на структуру курсовых данных

float ky_pel(float, float);

int idnt(arxs *parx, t_str_pchs_schs *asgn, int adl,
         int sdl, int kols, int klz)
{
    int i, j, k, ll;
    int min_srv; // минимальная разница пеленгов
    int frst;
    int min_z; // ближайшая отметка для трассы
    int pr_perex1[100][100], pr_perex2[100][100];
    int pr_idnt;
    int pr_tras1, pr_tras2; // признаки трасс
    int pr_strpel;          // признак вычисления размера строба по пеленгу
    float strob;
    float ky_tr1, ky_tr2;
    double srv, ssrv;     // разность прогнозируемого пеленга и измеренного
    double spsrv, sspsrv; // квадрат ошибки пеленга
    t_str_pchs_schs *sgn; // набор отметок целей
    arxs *arx;            // банк трасс
    int noms[DLS];        // кол-во архивных сигналов для данной нерпрерывной области
    int siden[DLA];       // признак обработки архивного сигнала (0 или 1), индекс соотв. индексу в массиве банка трасс
    int hnm[NM];
    int imas[DLA]; // номера трасс, проидентифицированных с отметками, znomer[i] != 0 => запись номера i-ой трассы (i+1), используется для завершения идентификации

    cnst[0] = 1250. / 1024.;
    cnst[1] = 2500. / 1024.;
    cnst[2] = 10000. / 1024.;
    // обнуление массивов
    for (i = 0; i < adl; i++)
    {
        znomer[i] = 0;
        iden[i] = 0;
        siden[i] = 0;
        imas[i] = 0;
    }
    for (i = 0; i < adl; i++)
    {
        for (j = 0; j < NM; j++)
        {
            prbrt[i][j] = 0;
            progbrt[i][j] = 0;
            zmas[i][j] = 0;
        }
    }
    for (i = 0; i < sdl; i++)
        for (j = 0; j < NM; j++)
            mass_bn_a[i][j] = 0;

    for (i = 0; i < NM; i++)
    {
        hnm[i] = 0;
        hhel[i] = 0;
    }

    for (i = 0; i < sdl; i++)
    {
        noms[i] = 0;
        svn[i] = 0;
    }

    sgn = asgn;
    if (prizn_kon_man == 0)
    {
        for (k = 0; k < kols; k++, sgn++) // кол-во отметок
        {
            arx = parx;
            for (i = 0; i < adl; i++, arx++) // кол-во трасс в архиве
            {
                pr_strpel = 0;
                if (arx->anom != 0) // если трасса не пустая
                {
                    if (sgn->pr_pchs_schs != 1) // если ДС
                    {
                        if (fabs(arx->f - sgn->f) < CHSLDC * cnst[arx->chd - 1]) // вычисление строба по частоте
                        {
                            pr_strpel = 1;
                        }
                    }
                    else
                    {
                        pr_strpel = 1;
                    }

                    if (pr_strpel == 1)
                    {
                        if (arx->prbort == 3) // неоднозначность не устранена, рассматриваем обе трассы
                        {
                            pr_idnt = 0;
                            pr_perex1[i][k] = 0;
                            pr_perex2[i][k] = 0;
                            ky_tr1 = ky_pel(arx->PelProgn * 180. / M_PI, ykta->kyrs * 180. / M_PI);
                            ky_tr2 = ky_pel(arx->PelPrognZ * 180. / M_PI, ykta->kyrs * 180. / M_PI);

                            if (ky_tr1 >= 0.)
                                pr_tras1 = 1;
                            else
                                pr_tras1 = 2;
                            if (ky_tr2 >= 0.)
                                pr_tras2 = 1;
                            else
                                pr_tras2 = 2;

                            if (pr_tras1 == 1)
                            {
                                srv = arx->PelProgn - sgn->pel;
                                if (fabs(srv) > PISR)
                                {
                                    pr_perex1[i][k] = 1;
                                    srv = PPI - fabs(srv);
                                }
                                else
                                    srv = fabs(srv);
                                spsrv = arx->sigmaPP4;
                                strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel); // размер строба по пеленгу
                                // проверка строба if (strob < min) strob = min; if (strob > max) strob = max;
                                if (strob < STRB1 * MNJ)
                                    strob = STRB1 * MNJ;
                                if (strob > STRB1 * MNJ)
                                    strob = STRB1 * MNJ;

                                if (srv <= strob) // проверка попадания в строб отметки
                                {
                                    // добавить отметку в массив отметок, увеличить счётчик отметок сигнала
                                    zmas[i][znomer[i]] = sgn->notm;
                                    znomer[i]++;
                                    progbrt[i][k] = 1;
                                    prbrt[i][k] = 1;
                                    pr_idnt = 1;
                                }
                                else
                                    pr_idnt = 0;
                            }
                            else
                            {
                                if (pr_tras1 == 2)
                                {
                                    srv = arx->PelProgn - sgn->zpel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex1[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);
                                    spsrv = arx->sigmaPP4;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 1;
                                        progbrt[i][k] = 2;
                                        pr_idnt = 1;
                                    }
                                    else
                                        pr_idnt = 0;
                                }
                            }

                            if (pr_idnt == 0) // если не проидентифицировалась первая трасса
                            {
                                if (pr_tras2 == 2)
                                {
                                    srv = arx->PelPrognZ - sgn->zpel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex2[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);
                                    spsrv = arx->sigmaPP4Z;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 2;
                                        progbrt[i][k] = 2;
                                        pr_idnt = 1;
                                    }
                                    else
                                        pr_idnt = 0;
                                }
                                else
                                {
                                    if (pr_tras2 == 1)
                                    {
                                        srv = arx->PelPrognZ - sgn->pel;
                                        if (fabs(srv) > PISR)
                                        {
                                            pr_perex2[i][k] = 1;
                                            srv = PPI - fabs(srv);
                                        }
                                        else
                                            srv = fabs(srv);
                                        spsrv = arx->sigmaPP4Z;
                                        strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                        if (strob < STRB1 * MNJ)
                                            strob = STRB1 * MNJ;
                                        if (strob > STRB1 * MNJ)
                                            strob = STRB1 * MNJ;

                                        if (srv <= strob)
                                        {
                                            zmas[i][znomer[i]] = sgn->notm;
                                            znomer[i]++;
                                            prbrt[i][k] = 2;
                                            progbrt[i][k] = 1;
                                            pr_idnt = 1;
                                        }
                                        else
                                            pr_idnt = 0;
                                    }
                                }
                            }
                        }

                        else
                        {
                            if (arx->prbort == 1) // рассматриваем первую трассу
                            {
                                pr_perex1[i][k] = 0;
                                ky_tr1 = ky_pel(arx->PelProgn * 180. / M_PI, ykta->kyrs * 180. / M_PI);

                                if (ky_tr1 >= 0.)
                                    pr_tras1 = 1;
                                else
                                    pr_tras1 = 2;

                                if (pr_tras1 == 1)
                                {
                                    srv = arx->PelProgn - sgn->pel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex1[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);
                                    spsrv = arx->sigmaPP4;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel); // размер строба по пеленгу
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 1;
                                        progbrt[i][k] = 1;
                                    }
                                }
                                else
                                {
                                    srv = arx->PelProgn - sgn->zpel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex1[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);
                                    spsrv = arx->sigmaPP4;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 1;
                                        progbrt[i][k] = 2;
                                    }
                                }
                            }

                            else if (arx->prbort == 2) // рассматриваем вторую трассу
                            {
                                pr_perex2[i][k] = 0;
                                ky_tr2 = ky_pel(arx->PelPrognZ * 180. / M_PI, ykta->kyrs * 180. / M_PI);

                                if (ky_tr2 >= 0.)
                                    pr_tras2 = 1;
                                else
                                    pr_tras2 = 2;

                                if (pr_tras2 == 2)
                                {
                                    srv = arx->PelPrognZ - sgn->zpel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex2[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);
                                    spsrv = arx->sigmaPP4;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 2;
                                        progbrt[i][k] = 2;
                                    }
                                }
                                else
                                {
                                    srv = arx->PelPrognZ - sgn->pel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex2[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);
                                    spsrv = arx->sigmaPP4Z;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 2;
                                        progbrt[i][k] = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        prizn_kon_man = 2;
        for (k = 0; k < kols; k++, sgn++)
        {
            arx = parx;
            for (i = 0; i < adl; i++, arx++)
            {
                pr_strpel = 0;
                if (arx->anom != 0)
                {
                    if (sgn->pr_pchs_schs != 1)
                    {
                        if ((fabs(arx->f - sgn->f)) <= CHSLDC * cnst[arx->chd - 1])
                            pr_strpel = 1;
                    }
                    else
                        pr_strpel = 1;

                    if (pr_strpel == 1)
                    {
                        if (arx->prbort == 3)
                        {
                            pr_idnt = 0;
                            pr_perex1[i][k] = 0;
                            pr_perex2[i][k] = 0;
                            ky_tr1 = ky_pel(arx->PelProgn * 180. / M_PI, ykta->kyrs * 180. / M_PI);
                            ky_tr2 = ky_pel(arx->PelPrognZ * 180. / M_PI, ykta->kyrs * 180. / M_PI);
                            if (ky_tr1 >= 0)
                                pr_tras1 = 1;
                            else
                                pr_tras1 = 2;

                            if (ky_tr2 >= 0.)
                                pr_tras2 = 1;
                            else
                                pr_tras2 = 2;

                            if (pr_tras1 == 1)
                            {
                                srv = arx->PelProgn - sgn->pel;
                                if (fabs(srv) > PISR)
                                {
                                    pr_perex1[i][k] = 1;
                                    srv = PPI - fabs(srv);
                                }
                                else
                                    srv = fabs(srv);

                                spsrv = arx->sigmaPP4;
                                strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                if (strob < STRB1 * MNJ)
                                    strob = STRB1 * MNJ;
                                if (strob > STRB1 * MNJ)
                                    strob = STRB1 * MNJ;

                                if (srv <= strob)
                                {
                                    zmas[i][znomer[i]] = sgn->notm; // номер входного сигнала
                                    znomer[i]++;                    // номер сигнала в массиве
                                    prbrt[i][k] = 1;
                                    progbrt[i][k] = 1;

                                    srv = fabs(arx->PelPrognZ - sgn->zpel);
                                    spsrv = arx->sigmaPP4Z;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        pr_idnt = 1;
                                    }
                                    else
                                    {
                                        pr_idnt = 1;
                                        // nom_tr[i][k] = 1;
                                    }
                                }
                            }
                            else if (pr_tras1 == 2)
                            {
                                srv = arx->PelProgn - sgn->zpel;
                                if (fabs(srv) > PISR)
                                {
                                    pr_perex1[i][k] = 1;
                                    srv = PPI - fabs(srv);
                                }
                                else
                                    srv = fabs(srv);

                                spsrv = arx->sigmaPP4;
                                strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                if (strob < STRB1 * MNJ)
                                    strob = STRB1 * MNJ;
                                if (strob > STRB1 * MNJ)
                                    strob = STRB1 * MNJ;

                                if (srv <= strob)
                                {
                                    zmas[i][znomer[i]] = sgn->notm; // номер входного сигнала
                                    znomer[i]++;                    // номер сигнала в массиве
                                    prbrt[i][k] = 1;
                                    progbrt[i][k] = 2;
                                    pr_idnt = 1;
                                    // nom_tr[i][k] = 1;
                                }
                                else
                                    pr_idnt = 0;
                            }

                            if (pr_idnt == 0) // если первая трасса не проидентифицировалась
                            {
                                if (pr_tras2 == 2)
                                {
                                    srv = arx->PelPrognZ - sgn->zpel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex2[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);

                                    spsrv = arx->sigmaPP4Z;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm; // номер входного сигнала
                                        znomer[i]++;                    // номер сигнала в массиве
                                        prbrt[i][k] = 2;
                                        progbrt[i][k] = 2;
                                        pr_idnt = 1;
                                        // nom_tr[i][k] = 2;
                                    }
                                    else
                                        pr_idnt = 0;
                                }
                                else if (pr_tras2 == 1)
                                {
                                    srv = arx->PelPrognZ - sgn->pel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex2[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);

                                    spsrv = arx->sigmaPP4Z;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 2;
                                        progbrt[i][k] = 1;
                                        pr_idnt = 1;
                                        // nom_tr[i][k] = 2;
                                    }
                                    else
                                        pr_idnt = 0;
                                }
                            }
                        }
                        else // если трасса определена
                        {
                            if (arx->prbort == 1)
                            {

                                pr_perex1[i][k] = 0;
                                ky_tr1 = ky_pel(arx->PelProgn * 180. / M_PI, ykta->kyrs * 180. / M_PI);
                                if (ky_tr1 >= 0.)
                                    pr_tras1 = 1;
                                else
                                    pr_tras1 = 2;

                                if (pr_tras1 == 1)
                                {
                                    srv = arx->PelProgn - sgn->pel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex1[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);

                                    spsrv = arx->sigmaPP4;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel); // размер строба по пеленгу
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 1;
                                        progbrt[i][k] = 1;
                                    }
                                }
                                else if (pr_tras1 == 2)
                                {
                                    srv = arx->PelProgn - sgn->zpel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex1[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);

                                    spsrv = arx->sigmaPP4;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 1;
                                        progbrt[i][k] = 2;
                                    }
                                }
                            }
                            else if (arx->prbort == 2)
                            {

                                pr_perex2[i][k] = 0;
                                ky_tr2 = ky_pel(arx->PelPrognZ * 180. / M_PI, ykta->kyrs * 180. / M_PI);
                                if (ky_tr2 >= 0.)
                                    pr_tras2 = 1;
                                else
                                    pr_tras2 = 2;

                                if (pr_tras2 == 1)
                                {
                                    srv = arx->PelPrognZ - sgn->pel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex1[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);

                                    spsrv = arx->sigmaPP4Z;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (strob < STRB1 * MNJ)
                                        strob = STRB1 * MNJ;
                                    if (strob > STRB1 * MNJ)
                                        strob = STRB1 * MNJ;

                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 2;
                                        progbrt[i][k] = 1;
                                    }
                                }
                                else if (pr_tras2 == 2)
                                {
                                    srv = arx->PelPrognZ - sgn->zpel;
                                    if (fabs(srv) > PISR)
                                    {
                                        pr_perex2[i][k] = 1;
                                        srv = PPI - fabs(srv);
                                    }
                                    else
                                        srv = fabs(srv);

                                    spsrv = arx->sigmaPP4Z;
                                    strob = lz * sqrt(spsrv + sgn->cko_pel * sgn->cko_pel);
                                    if (srv <= strob)
                                    {
                                        zmas[i][znomer[i]] = sgn->notm;
                                        znomer[i]++;
                                        prbrt[i][k] = 2;
                                        progbrt[i][k] = 2;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (i = 0; i < adl; i++)
    {
        if (znomer[i] != 0) // если отметка проидентифицировалась
            for (j = 0; j < znomer[i]; j++)
            {
                mass_bn_a[zmas[i][j] - 1][noms[zmas[i][j] - 1]] = i + 1; // архивные номера проидентифицированных сигналов для каждой отметки
                noms[zmas[i][j] - 1]++;                                  // количество проидентифицированных архивных сигналов для каждой отметки
            }
    }
    for (i = 0, j = 0; i < adl; i++)
    {
        if (znomer[i] != 0) // если сигнал проидентифицировался хотя бы с одной омтеткой
        {
            imas[j] = i + 1; // архивные номера проидентифицированных сигналов
            j++;
        }
    }

    sgn = asgn;
    arx = parx;
    int chs = 0;
    int perem = 0;
    while (perem)
    {
        for (k = 0; k < kols; k++)
        {
            if (svn[k] == 0)
            {
                if (noms[k] == 0) // если не проидент. цели
                {
                    hhel[chs] = k + 1; // номера отметок,которые не проидентифицировались
                    chs++;
                    svn[k] = 1; // отметка обработана
                }
                else
                {
                    for (ll = 0, j = 0; j < noms[k]; j++) // кол-во сигналов для каждой отметки
                    {
                        min_srv = 0;
                        min_z = 0;
                        frst = 0;
                        if (siden[mass_bn_a[k][j] - 1] != 1)
                        {
                            for (i = 0; i < znomer[mass_bn_a[k][j] - 1]; i++)
                            {
                                if (svn[zmas[mass_bn_a[k][j] - 1][i] - 1] != 1)
                                {
                                    if (prbrt[mass_bn_a[k][j] - 1][zmas[mass_bn_a[k][j] - 1][i] - 1] == 1)
                                    {
                                        spsrv = (arx + (mass_bn_a[k][j] - 1))->PelProgn; // прогнозируемый пеленг
                                    }
                                    else
                                    {
                                        spsrv = (arx + (mass_bn_a[k][j] - 1))->PelPrognZ;
                                    }

                                    if (progbrt[mass_bn_a[k][j] - 1][zmas[mass_bn_a[k][j] - 1][i] - 1] == 1)
                                    {
                                        sspsrv = (sgn + (zmas[mass_bn_a[k][j] - 1][i] - 1))->pel; // измеренный пеленг
                                    }
                                    else
                                    {
                                        sspsrv = (sgn + (zmas[mass_bn_a[k][j] - 1][i] - 1))->zpel;
                                    }

                                    srv = fabs(spsrv - sspsrv);
                                    if (frst == 0)
                                    {
                                        frst = 1;
                                        min_srv = srv;
                                        min_z = zmas[mass_bn_a[k][j] - 1][i];
                                    }
                                    else
                                    {
                                        if (srv < min_srv)
                                        {
                                            min_srv = srv;
                                            min_z = zmas[mass_bn_a[k][j] - 1][i];
                                        }
                                    }
                                } // конец обработки отметки
                            }
                        }
                        // сраниваем ближайшую отметку с текущей
                        if (min_z == k + 1)
                        // если совпадает
                        {
                            hnm[ll] = mass_bn_a[k][j]; // в массив проидентифицированных сигналов для отметки
                            ll++;                      // кол-во сигналов, для которых отметка минимальная
                        }
                    }
                    if (ll != 0) // если прориденифицировалось несколько сигналов
                    {
                        if (ll == 1) // если отметка минимальна для одного сигнала
                        {
                            svn[k] = 1;
                            siden[hnm[0] - 1] == 1;
                            iden[hnm[0] - 1] = k + 1; // номер отметки для данного сигнала
                        }
                        else // выбор архивного сигнала с минимальным расстоянием по пеленгу
                        {
                            frst = 0;
                            min_srv = 0;
                            min_z = 0;
                            for (i = 0; i < ll; i++)
                            {
                                if (prbrt[hnm[i] - 1][k] == 1)
                                {
                                    spsrv = (arx + (hnm[i] - 1))->PelProgn;
                                }
                                else
                                {
                                    spsrv = (arx + (hnm[i] - 1))->PelPrognZ;
                                }

                                if (progbrt[hnm[i] - 1][k] == 1)
                                {
                                    sspsrv = (sgn + k)->pel;
                                }
                                else
                                {
                                    sspsrv = (sgn + k)->zpel;
                                }
                                srv = fabs(spsrv - sspsrv);
                                if (frst == 0)
                                {
                                    frst = 1;
                                    min_z = hnm[i];
                                    min_srv = srv;
                                }

                                if (srv < min_srv)
                                {
                                    min_srv = srv;
                                    min_z = hnm[i];
                                }
                            }
                            iden[min_z - 1] = k + 1;
                            svn[k] = 1;
                            siden[min_z - 1] = 1;
                        }
                    }
                }
            }
        }
    }
}
// функция вычисления курсового угла
float ky_pel(float pel, float kurs)
{
    float ky;

    if ((pel - kurs) > 0.)
        if ((pel - kurs) < 180.)
            ky = pel - kurs;
        else
            ky = pel - kurs - 360.;
    else if ((kurs - pel) < 180.)
        ky = pel - kurs;
    else
        ky = 360. - (kurs - pel);
    return (ky);
}
