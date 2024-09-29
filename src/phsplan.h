#ifndef PHSPLAN_H
#define PHSPLAN_H
#include "integers.h"

/* This include file has all of the global phaseplane stuff.
 * This is not where it is defined */

#include "xpplim.h"

extern int32 NEQ;
extern int32 NODE;

extern int32 PLOT_3D;

extern double **storage;
extern int32 storind;
extern int32 MAXSTOR;
extern int32 INFLAG;
extern int32 MY_STOR;
extern int32 STORFLAG;

extern int32 forever;

extern int32 end_sing;

extern int32 METHOD;
extern int32 NJMP;
extern double HMIN;
extern double HMAX;
extern double TOLER;
extern double atoler;
extern double bound;
extern double delay;
extern double NULL_ERR;
extern double evec_err;
extern double NEWT_ERR;
extern int32 evec_iter;
extern int32 NMESH;
extern int32 NC_ITER;

extern double *fft_data;
extern int32 FFT;

extern double *hist_data;
extern int32 HIST;
extern int32 HVAR;
extern int32 hist_ind;

extern double TEND;
extern double delta_t;
extern double T0;
extern double TRANS;

extern double *WORK;

extern int32 TORUS;
extern int32 itor[MAX_ODE];
extern double TOR_PERIOD;

extern int32 POIMAP;
extern int32 POISGN;
extern int32 POIEXT;
extern int32 SOS;
extern int32 POIVAR;
extern double POIPLN;

extern int32 NULL_HERE;

#endif
