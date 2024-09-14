#ifndef phsplan_h
#define phsplan_h
#include "integers.h"

/*          This include file has all of the global phaseplane
            stuff.
            This is not where it is defined

*/

#include "xpplim.h"

extern int32 NEQ, NODE;

extern int32 PLOT_3D;

extern float **storage;
extern int32 storind, MAXSTOR, INFLAG, MY_STOR, STORFLAG;

extern int32 FOREVER;

extern int32 ENDSING, PAUSER;

/*  extern GRAPH *MyGraph; */

extern int32 METHOD, NJMP;
extern double HMIN, HMAX, TOLER, ATOLER, BOUND, DELAY;
extern double NULL_ERR, EVEC_ERR, NEWT_ERR;
extern int32 EVEC_ITER, NMESH, NC_ITER;

extern float *fft_data;
extern int32 FFT;

extern float *hist_data;
extern int32 HIST, HVAR, hist_ind;

extern double TEND, DELTA_T, T0, TRANS;

extern double *WORK;
extern int32 IWORK[1000];

extern int32 TORUS, itor[MAX_ODE];
extern double TOR_PERIOD;

extern int32 POIMAP, POISGN, POIEXT, SOS, POIVAR;
extern double POIPLN;

extern int32 NULL_HERE;
extern float *X_n, *Y_n;

extern char uvar_names[MAX_ODE][12];

#endif
