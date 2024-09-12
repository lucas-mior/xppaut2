#ifndef _autevd_h_
#define _autevd_h_
#include "integers.h"

#include "autlim.h"
#include "auto_f2c.h"
#include "auto_c.h"
/* autevd.c */

/* typedef struct {double r,i;} doublecomplex; */

typedef struct {
    int32 pt, br;
    double evr[NAUTO], evi[NAUTO];
} EIGVAL;

void send_eigen(int32 ibr, int32 ntot, int32 n, doublecomplex *ev);
void send_mult(int32 ibr, int32 ntot, int32 n, doublecomplex *ev);
int32 get_bif_type(int32 ibr, int32 ntot, int32 lab);
void addbif(iap_type *iap, rap_type *rap, int64 ntots, int64 ibrs, double *par,
            int64 *icp, int32 labw, double *a, double *uhigh, double *ulow,
            double *u0, double *ubar);
double etime_(double *z);
int32 eigrf_(double *a, int32 *n, int32 *m, doublecomplex *ecv, double *work,
           int32 *ier);
void init_auto(int32 ndim, int32 nicp, int32 nbc, int32 ips, int32 irs, int32 ilp, int32 ntst,
               int32 isp, int32 isw, int32 nmx, int32 npr, double ds, double dsmin,
               double dsmax, double rl0, double rl1, double a0, double a1,
               int32 ip1, int32 ip2, int32 ip3, int32 ip4, int32 ip5, int32 nuzr,
               double epsl, double epsu, double epss, int32 ncol);

#endif
