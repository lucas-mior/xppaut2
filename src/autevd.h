#ifndef autevd_h_
#define autevd_h_

#include "autlim.h"
#include "auto_f2c.h"
#include "auto_c.h"
#include "integers.h"

typedef struct {
    int32 pt, br;
    double evr[NAUTO], evi[NAUTO];
} EIGVAL;

void send_eigen(int32 ibr, int32 ntot, int32 n, doublecomplex *ev);
void send_mult(int32 ibr, int32 ntot, int32 n, doublecomplex *ev);
int32 get_bif_type(int32 ibr, int32 ntot);
void addbif(iap_type *iap, int64 ntots, int64 ibrs, double *par, int64 *icp,
            int32 labw, double *a, double *uhigh, double *ulow, double *u0,
            double *ubar);
int32 eigrf_(double *a, int32 *n, doublecomplex *ecv, double *work, int32 *ier);
void init_auto(int32 ndim, int32 nicp, int32 ips, int32 irs, int32 ilp,
               int32 ntst, int32 isp, int32 isw, int32 nmx, int32 npr,
               double ds, double dsmin, double dsmax, double rl0, double rl1,
               double a0, double a1, int32 ip1, int32 ip2, int32 ip3, int32 ip4,
               int32 ip5, double epsl, double epsu, double epss, int32 ncol);

#endif
