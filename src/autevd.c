#include <math.h>
#include <stdlib.h>

#include "functions.h"
#include "autevd.h"
#include "autlim.h"
#include "auto_c.h"

#include "auto_nox.h"
#include "integers.h"
#include "x_auto.h"

#define SPER 3
#define UPER 4
#define SEQ 1
#define UEQ 2

XAuto x_auto;

int32 DiagFlag = 0;

static struct {
    int32 pt;
    int32 br;
    double evr[NAUTO], evi[NAUTO];
} my_ev;

void
init_auto(int32 ndim, int32 nicp, int32 ips, int32 irs, int32 ilp, int32 ntst,
          int32 isp, int32 isw, int32 nmx, int32 npr, double ds, double dsmin,
          double dsmax, double rl0, double rl1, double a0, double a1, int32 ip1,
          int32 ip2, int32 ip3, int32 ip4, int32 ip5, double epsl, double epsu,
          double epss, int32 ncol) {
    /* here are the constants that we do not allow the user to change */
    int32 nnbc;
    int32 i;
    x_auto.iad = aauto.iad;
    x_auto.iplt = 0;
    x_auto.mxbf = aauto.mxbf;
    x_auto.iid = aauto.iid;
    x_auto.itmx = aauto.itmx;
    x_auto.itnw = aauto.itnw;
    x_auto.nwtn = aauto.nwtn;
    x_auto.jac = 0;
    x_auto.iads = aauto.iads;
    x_auto.nthl = 1;
    x_auto.ithl[0] = 10;
    x_auto.thl[0] = 0.0;
    x_auto.nint = 0;

    if (ips == 4)
        nnbc = ndim;
    else
        nnbc = 0;
    x_auto.ndim = ndim;
    x_auto.nbc = nnbc;
    x_auto.ips = ips;
    x_auto.irs = irs;
    x_auto.ilp = ilp;
    x_auto.nicp = nicp;
    x_auto.icp[0] = ip1;
    x_auto.icp[1] = ip2;
    x_auto.icp[2] = ip3;
    x_auto.icp[3] = ip4;
    x_auto.icp[4] = ip5;
    x_auto.ntst = ntst;
    x_auto.ncol = ncol;
    x_auto.isp = isp;
    x_auto.isw = isw;
    x_auto.nmx = nmx;
    x_auto.rl0 = rl0;
    x_auto.rl1 = rl1;
    x_auto.a0 = a0;
    x_auto.a1 = a1;
    x_auto.npr = npr;

    x_auto.epsl = epsl;
    x_auto.epss = epss;
    x_auto.epsu = epsu;
    x_auto.ds = ds;
    x_auto.dsmax = dsmax;
    x_auto.dsmin = dsmin;

    x_auto.nuzr = NAutoUzr;
    for (i = 0; i < NAutoUzr; i++) {
        x_auto.iuz[i] = (int32)UzrPar[i];
        x_auto.vuz[i] = outperiod[i];
    }
    return;
}

void
send_eigen(int32 ibr, int32 ntot, int32 n, doublecomplex *ev) {
    int32 i;
    double er, cs, sn;
    my_ev.pt = abs(ntot);
    my_ev.br = abs(ibr);
    for (i = 0; i < n; i++) {
        er = exp((ev + i)->r);
        cs = cos((ev + i)->i);
        sn = sin((ev + i)->i);
        my_ev.evr[i] = er*cs;
        my_ev.evi[i] = er*sn;
    }
    return;
}

void
send_mult(int32 ibr, int32 ntot, int32 n, doublecomplex *ev) {
    int32 i;
    my_ev.pt = abs(ntot);
    my_ev.br = abs(ibr);
    for (i = 0; i < n; i++) {
        my_ev.evr[i] = (ev + i)->r;
        my_ev.evi[i] = (ev + i)->i;
    }
    return;
}

/* Only unit 8,3 or q.prb is important; all others are unnecesary */

int32
get_bif_type(int32 ibr, int32 ntot) {
    int32 type = SEQ;

    if (ibr < 0 && ntot < 0)
        type = SPER;
    if (ibr < 0 && ntot > 0)
        type = UPER;
    if (ibr > 0 && ntot > 0)
        type = UEQ;
    if (ibr > 0 && ntot < 0)
        type = SEQ;
    return type;
}

void
addbif(iap_type *iap, int64 ntots, int64 ibrs, double *par, int64 *icp,
       int32 lab, double *a, double *uhigh, double *ulow, double *u0,
       double *ubar) {
    int32 type;
    int32 icp1 = (int32)icp[0], icp2 = (int32)icp[1], icp3 = (int32)icp[2],
          icp4 = (int32)icp[3];
    double per = par[10];
    type = get_bif_type((int32)ibrs, (int32)ntots);

    /*if(my_ev.br==abs(*ibr)&&my_ev.pt==abs(*ntot)){evflag=1;}*/
    if (iap->ntot == 1) {
        add_point(par, per, uhigh, ulow, ubar, *a, type, 0, lab, icp1, icp2,
                  AutoTwoParam, my_ev.evr, my_ev.evi);
    } else {
        add_point(par, per, uhigh, ulow, ubar, *a, type, 1, lab, icp1, icp2,
                  AutoTwoParam, my_ev.evr, my_ev.evi);
    }

    if (DiagFlag == 0) {
        /* diagram_start(*ndim); */
        edit_start((int32)ibrs, (int32)ntots, (int32)iap->itp, lab,
                   (int32)iap->nfpr, *a, uhigh, ulow, u0, ubar, par, per,
                   (int32)iap->ndim, icp1, icp2, icp3, icp4, my_ev.evr,
                   my_ev.evi);
        DiagFlag = 1;
        return;
    }
    add_diagram((int32)ibrs, (int32)ntots, (int32)iap->itp, lab,
                (int32)iap->nfpr, *a, uhigh, ulow, u0, ubar, par, per,
                (int32)iap->ndim, icp1, icp2, icp3, icp4, AutoTwoParam,
                my_ev.evr, my_ev.evi);
    return;
}

int32
eigrf_(double *a, int32 *n, doublecomplex *ecv, double *work, int32 *ier) {
    double ev[400];
    gear_eigen(*n, a, ev, work, ier);
    for (int32 i = 0; i < *n; i++) {
        (ecv + i)->r = ev[2*i];
        (ecv + i)->i = ev[2*i + 1];
    }
    return 0;
}
