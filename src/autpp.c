#include "autlim.h"
#include "auto_nox.h"
#include "functions.h"
#include "integers.h"
#include "x_auto.h"

void getjactrans(double *x, double *y, double *yp, double *xp, double eps,
                 double *d, int32 n);
extern XAuto x_auto;

/*    Hooks to xpp RHS     */
extern int32 (*rhs)(double t, double *y, double *ydot, int32 neq);
extern double constants[], last_ic[];

extern int32 Auto_index_to_array[8];
extern int32 NewPeriodFlag;
extern int32 AutoTwoParam, NAutoPar;
extern int32 HomoFlag;
extern double homo_l[100], homo_r[100];
extern int32 METHOD, NJMP;
extern double outperiod[];
extern int32 UzrPar[], NAutoUzr;

extern double NEWT_ERR;
int32
func(int64 ndim, double *u, int64 *icp, double *par, int64 ijac, double *f,
     double *dfdu, double *dfdp) {
    int32 i, j;
    double zz[NAUTO];
    double y[NAUTO], yp[NAUTO], xp[NAUTO];
    for (i = 0; i < NAutoPar; i++) {
        constants[Auto_index_to_array[i]] = par[i];
    }
    evaluate_derived();
    redo_all_fun_tables();
    rhs(0.0, u, f, ndim);
    if (ijac == 1) {
        getjactrans(u, y, yp, xp, NEWT_ERR, dfdu, ndim);
    }
    if (METHOD > 0 || NJMP == 1)
        return 0;
    for (i = 1; i < NJMP; i++) {
        for (j = 0; j < ndim; j++)
            zz[j] = f[j];
        rhs(0.0, zz, f, ndim);
    }

    return 0;

} /* func_ */

int32
stpnt(int64 ndim, double t, double *u, double *par) {
    int32 i;

    double p;

    for (i = 0; i < NAutoPar; i++)
        par[i] = constants[Auto_index_to_array[i]];

    if (NewPeriodFlag == 0) {
        for (i = 0; i < ndim; i++)
            u[i] = last_ic[i];
        return 0;
    }

    get_start_period(&p);
    par[10] = p;
    if (HomoFlag != 1)
        get_start_orbit(u, t, ndim);
    /*  printf("%d %d %g %g %g %g \n",ndim,HomoFlag,t,u[0],u[1],p); */
    if (HomoFlag == 1) {

        get_shifted_orbit(u, t, p, ndim);
        for (i = 0; i < ndim; i++) {
            par[11 + i] = homo_l[i];
        }
    }
    if (HomoFlag == 2) { /* heteroclinic */
        for (i = 0; i < ndim; i++) {
            par[11 + i] = homo_l[i];
            par[11 + i + ndim] = homo_r[i];
        }
    }
    return 0;

} /* stpnt_ */

int32
bcnd(int64 ndim, double *par, int64 *icp, int64 nbc, double *u0, double *u1,
     int64 ijac, double *fb, double *dbc) {
    int32 i;
    /* Hooks to the XPP bc parser!! */

    for (i = 0; i < NAutoPar; i++) {
        constants[Auto_index_to_array[i]] = par[i];
    }

    evaluate_derived();
    redo_all_fun_tables();
    do_bc(u0, 0.0, u1, 1.0, fb, nbc);

    return 0;
} /* bcnd_ */

int32
icnd(int64 ndim, double *par, int64 *icp, int64 *nint, double *u, double *uold,
     double *udot, double *upold, double *fi, int64 *ijac, double *dint) {
    int32 i;
    double dum = 0.0;
    /*
   for(i=0;i<Homo_n;i++)
     dum+=upold[i]*(u[i]-uold[i]);
   fi[0]=dum;
    */
    return 0;
} /* icnd_ */

int32
fopt(int64 *ndim, double *u, int64 *icp, double *par, int64 *ijac, double *fs,
     double *dfdu, double *dfdp) {
    /*     ---------- ---- */
    return 0;
} /* fopt_ */

/*  Not sure what to do here; I think  do nothing  since IEQUIB is always
    -2
*/
int32
pvls(int64 ndim, const double *u, double *par) {
    return 0;
}
