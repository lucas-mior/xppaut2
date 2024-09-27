#include "autlim.h"
#include "auto_nox.h"
#include "functions.h"
#include "parserslow.h"
#include "integers.h"
#include "x_auto.h"
#include "auto_c.h"

/*    Hooks to xpp RHS     */

int32
func(int64 ndim, double *u, int64 *icp, double *par, int64 ijac, double *f,
     double *dfdu, double *dfdp) {
    double zz[NAUTO];
    double y[NAUTO];
    double yp[NAUTO];
    double xp[NAUTO];

    (void)icp;
    (void)dfdu;
    (void)dfdp;

    for (int32 i = 0; i < NAutoPar; i++) {
        constants[Auto_index_to_array[i]] = par[i];
    }
    derived_evaluate();
    tabular_redo_all_fun_tables();
    rhs_function(0.0, u, f, (int32)ndim);
    if (ijac == 1) {
        gear_jac_trans(u, y, yp, xp, NEWT_ERR, dfdu, (int32)ndim);
    }
    if (METHOD > 0 || NJMP == 1) {
        return 0;
    }
    for (int32 i = 1; i < NJMP; i++) {
        for (int32 j = 0; j < ndim; j++) {
            zz[j] = f[j];
        }
        rhs_function(0.0, zz, f, (int32)ndim);
    }

    return 0;
}

int32
stpnt(int64 ndim, double t, double *u, double *par) {
    double p;

    for (int32 i = 0; i < NAutoPar; i++) {
        par[i] = constants[Auto_index_to_array[i]];
    }

    if (NewPeriodFlag == 0) {
        for (int32 i = 0; i < ndim; i++) {
            u[i] = last_ic[i];
        }
        return 0;
    }

    auto_nox_get_start_period(&p);
    par[10] = p;
    if (HomoFlag != 1) {
        auto_nox_get_start_orbit(u, t, (int32)ndim);
    }
    /*  printf("%d %d %g %g %g %g \n",ndim,HomoFlag,t,u[0],u[1],p); */
    if (HomoFlag == 1) {
        auto_nox_get_shifted_orbit(u, t, p, (int32)ndim);
        for (int32 i = 0; i < ndim; i++) {
            par[11 + i] = homo_l[i];
        }
    }
    if (HomoFlag == 2) { /* heteroclinic */
        for (int32 i = 0; i < ndim; i++) {
            par[11 + i] = homo_l[i];
            par[11 + i + ndim] = homo_r[i];
        }
    }
    return 0;
}

int32
bcnd(int64 ndim, double *par, int64 *icp, int64 nbc, double *u0, double *u1,
     int64 ijac, double *fb, double *dbc) {
    (void)dbc;
    (void)ijac;
    (void)icp;
    (void)ndim;
    /* Hooks to the XPP bc parser!! */

    for (int32 i = 0; i < NAutoPar; i++) {
        constants[Auto_index_to_array[i]] = par[i];
    }

    derived_evaluate();
    tabular_redo_all_fun_tables();
    pp_shoot_do_bc(u0, 0.0, u1, 1.0, fb, (int32)nbc);

    return 0;
}

int32
icnd(int64 ndim, double *par, int64 *icp, int64 nint, double *u, double *uold,
     double *udot, double *upold, int64 ijac, double *fi, double *dint) {
    (void)dint;
    (void)ijac;
    (void)fi;
    (void)upold;
    (void)udot;
    (void)uold;
    (void)u;
    (void)nint;
    (void)icp;
    (void)par;
    (void)ndim;
    return 0;
}

int32
fopt(int64 ndim, double *u, int64 *icp, double *par, int64 ijac, double *fs,
     double *dfdu, double *dfdp) {
    (void)dfdp;
    (void)dfdu;
    (void)fs;
    (void)ijac;
    (void)par;
    (void)icp;
    (void)u;
    (void)ndim;
    return 0;
}

/* Not sure what to do here; I think  do nothing  since IEQUIB is always -2 */
int32
pvls(int64 ndim, double *u, double *par) {
    (void)par;
    (void)u;
    (void)ndim;
    return 0;
}
