/* autlib3.f -- translated by f2c (version 19970805).
 * You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

#include "auto_f2c.h"
#include "auto_c.h"
#include "integers.h"
#include "autlib.h"
#include "xmalloc.h"

/* ----------------------------------------------------------------------- */
/*  Subroutines for the Continuation of Folds (Algebraic Problems) */
/* ----------------------------------------------------------------------- */

int32
fnlp(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double rtmp;
    double ep;
    int64 ndm;
    double umx;

    // Generates the equations for the 2-par continuation of folds.

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;

    // Generate the function.

    fflp(iap, rap, ndim, u, uold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        fflp(iap, rap, ndim, global_scratch.uu1, uold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        fflp(iap, rap, ndim, global_scratch.uu2, uold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    par[icp[0]] += ep;

    fflp(iap, rap, ndim, u, uold, icp, par, global_scratch.ff1, ndm, global_scratch.dfu,
         global_scratch.dfp);

    for (int64 j = 0; j < ndim; ++j) {
        ARRAY2D(dfdp, j, (icp[0])) = (global_scratch.ff1[j] - f[j]) / ep;
    }

    par[icp[0]] -= ep;

    return 0;
}

int32
fflp(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 ips;

    // Parameter adjustments
    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    ips = iap->ips;

    par[icp[1]] = u[-1 + ndim];
    if (ips == -1) {
        fnds(iap, rap, ndm, u, uold, icp, par, 1, f, dfdu, dfdp);
    } else {
        funi(iap, rap, ndm, u, uold, icp, par, 1, f, dfdu, dfdp);
    }

    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = 0.;
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm + i] += ARRAY2D(dfdu, i, j)*u[ndm + j];
        }
    }

    f[-1 + ndim] = -1.;

    for (int64 i = 0; i < ndm; ++i) {
        f[-1 + ndim] += u[ndm + i]*u[ndm + i];
    }

    return 0;
}

int32
stpnlp(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *u) {
    int64 ndim;

    double uold;
    int64 nfpr1;
    double *f;
    double *v;
    int64 found;

    int64 ndm;
    int64 ips;
    int64 irs;

    f = xmalloc(sizeof(*f)*(usize)(iap->ndim));
    v = xmalloc(sizeof(*v)*(usize)(iap->ndim));
    // Generates starting data for the continuation of folds.

    // Local

    // Parameter adjustments

    ndim = iap->ndim;
    ips = iap->ips;
    irs = iap->irs;
    ndm = iap->ndm;

    findlb(iap, rap, irs, &nfpr1, &found);
    readlb(u, par);

    if (ips == -1) {
        fnds(iap, rap, ndm, u, &uold, icp, par, 1, f, global_scratch.dfu, global_scratch.dfp);
    } else {
        funi(iap, rap, ndm, u, &uold, icp, par, 1, f, global_scratch.dfu, global_scratch.dfp);
    }
    nlvc(ndm, ndm, 1, global_scratch.dfu, v);
    nrmlz(&ndm, v);
    for (int64 i = 0; i < ndm; ++i) {
        u[ndm + i] = v[i];
    }
    u[-1 + ndim] = par[icp[1]];

    free(f);
    free(v);
    return 0;
}

/* ----------------------------------------------------------------------- */
/*     Subroutines for the Optimization of Algebraic Systems */
/* ----------------------------------------------------------------------- */

int32
fnc1(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double ddp[NPARX], *ddu;
    int64 ndm;

    ddu = xmalloc(sizeof(*ddu)*(usize)(iap->ndim));
    // Generate the equations for the continuation scheme used for
    // the optimization of algebraic systems (one parameter).

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;

    par[icp[1]] = u[-1 + ndim];
    funi(iap, rap, ndm, u, uold, icp, par, ijac, f, dfdu, dfdp);

    // Rearrange (Since dimensions in FNC1 and FUNI differ).

    if (ijac != 0) {
        for (int64 j = ndm - 1; j >= 0; --j) {
            for (int64 i = ndm - 1; i >= 0; --i) {
                ARRAY2D(dfdu, i, j) = dfdu[j*ndm + i];
            }
        }

        for (int64 j = NPARX - 1; j >= 0; --j) {
            for (int64 i = ndm - 1; i >= 0; --i) {
                ARRAY2D(dfdp, i, j) = dfdp[j*ndm + i];
            }
        }
    }

    fopi(iap, rap, ndm, u, icp, par, ijac, &f[-1 + ndim], ddu, ddp);
    f[-1 + ndim] = par[icp[0]] - f[-1 + ndim];

    if (ijac != 0) {
        for (int64 i = 0; i < ndm; ++i) {
            ARRAY2D(dfdu, (ndim - 1), i) = -ddu[i];
            ARRAY2D(dfdu, i, (ndim - 1)) = ARRAY2D(dfdp, i, (icp[1]));
            ARRAY2D(dfdp, i, (icp[0])) = 0.;
        }
        ARRAY2D(dfdu, (ndim - 1), (ndim - 1)) = -ddp[icp[1]];
        ARRAY2D(dfdp, (ndim - 1), (icp[0])) = 1.;
    }
    free(ddu);
    return 0;
}

int32
stpnc1(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *u) {
    int64 ndim;

    int64 nfpr;

    int64 ndm;
    double fop;
    double dum;

    // Generate starting data for optimization problems (one parameter).

    // Parameter adjustments

    ndim = iap->ndim;
    ndm = iap->ndm;

    stpnt(ndim, 0.0, u, par);
    nfpr = 2;
    iap->nfpr = nfpr;
    fopi(iap, rap, ndm, u, icp, par, 0, &fop, &dum, &dum);
    par[icp[0]] = fop;
    u[-1 + ndim] = par[icp[1]];

    return 0;
}

int32
fnc2(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double rtmp;
    double ep;
    int64 ndm;
    double umx;

    // Generate the equations for the continuation scheme used for the
    // optimization of algebraic systems (more than one parameter).

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;

    // Generate the function.

    ffc2(iap, rap, ndim, u, uold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        ffc2(iap, rap, ndim, global_scratch.uu1, uold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        ffc2(iap, rap, ndim, global_scratch.uu2, uold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < ndim; ++i) {
        ARRAY2D(dfdp, i, (icp[0])) = 0.;
    }
    ARRAY2D(dfdp, (ndim - 1), (icp[0])) = 1.;

    return 0;
}

int32
ffc2(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 icpm;

    int64 nfpr;
    double ddp[NPARX], *ddu, fop;
    int64 ndm2;

    ddu = xmalloc(sizeof(*ddu)*(usize)(iap->ndim));
    // Local

    // Parameter adjustments
    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    nfpr = iap->nfpr;

    for (int64 i = 1; i < nfpr; ++i) {
        par[icp[i]] = u[(ndm*2) + i];
    }
    funi(iap, rap, ndm, u, uold, icp, par, 2, f, dfdu, dfdp);
    fopi(iap, rap, ndm, u, icp, par, 2, &fop, ddu, ddp);

    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = ddu[i]*u[(ndm*2)];
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm + i] += ARRAY2D(dfdu, j, i)*u[ndm + j];
        }
    }

    ndm2 = ndm*2;
    icpm = nfpr - 2;
    for (int64 i = 0; i < icpm; ++i) {
        f[ndm2 + i] = ddp[icp[i + 1]]*u[ndm2];
    }

    for (int64 i = 0; i < icpm; ++i) {
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm2 + i] += u[ndm + j]*ARRAY2D(dfdp, j, (icp[i + 1]));
        }
    }

    f[ndim - 2] = u[ndm2]*u[ndm2] - 1;
    for (int64 j = 0; j < ndm; ++j) {
        f[ndim - 2] += u[ndm + j]*u[ndm + j];
    }
    f[-1 + ndim] = par[icp[0]] - fop;

    free(ddu);
    return 0;
}

int32
stpnc2(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *u) {
    int64 ndim;

    double uold;
    int64 nfpr;
    double *f;
    double *v;
    int64 found;

    double *dd;
    double dp[NPARX], *du;

    int64 ndm;
    double fop;
    int64 irs;

    f = xmalloc(sizeof(*f)*(usize)(iap->ndim));
    v = xmalloc(sizeof(*v)*(usize)(iap->ndim));
    dd = xmalloc(sizeof(*dd)*(usize)((iap->ndim)*(iap->ndim)));
    du = xmalloc(sizeof(*du)*(usize)(iap->ndim));
    // Generates starting data for the continuation equations for
    // optimization of algebraic systems (More than one parameter).

    // Local

    // Parameter adjustments
    ndim = iap->ndim;
    irs = iap->irs;
    ndm = iap->ndm;

    findlb(iap, rap, irs, &nfpr, &found);
    ++nfpr;
    iap->nfpr = nfpr;
    readlb(u, par);

    if (nfpr == 3) {
        funi(iap, rap, ndm, u, &uold, icp, par, 2, f, global_scratch.dfu, global_scratch.dfp);
        fopi(iap, rap, ndm, u, icp, par, 2, &fop, du, dp);
        //       TRANSPOSE
        for (int64 i = 0; i < ndm; ++i) {
            for (int64 j = 0; j < ndm; ++j) {
                dd[i + j*ndim] = global_scratch.dfu[i*ndm + j];
            }
        }
        for (int64 i = 0; i < ndm; ++i) {
            dd[i + ndm*ndim] = du[i];
            dd[ndm + i*ndim] = global_scratch.dfp[(icp[1])*ndm + i];
        }
        dd[ndm + ndm*ndim] = dp[icp[1]];
        nlvc(ndm + 1, ndim, 1, dd, v);
        {
            int64 tmp = ndm + 1;
            nrmlz(&tmp, v);
        }
        for (int64 i = 0; i < ndm + 1; ++i) {
            u[ndm + i] = v[i];
        }
        par[icp[0]] = fop;
    }

    for (int64 i = 0; i < nfpr - 1; ++i) {
        u[ndim - nfpr + 1 + i] = par[icp[i + 1]];
    }

    free(f);
    free(v);
    free(dd);
    free(du);
    return 0;
}

/* ----------------------------------------------------------------------- */
/*        Subroutines for Discrete Dynamical Systems */
/* ----------------------------------------------------------------------- */

int32
fnds(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    // Generate the equations for continuing fixed points.

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    funi(iap, rap, ndim, u, uold, icp, par, ijac, f, dfdu, dfdp);

    for (int64 i = 0; i < ndim; ++i) {
        f[i] -= u[i];
    }

    if (ijac == 0) {
        return 0;
    }

    for (int64 i = 0; i < ndim; ++i) {
        --ARRAY2D(dfdu, i, i);
    }

    return 0;
}

/* ----------------------------------------------------------------------- */
/*        Subroutines for Time Integration of ODEs */
/* ----------------------------------------------------------------------- */

int32
fnti(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double told;
    double dt;

    // Generate the equations for continuing fixed points.

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    funi(iap, rap, ndim, u, uold, icp, par, ijac, f, dfdu, dfdp);

    told = rap->tivp;
    dt = par[icp[0]] - told;

    for (int64 i = 0; i < ndim; ++i) {
        ARRAY2D(dfdp, i, (icp[0])) = f[i];
        f[i] = dt*f[i] - u[i] + uold[i];
    }

    if (ijac == 0) {
        return 0;
    }

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, i, j) = dt*ARRAY2D(dfdu, i, j);
        }
        ARRAY2D(dfdu, i, i) += -1.;
    }

    return 0;
}

/* ----------------------------------------------------------------------- */
/*     Subroutines for the Continuation of Hopf Bifurcation Points (Maps) */
/* ----------------------------------------------------------------------- */

int32
fnhd(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double rtmp;
    double ep;
    int64 ndm;
    double umx;

    // Generates the equations for the 2-parameter continuation of Hopf
    // bifurcation points for maps.

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;

    // Generate the function.

    ffhd(iap, rap, ndim, u, uold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        ffhd(iap, rap, ndim, global_scratch.uu1, uold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        ffhd(iap, rap, ndim, global_scratch.uu2, uold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    par[icp[0]] += ep;

    ffhd(iap, rap, ndim, u, uold, icp, par, global_scratch.ff1, ndm, global_scratch.dfu,
         global_scratch.dfp);

    for (int64 j = 0; j < ndim; ++j) {
        ARRAY2D(dfdp, j, icp[0]) = (global_scratch.ff1[j] - f[j]) / ep;
    }

    par[icp[0]] -= ep;

    return 0;
}

int32
ffhd(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double thta;

    double c1;
    double s1;
    int64 ndm2;

    // Parameter adjustments

    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    ndm2 = ndm*2;

    thta = u[-1 + ndim - 1];
    s1 = sin(thta);
    c1 = cos(thta);
    par[icp[1]] = u[-1 + ndim];
    funi(iap, rap, ndm, u, uold, icp, par, 1, f, dfdu, dfdp);
    for (int64 i = 0; i < ndm; ++i) {
        f[i] -= u[i];
        ARRAY2D(dfdu, i, i) -= c1;
    }

    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = s1*u[ndm2 + i];
        f[ndm2 + i] = -s1*u[ndm + i];
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm + i] += ARRAY2D(dfdu, i, j)*u[ndm + j];
            f[ndm2 + i] += ARRAY2D(dfdu, i, j)*u[ndm2 + j];
        }
    }

    f[ndim - 2] = -1.;

    for (int64 i = 0; i < ndm; ++i) {
        f[ndim - 2] = f[ndim - 2] + u[ndm + i]*u[ndm + i] + u[ndm2 + i]*u[ndm2 + i];
    }

    f[-1 + ndim] = 0.;

    for (int64 i = 0; i < ndm; ++i) {
        f[-1 + ndim] = f[-1 + ndim] + uold[ndm2 + i]*u[ndm + i] - uold[ndm + i]*u[ndm2 + i];
    }

    return 0;
}

int32
stpnhd(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *u) {
    int64 ndim;
    double thta;

    double uold;
    double *smat;

    int64 nfpr1;
    double *f;
    double *v;
    int64 found;
    double c1;

    double s1;

    int64 ndm;
    int64 irs;
    int64 ndm2;

    f = xmalloc(sizeof(*f)*(usize)(iap->ndim));
    v = xmalloc(sizeof(*v)*(usize)(iap->ndim));
    smat = xmalloc(sizeof(*smat)*(usize)((iap->ndim*2)*(iap->ndim*2)));
    // Generates starting data for the continuation of Hopf bifurcation
    // points for maps.

    // Local

    // Parameter adjustments

    ndim = iap->ndim;
    irs = iap->irs;
    ndm = iap->ndm;

    findlb(iap, rap, irs, &nfpr1, &found);
    readlb(u, par);

    thta = pi(2.0) / par[10];
    s1 = sin(thta);
    c1 = cos(thta);
    funi(iap, rap, ndm, u, &uold, icp, par, 1, f, global_scratch.dfu, global_scratch.dfp);

    ndm2 = ndm*2;
    for (int64 i = 0; i < ndm2; ++i) {
        for (int64 j = 0; j < ndm2; ++j) {
            smat[i + j*(ndim*2)] = 0.;
        }
    }

    for (int64 i = 0; i < ndm; ++i) {
        smat[i + (ndm + i)*(ndim*2)] = s1;
    }

    for (int64 i = 0; i < ndm; ++i) {
        smat[ndm + i + i*(ndim*2)] = -s1;
    }

    for (int64 i = 0; i < ndm; ++i) {
        for (int64 j = 0; j < ndm; ++j) {
            smat[i + j*(ndim*2)] = global_scratch.dfu[j*ndm + i];
            smat[ndm + i + (ndm + j)*(ndim*2)] = global_scratch.dfu[j*ndm + i];
        }
        smat[i + i*(ndim*2)] -= c1;
        smat[ndm + i + (ndm + i)*(ndim*2)] -= c1;
    }
    {
        int64 tmp = (ndim*2);
        nlvc(ndm2, tmp, 2, smat, v);
    }
    nrmlz(&ndm2, v);

    for (int64 i = 0; i < ndm2; ++i) {
        u[ndm + i] = v[i];
    }

    u[-1 + ndim - 1] = thta;
    u[-1 + ndim] = par[icp[1]];

    free(smat);
    free(f);
    free(v);
    return 0;
}

/* ----------------------------------------------------------------------- */
/*     Subroutines for the Continuation of Hopf Bifurcation Points (ODE) */
/* ----------------------------------------------------------------------- */

int32
fnhb(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double rtmp;
    double ep;
    int64 ndm;
    double umx;

    // Generates the equations for the 2-parameter continuation of Hopf
    // bifurcation points in ODE.

    // Local

    // Parameter adjustments

    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;

    // Generate the function.

    ffhb(iap, rap, ndim, u, uold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        ffhb(iap, rap, ndim, global_scratch.uu1, uold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        ffhb(iap, rap, ndim, global_scratch.uu2, uold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    par[icp[0]] += ep;

    ffhb(iap, rap, ndim, u, uold, icp, par, global_scratch.ff1, ndm, global_scratch.dfu,
         global_scratch.dfp);

    for (int64 j = 0; j < ndim; ++j) {
        ARRAY2D(dfdp, j, icp[0]) = (global_scratch.ff1[j] - f[j]) / ep;
    }

    par[icp[0]] -= ep;
    return 0;
}

int32
ffhb(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double rom;
    int64 ndm2;

    // Parameter adjustments

    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    ndm2 = ndm*2;

    rom = u[ndim - 2];
    par[10] = rom*pi(2.0);
    par[icp[1]] = u[-1 + ndim];
    funi(iap, rap, ndm, u, uold, icp, par, 1, f, dfdu, dfdp);

    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = u[ndm2 + i];
        f[ndm2 + i] = -u[ndm + i];
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm + i] += rom*ARRAY2D(dfdu, i, j)*u[ndm + j];
            f[ndm2 + i] += rom*ARRAY2D(dfdu, i, j)*u[ndm2 + j];
        }
    }

    f[ndim - 2] = -1.;

    for (int64 i = 0; i < ndm; ++i) {
        f[ndim - 2] = f[ndim - 2] + u[ndm + i]*u[ndm + i] + u[ndm2 + i]*u[ndm2 + i];
    }

    f[-1 + ndim] = 0.;

    for (int64 i = 0; i < ndm; ++i) {
        f[-1 + ndim] = f[-1 + ndim] + uold[ndm2 + i]*(u[ndm + i] - uold[ndm + i]) -
                       uold[ndm + i]*(u[ndm2 + i] - uold[ndm2 + i]);
    }

    return 0;
}

int32
stpnhb(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *u) {
    int64 ndim;

    double uold;
    double *smat;
    int64 nfpr1;
    double *f;
    double *v;
    int64 found;

    double period;
    int64 ndm;
    int64 irs;
    double rom;
    int64 ndm2;
    smat = xmalloc(sizeof(*smat)*(usize)((iap->ndim*2)*(iap->ndim*2)));
    f = xmalloc(sizeof(*f)*(usize)(iap->ndim));
    v = xmalloc(sizeof(*v)*(usize)(iap->ndim));
    // Generates starting data for the 2-parameter continuation of
    // Hopf bifurcation point (ODE).

    // Local

    // Parameter adjustments

    ndim = iap->ndim;
    irs = iap->irs;
    ndm = iap->ndm;

    findlb(iap, rap, irs, &nfpr1, &found);
    readlb(u, par);

    period = par[10];
    rom = period / pi(2.0);
    funi(iap, rap, ndm, u, &uold, icp, par, 1, f, global_scratch.dfu, global_scratch.dfp);

    ndm2 = ndm*2;
    for (int64 i = 0; i < ndm2; ++i) {
        for (int64 j = 0; j < ndm2; ++j) {
            smat[i + j*(ndim*2)] = 0.;
        }
    }

    for (int64 i = 0; i < ndm; ++i) {
        smat[i + (ndm + i)*(ndim*2)] = 1.;
    }

    for (int64 i = 0; i < ndm; ++i) {
        smat[ndm + i + i*(ndim*2)] = -1.;
    }

    for (int64 i = 0; i < ndm; ++i) {
        for (int64 j = 0; j < ndm; ++j) {
            smat[i + j*(ndim*2)] = rom*global_scratch.dfu[j*ndm + i];
            smat[ndm + i + (ndm + j)*(ndim*2)] = rom*global_scratch.dfu[j*ndm + i];
        }
    }
    {
        int64 tmp = (ndim*2);
        nlvc(ndm2, tmp, 2, smat, v);
    }
    nrmlz(&ndm2, v);

    for (int64 i = 0; i < ndm2; ++i) {
        u[ndm + i] = v[i];
    }

    u[ndim - 2] = rom;
    u[-1 + ndim] = par[icp[1]];
    free(smat);
    free(f);
    free(v);
    return 0;
}

/* ----------------------------------------------------------------------- */
/*   Subroutines for the Continuation of Hopf Bifurcation Points (Waves) */
/* ----------------------------------------------------------------------- */

int32
fnhw(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double rtmp;
    double ep;
    int64 ndm;
    double umx;

    // Generates the equations for the 2-parameter continuation of a
    // bifurcation to a traveling wave.

    // Local

    // Parameter adjustments

    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;

    // Generate the function.

    ffhw(iap, rap, ndim, u, uold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        ffhw(iap, rap, ndim, global_scratch.uu1, uold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        ffhw(iap, rap, ndim, global_scratch.uu2, uold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    par[icp[0]] += ep;

    ffhw(iap, rap, ndim, u, uold, icp, par, global_scratch.ff1, ndm, global_scratch.dfu,
         global_scratch.dfp);

    for (int64 j = 0; j < ndim; ++j) {
        ARRAY2D(dfdp, j, icp[0]) = (global_scratch.ff1[j] - f[j]) / ep;
    }

    par[icp[0]] -= ep;

    return 0;
}

int32
ffhw(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 ijac;

    double rom;
    int64 ndm2;

    // Parameter adjustments
    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    ndm2 = ndm*2;

    rom = u[-1 + ndim - 1];
    par[icp[1]] = u[-1 + ndim];
    ijac = 1;
    fnws(iap, rap, ndm, u, uold, icp, par, ijac, f, dfdu, dfdp);

    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = u[ndm2 + i];
        f[ndm2 + i] = -u[ndm + i];
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm + i] += rom*ARRAY2D(dfdu, i, j)*u[ndm + j];
            f[ndm2 + i] += rom*ARRAY2D(dfdu, i, j)*u[ndm2 + j];
        }
    }

    f[ndim - 2] = -1.;

    for (int64 i = 0; i < ndm; ++i) {
        f[ndim - 2] = f[ndim - 2] + u[ndm + i]*u[ndm + i] + u[ndm2 + i]*u[ndm2 + i];
    }

    f[-1 + ndim] = 0.;

    for (int64 i = 0; i < ndm; ++i) {
        f[-1 + ndim] = f[-1 + ndim] + uold[ndm2 + i]*(u[ndm + i] - uold[ndm + i]) -
                       uold[ndm + i]*(u[ndm2 + i] - uold[ndm2 + i]);
    }

    return 0;
}

int32
stpnhw(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *u) {
    int64 ijac;
    int64 ndim;

    double uold;
    double *smat;

    int64 nfpr1;
    double *f;
    double *v;
    int64 found;

    double period, *dfp, *dfu;
    int64 ndm;
    int64 irs;
    double rom;
    int64 ndm2;

    smat = xmalloc(sizeof(*smat)*(usize)((2*iap->ndim)*(2*iap->ndim)));
    f = xmalloc(sizeof(*f)*(usize)(iap->ndim));
    v = xmalloc(sizeof(*v)*(usize)(iap->ndim));
    dfp = xmalloc(sizeof(*dfp)*(usize)((iap->ndim)*NPARX));
    dfu = xmalloc(sizeof(*dfu)*(usize)((iap->ndim)*(iap->ndim)));

    // Generates starting data for the continuation of a bifurcation to a
    // traveling wave.

    // Local (Can't use BLLOC here.)

    // Parameter adjustments
    ndim = iap->ndim;
    irs = iap->irs;
    ndm = iap->ndm;

    findlb(iap, rap, irs, &nfpr1, &found);
    readlb(u, par);

    ijac = 1;
    period = par[10];
    rom = period / pi(2.0);
    fnws(iap, rap, ndm, u, &uold, icp, par, ijac, f, dfu, dfp);

    ndm2 = ndm*2;
    for (int64 i = 0; i < ndm2; ++i) {
        for (int64 j = 0; j < ndm2; ++j) {
            smat[i + j*(ndim*2)] = 0.;
        }
    }

    for (int64 i = 0; i < ndm; ++i) {
        smat[i + (ndm + i)*(ndim*2)] = 1.;
    }

    for (int64 i = 0; i < ndm; ++i) {
        smat[ndm + i + i*(ndim*2)] = -1.;
    }

    for (int64 i = 0; i < ndm; ++i) {
        for (int64 j = 0; j < ndm; ++j) {
            smat[i + j*(ndim*2)] = rom*global_scratch.dfu[j*ndm + i];
            smat[ndm + i + (ndm + j)*(ndim*2)] = rom*global_scratch.dfu[j*ndm + i];
        }
    }
    {
        int64 tmp = (ndim*2);
        nlvc(ndm2, tmp, 2, smat, v);
    }
    nrmlz(&ndm2, v);

    for (int64 i = 0; i < ndm2; ++i) {
        u[ndm + i] = v[i];
    }

    u[ndim - 2] = rom;
    u[-1 + ndim] = par[icp[1]];

    free(smat);
    free(f);
    free(v);
    free(dfp);
    free(dfu);

    return 0;
}

/* ----------------------------------------------------------------------- */
/*          Periodic Solutions and Fixed Period Orbits */
/* ----------------------------------------------------------------------- */

int32
fnps(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double period;

    // Generates the equations for the continuation of periodic orbits.

    // Generate the function.

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    if (icp[1] == 10) {
        //          **Variable period continuation
        funi(iap, rap, ndim, u, uold, icp, par, ijac, f, dfdu, dfdp);
        period = par[10];
        for (int64 i = 0; i < ndim; ++i) {
            ARRAY2D(dfdp, i, 10) = f[i];
            f[i] = period*ARRAY2D(dfdp, i, 10);
        }
        if (ijac == 0) {
            return 0;
        }
        //          **Generate the Jacobian.
        for (int64 i = 0; i < ndim; ++i) {
            for (int64 j = 0; j < ndim; ++j) {
                ARRAY2D(dfdu, i, j) = period*ARRAY2D(dfdu, i, j);
            }
            ARRAY2D(dfdp, i, (icp[0])) = period*ARRAY2D(dfdp, i, (icp[0]));
        }
    } else {
        //          **Fixed period continuation
        period = par[10];
        funi(iap, rap, ndim, u, uold, icp, par, ijac, f, dfdu, dfdp);
        for (int64 i = 0; i < ndim; ++i) {
            f[i] = period*f[i];
        }
        if (ijac == 0) {
            return 0;
        }
        //          **Generate the Jacobian.
        for (int64 i = 0; i < ndim; ++i) {
            for (int64 j = 0; j < ndim; ++j) {
                ARRAY2D(dfdu, i, j) = period*ARRAY2D(dfdu, i, j);
            }
            for (int64 j = 0; j < 2; ++j) {
                ARRAY2D(dfdp, i, icp[j]) = period*ARRAY2D(dfdp, i, icp[j]);
            }
        }
    }

    return 0;
}

int32
bcps(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nbc, double *u0,
     double *u1, double *f, int64 ijac, double *dbc) {
    int64 dbc_dim1;

    int64 jtmp;
    int64 nn;

    (void)iap;
    (void)rap;
    (void)icp;

    // Parameter adjustments

    dbc_dim1 = nbc;

    for (int64 i = 0; i < ndim; ++i) {
        f[i] = u0[i] - u1[i];
    }

    // Rotations
    if (global_rotations.irtn != 0) {
        for (int64 i = 0; i < ndim; ++i) {
            if (global_rotations.nrtn[i] != 0) {
                f[i] += par[18]*(double)global_rotations.nrtn[i];
            }
        }
    }

    if (ijac == 0) {
        return 0;
    }

    jtmp = NPARX;
    nn = (ndim*2) + jtmp;
    for (int64 i = 0; i < nbc; ++i) {
        for (int64 j = 0; j < nn; ++j) {
            ARRAY2D(dbc, i, j) = 0.;
        }
    }

    for (int64 i = 0; i < ndim; ++i) {
        ARRAY2D(dbc, i, i) = 1.;
        ARRAY2D(dbc, i, (ndim + i)) = -1.;
    }

    return 0;
}

int32
icps(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, double *u,
     double *uold, double *udot, double *upold, double *f, int64 ijac, double *dint) {
    int64 dint_dim1;

    int64 jtmp;
    int64 nn;

    (void)iap;
    (void)rap;
    (void)par;
    (void)icp;
    (void)uold;
    (void)udot;

    // Parameter adjustments

    dint_dim1 = nint;

    f[0] = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (global_rotations.nrtn[i] == 0) {
            f[0] += u[i]*upold[i];
        }
    }

    if (ijac == 0) {
        return 0;
    }

    jtmp = NPARX;
    nn = ndim + jtmp;
    for (int64 i = 0; i < nn; ++i) {
        ARRAY2D(dint, 0, i) = 0.;
    }

    for (int64 i = 0; i < ndim; ++i) {
        if (global_rotations.nrtn[i] == 0) {
            ARRAY2D(dint, 0, i) = upold[i];
        } else {
            ARRAY2D(dint, 0, i) = 0.;
        }
    }

    return 0;
}

int32
pdble(iap_type *iap, rap_type *rap, int64 *ndim, int64 *ntst, int64 *ncol, int64 *ndxloc,
      double *ups, double *udotps, double *tm, double *par) {
    int64 ups_dim1;
    int64 udotps_dim1;

    int64 i;
    int64 i1;
    int64 i2;

    (void)iap;
    (void)rap;

    // Preprocesses restart data for switching branches at a period doubling

    // Parameter adjustments
    udotps_dim1 = *ndxloc;
    ups_dim1 = *ndxloc;

    par[10] *= 2.;
    if (global_rotations.irtn != 0) {
        par[18] *= 2.;
    }

    for (i = 0; i < *ntst; ++i) {
        tm[i] *= .5;
        tm[*ntst + i] = tm[i] + .5;
    }

    tm[(*ntst*2)] = 1.;

    for (int64 j = 0; j < *ntst + 1; ++j) {
        for (i1 = 0; i1 < *ndim; ++i1) {
            for (i2 = 0; i2 < *ncol; ++i2) {
                i = i2**ndim + i1;
                ARRAY2D(ups, *ntst + j, i) =
                    ARRAY2D(ups, *ntst, i1) + ARRAY2D(ups, j, i) - ARRAY2D(ups, 0, i1);
                ARRAY2D(udotps, *ntst + j, i) =
                    ARRAY2D(udotps, *ntst, i1) + ARRAY2D(udotps, j, i) - ARRAY2D(udotps, 0, i);
            }
        }
    }

    *ntst *= 2;

    return 0;
}

int32
stpnps(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr, int64 *ncolrs,
       double *rlcur, double *rldot, int64 *ndxloc, double *ups, double *udotps, double *upoldp,
       double *tm, double *dtm, int64 *nodir, double *thl, double *thu) {
    int64 ups_dim1;
    int64 udotps_dim1;
    int64 upoldp_dim1;

    int64 ndim;
    int64 ncol;

    double uold;
    double *smat;
    int64 nfpr;
    int64 ntst;
    int64 ndim2;
    int64 nfpr1;
    double c;
    double *f;
    double s, t, *u, rimhb;
    int64 found;
    int64 k1;
    double *rnllv;

    double dt;

    double period;

    double tpi;
    int64 irs;

    smat = xmalloc(sizeof(*smat)*(usize)((iap->ndim*2)*(iap->ndim*2)));
    rnllv = xmalloc(sizeof(*rnllv)*(usize)((iap->ndim*2)*(iap->ndim*2)));
    f = xmalloc(sizeof(*f)*(usize)(iap->ndim));
    u = xmalloc(sizeof(*u)*(usize)(iap->ndim));
    // Generates starting data for the continuation of a branch of periodic
    // solutions from a Hopf bifurcation point.

    // Local

    // Parameter adjustments
    upoldp_dim1 = *ndxloc;
    udotps_dim1 = *ndxloc;
    ups_dim1 = *ndxloc;

    ndim = iap->ndim;
    irs = iap->irs;
    ntst = iap->ntst;
    ncol = iap->ncol;
    nfpr = iap->nfpr;

    findlb(iap, rap, irs, &nfpr1, &found);
    readlb(u, par);

    for (int64 i = 0; i < nfpr; ++i) {
        rlcur[i] = par[icp[i]];
    }

    period = par[10];
    tpi = pi(2.0);
    rimhb = tpi / period;
    *ntsr = ntst;
    *ncolrs = ncol;

    ndim2 = ndim*2;
    for (int64 i = 0; i < ndim2; ++i) {
        for (int64 j = 0; j < ndim2; ++j) {
            smat[i + j*(ndim*2)] = 0.;
        }
    }

    for (int64 i = 0; i < ndim; ++i) {
        smat[i + i*(ndim*2)] = -rimhb;
        smat[ndim + i + (ndim + i)*(ndim*2)] = rimhb;
    }

    funi(iap, rap, ndim, u, &uold, icp, par, 1, f, global_scratch.dfu, global_scratch.dfp);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            smat[i + (ndim + j)*(ndim*2)] = global_scratch.dfu[j*ndim + i];
            smat[ndim + i + j*(ndim*2)] = global_scratch.dfu[j*ndim + i];
        }
    }

    {
        int64 tmp = (ndim*2);
        nlvc(ndim2, tmp, 2, smat, rnllv);
    }
    nrmlz(&ndim2, rnllv);

    // Generate the (initially uniform) mesh.

    msh(iap, tm);
    dt = 1. / (double)ntst;

    for (int64 j = 0; j < ntst + 1; ++j) {
        t = tm[j];
        s = sin(tpi*t);
        c = cos(tpi*t);
        for (int64 k = 0; k < ndim; ++k) {
            ARRAY2D(udotps, j, k) = s*rnllv[k] + c*rnllv[ndim + k];
            ARRAY2D(upoldp, j, k) = c*rnllv[k] - s*rnllv[ndim + k];
            ARRAY2D(ups, j, k) = u[k];
        }
    }

    for (int64 i = 0; i < ncol - 1; ++i) {
        for (int64 j = 0; j < ntst; ++j) {
            t = tm[j] + (double)(i + 1)*(tm[j + 1] - tm[j]) / (double)ncol;
            s = sin(tpi*t);
            c = cos(tpi*t);
            for (int64 k = 0; k < ndim; ++k) {
                k1 = (i + 1)*ndim + k;
                ARRAY2D(udotps, j, k1) = s*rnllv[k] + c*rnllv[ndim + k];
                ARRAY2D(upoldp, j, k1) = c*rnllv[k] - s*rnllv[ndim + k];
                ARRAY2D(ups, j, k1) = u[k];
            }
        }
    }

    rldot[0] = 0.;
    rldot[1] = 0.;

    for (int64 i = 0; i < ntst; ++i) {
        dtm[i] = dt;
    }

    scaleb(iap, icp, ndxloc, udotps, rldot, dtm, thl, thu);

    *nodir = -1;
    free(smat);
    free(rnllv);
    free(f);
    free(u);

    return 0;
}

/* ----------------------------------------------------------------------- */
/*          Travelling Wave Solutions to Parabolic PDEs */
/* ----------------------------------------------------------------------- */

int32
fnws(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 ndm;
    int64 ndm2;

    // Sets up equations for the continuation of spatially homogeneous
    // solutions to parabolic systems, for the purpose of finding
    // bifurcations to travelling wave solutions.

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;

    // Generate the function.

    ndm2 = ndm / 2;
    ffws(iap, rap, ndim, u, uold, icp, par, ijac, f, dfdu, dfdp, ndm2, global_scratch.dfu,
         global_scratch.dfp);

    return 0;
}

int32
ffws(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp, int64 ndm, double *dfu, double *dfp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;
    int64 dfu_dim1;
    int64 dfp_dim1;

    int64 nfpr;
    double c;

    // Parameter adjustments

    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;
    dfp_dim1 = ndm;
    dfu_dim1 = ndm;

    nfpr = iap->nfpr;

    c = par[9];
    funi(iap, rap, ndm, u, uold, icp, par, ijac, f, dfu, dfp);

    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = -(c*u[ndm + i] + f[i]) / par[i + 14];
        f[i] = u[ndm + i];
    }

    if (ijac == 0) {
        return 0;
    }

    for (int64 i = 0; i < ndm; ++i) {
        for (int64 j = 0; j < ndm; ++j) {
            ARRAY2D(dfdu, i, j) = 0.;
            ARRAY2D(dfdu, i, (j + ndm)) = 0.;
            ARRAY2D(dfdu, i + ndm, j) = -ARRAY2D(dfu, i, j) / par[i + 14];
            ARRAY2D(dfdu, i + ndm, (j + ndm)) = 0.;
        }
        ARRAY2D(dfdu, i, (i + ndm)) = 1.;
        ARRAY2D(dfdu, i + ndm, (i + ndm)) = -c / par[i + 14];
        if (icp[0] < 9) {
            ARRAY2D(dfdp, i, (icp[0])) = 0.;
            ARRAY2D(dfdp, i + ndm, icp[0]) = -ARRAY2D(dfp, i, icp[0]) / par[i + 14];
        }
        if (nfpr > 1 && icp[1] < 9) {
            ARRAY2D(dfdp, i, (icp[1])) = 0.;
            ARRAY2D(dfdp, i + ndm, icp[1]) = -ARRAY2D(dfp, i, icp[1]) / par[i + 14];
        }
    }

    // Derivative with respect to the wave speed.

    for (int64 i = 0; i < ndm; ++i) {
        ARRAY2D(dfdp, i, 9) = 0.;
        ARRAY2D(dfdp, i + ndm, 9) = -u[ndm + i] / par[i + 14];
    }

    // Derivatives with respect to the diffusion coefficients.

    for (int64 j = 0; j < ndm; ++j) {
        for (int64 i = 0; i < ndm; ++i) {
            ARRAY2D(dfdp, i, (j + 14)) = 0.;
            ARRAY2D(dfdp, i + ndm, (j + 14)) = 0.;
        }
        ARRAY2D(dfdp, j + ndm, (j + 14)) = -f[j + ndm] / par[j + 14];
    }

    return 0;
}

int32
fnwp(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double period;

    // Equations for the continuation of traveling waves.

    // Generate the function and Jacobian.

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    if (icp[1] == 10) {
        //          **Variable wave length
        fnws(iap, rap, ndim, u, uold, icp, par, ijac, f, dfdu, dfdp);
        period = par[10];
        for (int64 i = 0; i < ndim; ++i) {
            ARRAY2D(dfdp, i, 10) = f[i];
            f[i] = period*f[i];
        }
        if (ijac == 0) {
            return 0;
        }
        for (int64 i = 0; i < ndim; ++i) {
            for (int64 j = 0; j < ndim; ++j) {
                ARRAY2D(dfdu, i, j) = period*ARRAY2D(dfdu, i, j);
            }
        }
        for (int64 i = 0; i < ndim; ++i) {
            ARRAY2D(dfdp, i, (icp[0])) = period*ARRAY2D(dfdp, i, (icp[0]));
        }
    } else {
        //          **Fixed wave length
        fnws(iap, rap, ndim, u, uold, icp, par, ijac, f, dfdu, dfdp);
        period = par[10];
        for (int64 i = 0; i < ndim; ++i) {
            f[i] = period*f[i];
        }
        if (ijac == 0) {
            return 0;
        }
        for (int64 i = 0; i < ndim; ++i) {
            for (int64 j = 0; j < ndim; ++j) {
                ARRAY2D(dfdu, i, j) = period*ARRAY2D(dfdu, i, j);
            }
        }
        for (int64 i = 0; i < ndim; ++i) {
            for (int64 j = 0; j < 2; ++j) {
                ARRAY2D(dfdp, i, icp[j]) = period*ARRAY2D(dfdp, i, icp[j]);
            }
        }
    }

    return 0;
}

int32
stpnwp(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr, int64 *ncolrs,
       double *rlcur, double *rldot, int64 *ndxloc, double *ups, double *udotps, double *upoldp,
       double *tm, double *dtm, int64 *nodir, double *thl, double *thu) {
    int64 ups_dim1;
    int64 udotps_dim1;
    int64 upoldp_dim1;

    int64 ndim;
    int64 ncol;

    double uold;
    double *smat;
    int64 nfpr;

    int64 ntst;
    int64 ndim2;
    int64 nfpr1;
    double c;
    double *f;
    double s, t, *u, rimhb;
    int64 found;
    int64 k1;
    double *rnllv;

    double dt;

    double period, *dfp, *dfu;

    double tpi;
    int64 irs;

    smat = xmalloc(sizeof(*smat)*(usize)((2*iap->ndim)*(2*iap->ndim)));
    f = xmalloc(sizeof(*f)*(usize)(iap->ndim));
    u = xmalloc(sizeof(*u)*(usize)(iap->ndim));
    rnllv = xmalloc(sizeof(*rnllv)*(usize)(2*(iap->ndim)));
    dfp = xmalloc(sizeof(*dfp)*(usize)((iap->ndim)*NPARX));
    dfu = xmalloc(sizeof(*dfu)*(usize)((iap->ndim)*(iap->ndim)));

    // Generates starting data for the continuation of a branch of periodic
    // solutions starting from a Hopf bifurcation point (Waves).

    // Local (Can't use BLLOC here.)

    // Parameter adjustments
    upoldp_dim1 = *ndxloc;
    udotps_dim1 = *ndxloc;
    ups_dim1 = *ndxloc;

    ndim = iap->ndim;
    irs = iap->irs;
    ntst = iap->ntst;
    ncol = iap->ncol;
    nfpr = iap->nfpr;

    findlb(iap, rap, irs, &nfpr1, &found);
    readlb(u, par);

    for (int64 i = 0; i < nfpr; ++i) {
        rlcur[i] = par[icp[i]];
    }

    period = par[10];
    tpi = pi(2.0);
    rimhb = tpi / period;
    *ntsr = ntst;
    *ncolrs = ncol;

    ndim2 = ndim*2;
    for (int64 i = 0; i < ndim2; ++i) {
        for (int64 j = 0; j < ndim2; ++j) {
            smat[i + j*(ndim*2)] = 0.;
        }
    }

    for (int64 i = 0; i < ndim; ++i) {
        smat[i + i*(ndim*2)] = -rimhb;
        smat[ndim + i + (ndim + i)*(ndim*2)] = rimhb;
    }

    fnws(iap, rap, ndim, u, &uold, icp, par, 1, f, dfu, dfp);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            smat[i + (ndim + j)*(ndim*2)] = dfu[j*ndim + i];
            smat[ndim + i + j*(ndim*2)] = dfu[j*ndim + i];
        }
    }

    nlvc(ndim2, ndim*2, 2, smat, rnllv);
    nrmlz(&ndim2, rnllv);

    // Generate the (initially uniform) mesh.

    msh(iap, tm);
    dt = 1. / (double)ntst;

    for (int64 j = 0; j < ntst + 1; ++j) {
        t = tm[j];
        s = sin(tpi*t);
        c = cos(tpi*t);
        for (int64 k = 0; k < ndim; ++k) {
            ARRAY2D(udotps, j, k) = s*rnllv[k] + c*rnllv[ndim + k];
            ARRAY2D(upoldp, j, k) = c*rnllv[k] - s*rnllv[ndim + k];
            ARRAY2D(ups, j, k) = u[k];
        }
    }

    for (int64 i = 0; i < ncol - 1; ++i) {
        for (int64 j = 0; j < ntst; ++j) {
            t = tm[j] + (double)(i + 1)*(tm[j + 1] - tm[j]) / (double)ncol;
            s = sin(tpi*t);
            c = cos(tpi*t);
            for (int64 k = 0; k < ndim; ++k) {
                k1 = (i + 1)*ndim + k;
                ARRAY2D(udotps, j, k1) = s*rnllv[k] + c*rnllv[ndim + k];
                ARRAY2D(upoldp, j, k1) = c*rnllv[k] - s*rnllv[ndim + k];
                ARRAY2D(ups, j, k1) = u[k];
            }
        }
    }

    rldot[0] = 0.;
    rldot[1] = 0.;

    for (int64 i = 0; i < ntst; ++i) {
        dtm[i] = dt;
    }

    scaleb(iap, icp, ndxloc, udotps, rldot, dtm, thl, thu);

    *nodir = -1;

    free(smat);
    free(f);
    free(u);
    free(rnllv);
    free(dfp);
    free(dfu);

    return 0;
}

/* ----------------------------------------------------------------------- */
/*             Parabolic PDEs : Stationary States */
/* ----------------------------------------------------------------------- */

int32
fnsp(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 ndm;

    // Generates the equations for taking one time step (Implicit Euler).

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;

    // Generate the function and Jacobian.

    ffsp(iap, rap, ndim, u, uold, icp, par, ijac, f, dfdu, dfdp, ndm, global_scratch.dfu,
         global_scratch.dfp);

    return 0;
}

int32
ffsp(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp, int64 ndm, double *dfu, double *dfp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;
    int64 dfu_dim1;
    int64 dfp_dim1;

    double period;

    // Parameter adjustments

    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;
    dfp_dim1 = ndm;
    dfu_dim1 = ndm;

    funi(iap, rap, ndm, u, uold, icp, par, ijac, &f[ndm], dfu, dfp);

    period = par[10];
    for (int64 i = 0; i < ndm; ++i) {
        f[i] = period*u[ndm + i];
        f[ndm + i] = -period*f[ndm + i] / par[i + 14];
    }

    if (ijac == 0) {
        return 0;
    }

    for (int64 i = 0; i < ndm; ++i) {
        for (int64 j = 0; j < ndm; ++j) {
            ARRAY2D(dfdu, i, j) = 0.;
            ARRAY2D(dfdu, i, (j + ndm)) = 0.;
            ARRAY2D(dfdu, i + ndm, j) = -period*ARRAY2D(dfu, i, j) / par[i + 14];
            ARRAY2D(dfdu, i + ndm, (j + ndm)) = 0.;
        }
        ARRAY2D(dfdu, i, (i + ndm)) = period;
        if (icp[0] == 10) {
            ARRAY2D(dfdp, i, (icp[0])) = f[i] / period;
            ARRAY2D(dfdp, ndm + i, icp[0]) = f[ndm + i] / period;
        } else if (icp[0] == i + 13) {
            ARRAY2D(dfdp, i, (icp[0])) = 0.;
            ARRAY2D(dfdp, ndm + i, icp[0]) = -f[ndm + i] / par[i + 14];
        } else if (icp[0] != 10 && !(icp[0] > 13 && icp[0] <= ndm + 13)) {
            ARRAY2D(dfdp, i, (icp[0])) = 0.;
            ARRAY2D(dfdp, i + ndm, icp[0]) = -period*ARRAY2D(dfp, i, icp[0]) / par[i + 14];
        }
    }

    return 0;
}

/* ----------------------------------------------------------------------- */
/*            Time Evolution of Parabolic PDEs */
/* ----------------------------------------------------------------------- */

int32
fnpe(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 ndm;

    // Generates the equations for taking one time step (Implicit Euler).

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;

    // Generate the function and Jacobian.
    ffpe(iap, rap, ndim, u, uold, icp, par, ijac, f, dfdu, dfdp, ndm, global_scratch.dfu,
         global_scratch.dfp);

    return 0;
}

int32
ffpe(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp, int64 ndm, double *dfu, double *dfp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;
    int64 dfu_dim1;
    int64 dfp_dim1;

    double t;
    double dsmin;
    double rlold;
    double ds;
    double dt;
    double period;

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;
    dfp_dim1 = ndm;
    dfu_dim1 = ndm;

    ds = rap->ds;
    dsmin = rap->dsmin;

    period = par[10];
    t = par[icp[0]];
    rlold = rap->tivp;
    dt = t - rlold;
    if (fabs(dt) < dsmin) {
        dt = ds;
    }

    funi(iap, rap, ndm, u, uold, icp, par, ijac, &f[ndm], dfu, dfp);

    for (int64 i = 0; i < ndm; ++i) {
        f[i] = period*u[ndm + i];
        f[ndm + i] = period*((u[i] - uold[i]) / dt - f[ndm + i]) / par[i + 14];
    }

    if (ijac == 0) {
        return 0;
    }

    for (int64 i = 0; i < ndm; ++i) {
        for (int64 j = 0; j < ndm; ++j) {
            ARRAY2D(dfdu, i, j) = 0.;
            ARRAY2D(dfdu, i, (j + ndm)) = 0.;
            ARRAY2D(dfdu, i + ndm, j) = -period*ARRAY2D(dfu, i, j) / par[i + 14];
            ARRAY2D(dfdu, i + ndm, (j + ndm)) = 0.;
        }
        ARRAY2D(dfdu, i, (i + ndm)) = period;
        ARRAY2D(dfdu, i + ndm, i) += period / (dt*par[i + 14]);
        ARRAY2D(dfdp, i, (icp[0])) = 0.;
        // Computing 2nd power
        ARRAY2D(dfdp, i + ndm, icp[0]) = -period*(u[i] - uold[i]) / (dt*dt*par[i + 14]);
    }

    return 0;
}

int32
icpe(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, double *u,
     double *uold, double *udot, double *upold, double *f, int64 ijac, double *dint) {
    (void)iap;
    (void)rap;
    (void)ndim;
    (void)par;
    (void)icp;
    (void)nint;
    (void)u;
    (void)uold;
    (void)udot;
    (void)upold;
    (void)f;
    (void)ijac;
    (void)dint;

    // Dummy integral condition subroutine for parabolic systems.

    return 0;
}

/* ----------------------------------------------------------------------- */
/*    Subroutines for the Continuation of Folds for Periodic Solution */
/* ----------------------------------------------------------------------- */

int32
fnpl(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 nfpr;
    double rtmp;
    double ep;
    int64 ndm;
    double umx;

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;
    nfpr = iap->nfpr;

    // Generate the function.

    ffpl(iap, rap, ndim, u, uold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        ffpl(iap, rap, ndim, global_scratch.uu1, uold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        ffpl(iap, rap, ndim, global_scratch.uu2, uold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        ffpl(iap, rap, ndim, u, uold, icp, par, global_scratch.ff1, ndm, global_scratch.dfu,
             global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdp, j, icp[i]) = (global_scratch.ff1[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }

    return 0;
}

int32
ffpl(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double beta;

    double period;
    int64 ips;

    (void)ndim;

    // Parameter adjustments
    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    period = par[10];
    beta = par[11];
    funi(iap, rap, ndm, u, uold, icp, par, 2, f, dfdu, dfdp);

    ips = iap->ips;
    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = 0.;
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm + i] += ARRAY2D(dfdu, i, j)*u[ndm + j];
        }
        if (icp[2] == 10) {
            //            ** Variable period
            f[ndm + i] = period*f[ndm + i] + beta*f[i];
        } else {
            //            ** Fixed period
            f[ndm + i] = period*f[ndm + i] + beta*ARRAY2D(dfdp, i, (icp[1]));
        }
        f[i] = period*f[i];
    }

    return 0;
}

int32
bcpl(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nbc, double *u0,
     double *u1, double *f, int64 ijac, double *dbc) {
    int64 dbc_dim1;
    int64 jtmp;
    int64 nn;
    int64 ndm;

    (void)rap;
    (void)icp;

    // Boundary conditions for continuing folds (Periodic solutions)

    // Parameter adjustments
    dbc_dim1 = nbc;

    for (int64 i = 0; i < ndim; ++i) {
        f[i] = u0[i] - u1[i];
    }

    // Rotations
    if (global_rotations.irtn != 0) {
        ndm = iap->ndm;
        for (int64 i = 0; i < ndm; ++i) {
            if (global_rotations.nrtn[i] != 0) {
                f[i] += par[18]*(double)global_rotations.nrtn[i];
            }
        }
    }

    if (ijac == 0) {
        return 0;
    }

    jtmp = NPARX;
    nn = (ndim*2) + jtmp;
    for (int64 i = 0; i < nbc; ++i) {
        for (int64 j = 0; j < nn; ++j) {
            ARRAY2D(dbc, i, j) = 0.;
        }
    }

    for (int64 i = 0; i < ndim; ++i) {
        ARRAY2D(dbc, i, i) = 1.;
        ARRAY2D(dbc, i, (ndim + i)) = -1.;
    }

    return 0;
}

int32
icpl(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, double *u,
     double *uold, double *udot, double *upold, double *f, int64 ijac, double *dint) {
    int64 dint_dim1;

    int64 jtmp;
    int64 nn;
    int64 ndm;

    (void)udot;
    (void)rap;
    (void)icp;
    (void)uold;

    // Integral conditions for continuing folds (Periodic solutions)

    // Parameter adjustments
    dint_dim1 = nint;

    ndm = iap->ndm;

    f[0] = 0.;
    f[1] = 0.;
    // Computing 2nd power
    f[2] = par[11]*par[11] - par[12];

    for (int64 i = 0; i < ndm; ++i) {
        if (global_rotations.nrtn[i] == 0) {
            f[0] += u[i]*upold[i];
            f[1] += u[ndm + i]*upold[i];
        }
        f[2] += u[ndm + i]*u[ndm + i];
    }

    if (ijac == 0) {
        return 0;
    }

    jtmp = NPARX;
    nn = ndim + jtmp;
    for (int64 i = 0; i < nint; ++i) {
        for (int64 j = 0; j < nn; ++j) {
            ARRAY2D(dint, i, j) = 0.;
        }
    }

    for (int64 i = 0; i < ndm; ++i) {
        if (global_rotations.nrtn[i] == 0) {
            ARRAY2D(dint, 0, i) = upold[i];
            ARRAY2D(dint, 1, ndm + i) = upold[i];
        } else {
            ARRAY2D(dint, 0, i) = 0.;
            ARRAY2D(dint, 1, ndm + i) = 0.;
        }
        ARRAY2D(dint, 2, ndm + i) = u[ndm + i]*2.;
    }

    ARRAY2D(dint, 2, ndim + 11) = par[11]*2.;
    ARRAY2D(dint, 2, ndim + 12) = -1.;

    return 0;
}

int32
stpnpl(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr, int64 *ncolrs,
       double *rlcur, double *rldot, int64 *ndxloc, double *ups, double *udotps, double *upoldp,
       double *tm, double *dtm, int64 *nodir, double *thl, double *thu) {
    int64 ups_dim1;
    int64 udotps_dim1;
    int64 upoldp_dim1;

    int64 ndim;
    double temp[7];
    int64 nfpr;
    int64 nfpr1;
    int64 ntpl1;
    int64 nrsp1;
    int64 ntot1;
    int64 found;
    int64 icprs[NPARX];
    int64 nparr;
    int64 k1;
    int64 k2;
    int64 nskip1;

    double rd1;
    double rd2;
    int64 ibr;
    int64 ndm;
    int64 ips;
    int64 irs;
    int64 lab1;
    int64 nar1;
    int64 itp1;
    int64 isw1;

    (void)thu;
    (void)thl;
    (void)dtm;
    (void)upoldp;

    // Generates starting data for the 2-parameter continuation of folds
    // on a branch of periodic solutions.

    // Local

    // Parameter adjustments
    upoldp_dim1 = *ndxloc;
    udotps_dim1 = *ndxloc;
    ups_dim1 = *ndxloc;

    ndim = iap->ndim;
    ips = iap->ips;
    irs = iap->irs;
    ndm = iap->ndm;
    nfpr = iap->nfpr;
    ibr = iap->ibr;

    findlb(iap, rap, irs, &nfpr1, &found);
    fscanf(fp3, "%ld", &ibr);
    fscanf(fp3, "%ld", &ntot1);
    fscanf(fp3, "%ld", &itp1);
    fscanf(fp3, "%ld", &lab1);
    fscanf(fp3, "%ld", &nfpr1);
    fscanf(fp3, "%ld", &isw1);
    fscanf(fp3, "%ld", &ntpl1);
    fscanf(fp3, "%ld", &nar1);
    fscanf(fp3, "%ld", &nskip1);
    fscanf(fp3, "%ld", &(*ntsr));
    fscanf(fp3, "%ld", &(*ncolrs));
    fscanf(fp3, "%ld", &nparr);
    iap->ibr = ibr;
    nrsp1 = *ntsr + 1;

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim;
            k2 = k1 + ndm - 1;
            fscanf(fp3, "%lf", &temp[i]);
            for (int64 k = k1; k <= k2; ++k) {
                fscanf(fp3, "%lf", &ARRAY2D(ups, j, k));
            }
        }
        tm[j] = temp[0];
    }

    fscanf(fp3, "%lf", &tm[-1 + nrsp1]);
    for (int64 k = 0; k < ndm; ++k) {
        fscanf(fp3, "%lf", &ARRAY2D(ups, *ntsr, k));
    }

    fscanf(fp3, "%ld", icprs);
    fscanf(fp3, "%ld", &icprs[1]);
    fscanf(fp3, "%lf", &rd1);
    fscanf(fp3, "%lf", &rd2);

    // Read U-dot (derivative with respect to arclength).

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim;
            k2 = k1 + ndm - 1;
            for (int64 k = k1; k <= k2; ++k) {
                fscanf(fp3, "%lf", &ARRAY2D(udotps, j, k));
            }
        }
    }

    for (int64 k = 0; k < ndm; ++k) {
        fscanf(fp3, "%lf", &ARRAY2D(udotps, *ntsr, k));
    }

    // Read the parameter values.

    if (nparr > NPARX) {
        nparr = NPARX;
        printf("Warning : NPARX too small for restart data\n");
        printf("PAR(i) set to zero, fot i > %3ld\n", nparr);
    }
    for (int64 i = 0; i < nparr; ++i) {
        fscanf(fp3, "%lf", &par[i]);
    }

    // Complement starting data
    par[11] = 0.;
    par[12] = 0.;
    if (icp[2] == 10) {
        //          Variable period
        rldot[0] = rd1;
        rldot[1] = 0.;
        rldot[2] = rd2;
        rldot[3] = 0.;
        //          Variable period
    } else {
        //          Fixed period
        rldot[0] = rd1;
        rldot[1] = rd2;
        rldot[2] = 0.;
        rldot[3] = 0.;
    }

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim + ndm;
            k2 = (i + 1)*ndim - 1;
            for (int64 k = k1; k <= k2; ++k) {
                ARRAY2D(ups, j, k) = 0.;
                ARRAY2D(udotps, j, k) = 0.;
            }
        }
    }

    for (int64 k = ndm; k < ndim; ++k) {
        ARRAY2D(ups, nrsp1 - 1, k) = 0.;
        ARRAY2D(udotps, nrsp1 - 1, k) = 0.;
    }

    for (int64 i = 0; i < nfpr; ++i) {
        rlcur[i] = par[icp[i]];
    }

    *nodir = 0;

    return 0;
}

/* ----------------------------------------------------------------------- */
/*   Subroutines for the Continuation of Period Doubling Bifurcations */
/* ----------------------------------------------------------------------- */

int32
fnpd(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 nfpr;
    double rtmp;
    double ep;
    int64 ndm;
    double umx;

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;
    nfpr = iap->nfpr;

    // Generate the function.

    ffpd(iap, rap, ndim, u, uold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        ffpd(iap, rap, ndim, global_scratch.uu1, uold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        ffpd(iap, rap, ndim, global_scratch.uu2, uold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        ffpd(iap, rap, ndim, u, uold, icp, par, global_scratch.ff1, ndm, global_scratch.dfu,
             global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdp, j, icp[i]) = (global_scratch.ff1[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }

    return 0;
}

int32
ffpd(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double period;

    (void)ndim;

    // Parameter adjustments
    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    period = par[10];
    funi(iap, rap, ndm, u, uold, icp, par, 1, f, dfdu, dfdp);

    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = 0.;
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm + i] += ARRAY2D(dfdu, i, j)*u[ndm + j];
        }
        f[i] = period*f[i];
        f[ndm + i] = period*f[ndm + i];
    }

    return 0;
}

int32
bcpd(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nbc, double *u0,
     double *u1, double *f, int64 ijac, double *dbc) {
    int64 dbc_dim1;

    int64 jtmp;
    int64 nn;
    int64 ndm;

    (void)icp;
    (void)rap;

    // Generate boundary conditions for the 2-parameter continuation
    // of period doubling bifurcations.

    // Parameter adjustments
    dbc_dim1 = nbc;

    ndm = iap->ndm;

    for (int64 i = 0; i < ndm; ++i) {
        f[i] = u0[i] - u1[i];
        f[ndm + i] = u0[ndm + i] + u1[ndm + i];
    }

    // Rotations
    if (global_rotations.irtn != 0) {
        for (int64 i = 0; i < ndm; ++i) {
            if (global_rotations.nrtn[i] != 0) {
                f[i] += par[18]*(double)global_rotations.nrtn[i];
            }
        }
    }

    if (ijac == 0) {
        return 0;
    }

    jtmp = NPARX;
    nn = (ndim*2) + jtmp;
    for (int64 i = 0; i < nbc; ++i) {
        for (int64 j = 0; j < nn; ++j) {
            ARRAY2D(dbc, i, j) = 0.;
        }
    }

    for (int64 i = 0; i < ndim; ++i) {
        ARRAY2D(dbc, i, i) = 1.;
        if ((i + 1) <= ndm) {
            ARRAY2D(dbc, i, (ndim + i)) = -1.;
        } else {
            ARRAY2D(dbc, i, (ndim + i)) = 1.;
        }
    }

    return 0;
}

int32
icpd(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, double *u,
     double *uold, double *udot, double *upold, double *f, int64 ijac, double *dint) {
    int64 dint_dim1;

    int64 jtmp;
    int64 nn;
    int64 ndm;

    (void)udot;
    (void)uold;
    (void)icp;
    (void)rap;

    // Parameter adjustments
    dint_dim1 = nint;

    ndm = iap->ndm;

    f[0] = 0.;
    f[1] = -par[12];

    for (int64 i = 0; i < ndm; ++i) {
        if (global_rotations.nrtn[i] == 0) {
            f[0] += u[i]*upold[i];
        }
        f[1] += u[ndm + i]*u[ndm + i];
    }

    if (ijac == 0) {
        return 0;
    }

    jtmp = NPARX;
    nn = ndim + jtmp;
    for (int64 i = 0; i < nint; ++i) {
        for (int64 j = 0; j < nn; ++j) {
            ARRAY2D(dint, i, j) = 0.;
        }
    }

    for (int64 i = 0; i < ndm; ++i) {
        if (global_rotations.nrtn[i] == 0) {
            ARRAY2D(dint, 0, i) = upold[i];
        } else {
            ARRAY2D(dint, 0, i) = 0.;
        }
        ARRAY2D(dint, 1, ndm + i) = u[ndm + i]*2.;
    }

    ARRAY2D(dint, 1, ndim + 12) = -1.;

    return 0;
}

int32
stpnpd(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr, int64 *ncolrs,
       double *rlcur, double *rldot, int64 *ndxloc, double *ups, double *udotps, double *upoldp,
       double *tm, double *dtm, int64 *nodir, double *thl, double *thu) {
    int64 ups_dim1;
    int64 udotps_dim1;
    int64 upoldp_dim1;

    int64 ndim;
    double temp[7];
    int64 nfpr;
    int64 nfpr1;
    int64 ntpl1;
    int64 nrsp1;
    int64 ntot1;
    int64 found;
    int64 icprs[NPARX];
    int64 nparr;
    int64 k1;
    int64 k2;
    int64 nskip1;

    int64 ibr;
    int64 ndm;
    int64 irs;
    int64 lab1;
    int64 nar1;
    int64 itp1;
    int64 isw1;

    (void)thu;
    (void)thl;
    (void)dtm;
    (void)upoldp;

    // Generates starting data for the 2-parameter continuation of
    // period-doubling bifurcations on a branch of periodic solutions.

    // Local

    // Parameter adjustments
    upoldp_dim1 = *ndxloc;
    udotps_dim1 = *ndxloc;
    ups_dim1 = *ndxloc;

    ndim = iap->ndim;
    irs = iap->irs;
    ndm = iap->ndm;
    nfpr = iap->nfpr;
    ibr = iap->ibr;

    findlb(iap, rap, irs, &nfpr1, &found);
    fscanf(fp3, "%ld", &ibr);
    fscanf(fp3, "%ld", &ntot1);
    fscanf(fp3, "%ld", &itp1);
    fscanf(fp3, "%ld", &lab1);
    fscanf(fp3, "%ld", &nfpr1);
    fscanf(fp3, "%ld", &isw1);
    fscanf(fp3, "%ld", &ntpl1);
    fscanf(fp3, "%ld", &nar1);
    fscanf(fp3, "%ld", &nskip1);
    fscanf(fp3, "%ld", &(*ntsr));
    fscanf(fp3, "%ld", &(*ncolrs));
    fscanf(fp3, "%ld", &nparr);
    iap->ibr = ibr;
    nrsp1 = *ntsr + 1;

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim;
            k2 = k1 + ndm - 1;
            fscanf(fp3, "%lf", &temp[i]);
            for (int64 k = k1; k <= k2; ++k) {
                fscanf(fp3, "%lf", &ARRAY2D(ups, j, k));
            }
        }
        tm[j] = temp[0];
    }
    fscanf(fp3, "%lf", &tm[-1 + nrsp1]);
    for (int64 k = 0; k < ndm; ++k) {
        fscanf(fp3, "%lf", &ARRAY2D(ups, *ntsr, k));
    }

    fscanf(fp3, "%ld", icprs);
    fscanf(fp3, "%ld", &icprs[1]);
    fscanf(fp3, "%lf", rldot);
    fscanf(fp3, "%lf", &rldot[1]);

    // Read U-dot (derivative with respect to arclength).

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim;
            k2 = k1 + ndm - 1;
            for (int64 k = k1; k <= k2; ++k) {
                fscanf(fp3, "%lf", &ARRAY2D(udotps, j, k));
            }
        }
    }

    for (int64 k = 0; k < ndm; ++k) {
        fscanf(fp3, "%lf", &ARRAY2D(udotps, *ntsr, k));
    }

    // Read the parameter values.

    if (nparr > NPARX) {
        nparr = NPARX;
        printf("Warning : NPARX too small for restart data\n");
        printf("PAR(i) set to zero, fot i > %3ld\n", nparr);
    }
    for (int64 i = 0; i < nparr; ++i) {
        fscanf(fp3, "%lf", &par[i]);
    }

    // Complement starting data
    par[12] = 0.;
    rldot[2] = 0.;
    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim + ndm;
            k2 = (i + 1)*ndim - 1;
            for (int64 k = k1; k <= k2; ++k) {
                ARRAY2D(ups, j, k) = 0.;
                ARRAY2D(udotps, j, k) = 0.;
            }
        }
    }
    for (int64 k = ndm; k < ndim; ++k) {
        ARRAY2D(ups, *ntsr, k) = 0.;
        ARRAY2D(udotps, *ntsr, k) = 0.;
    }

    for (int64 i = 0; i < nfpr; ++i) {
        rlcur[i] = par[icp[i]];
    }

    *nodir = 0;

    return 0;
}

/* ----------------------------------------------------------------------- */
/*       Subroutines for the Continuation of Torus Bifurcations */
/* ----------------------------------------------------------------------- */

int32
fntr(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 nfpr;
    double rtmp;
    double ep;
    int64 ndm;
    double umx;

    // Generates the equations for the 2-parameter continuation of
    // torus bifurcations.

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;
    nfpr = iap->nfpr;

    // Generate the function.

    fftr(iap, rap, ndim, u, uold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        fftr(iap, rap, ndim, global_scratch.uu1, uold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        fftr(iap, rap, ndim, global_scratch.uu2, uold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        fftr(iap, rap, ndim, u, uold, icp, par, global_scratch.ff1, ndm, global_scratch.dfu,
             global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdp, j, icp[i]) = (global_scratch.ff1[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }

    return 0;
}

int32
fftr(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double period;
    int64 ndm2;

    (void)ndim;

    // Parameter adjustments
    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    period = par[10];
    funi(iap, rap, ndm, u, uold, icp, par, 1, f, dfdu, dfdp);

    ndm2 = ndm*2;
    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = 0.;
        f[ndm2 + i] = 0.;
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm + i] += ARRAY2D(dfdu, i, j)*u[ndm + j];
            f[ndm2 + i] += ARRAY2D(dfdu, i, j)*u[ndm2 + j];
        }
        f[ndm + i] = period*f[ndm + i];
        f[ndm2 + i] = period*f[ndm2 + i];
        f[i] = period*f[i];
    }

    return 0;
}

int32
bctr(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nbc, double *u0,
     double *u1, double *f, int64 ijac, double *dbc) {
    int64 dbc_dim1;

    int64 jtmp;
    double theta;
    double cs;
    int64 nn;
    double ss;
    int64 ndm;
    int64 ndm2;

    (void)icp;
    (void)rap;

    // Parameter adjustments
    dbc_dim1 = nbc;

    ndm = iap->ndm;

    ndm2 = ndm << 1;
    theta = par[11];

    ss = sin(theta);
    cs = cos(theta);

    for (int64 i = 0; i < ndm; ++i) {
        f[i] = u0[i] - u1[i];
        f[ndm + i] = u1[ndm + i] - cs*u0[ndm + i] + ss*u0[ndm2 + i];
        f[ndm2 + i] = u1[ndm2 + i] - cs*u0[ndm2 + i] - ss*u0[ndm + i];
    }

    // Rotations
    if (global_rotations.irtn != 0) {
        for (int64 i = 0; i < ndm; ++i) {
            if (global_rotations.nrtn[i] != 0) {
                f[i] += par[18]*(double)global_rotations.nrtn[i];
            }
        }
    }

    if (ijac == 0) {
        return 0;
    }

    jtmp = NPARX;
    nn = (ndim*2) + jtmp;
    for (int64 i = 0; i < nbc; ++i) {
        for (int64 j = 0; j < nn; ++j) {
            ARRAY2D(dbc, i, j) = 0.;
        }
    }

    for (int64 i = 0; i < ndm; ++i) {
        ARRAY2D(dbc, i, i) = 1.;
        ARRAY2D(dbc, i, (ndim + i)) = -1.;
        ARRAY2D(dbc, ndm + i, (ndm + i)) = -cs;
        ARRAY2D(dbc, ndm + i, (ndm2 + i)) = ss;
        ARRAY2D(dbc, ndm + i, (ndim + ndm + i)) = 1.;
        ARRAY2D(dbc, ndm + i, ((ndim*2) + 11)) = cs*u0[ndm2 + i] + ss*u0[ndm + i];
        ARRAY2D(dbc, ndm2 + i, (ndm + i)) = -ss;
        ARRAY2D(dbc, ndm2 + i, (ndm2 + i)) = -cs;
        ARRAY2D(dbc, ndm2 + i, (ndim + ndm2 + i)) = 1.;
        ARRAY2D(dbc, ndm2 + i, ((ndim*2) + 11)) = ss*u0[ndm2 + i] - cs*u0[ndm + i];
    }

    return 0;
}

int32
ictr(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, double *u,
     double *uold, double *udot, double *upold, double *f, int64 ijac, double *dint) {
    int64 dint_dim1;

    int64 jtmp;
    int64 nn;
    int64 ndm;
    int64 ndm2;

    (void)udot;
    (void)uold;
    (void)icp;
    (void)rap;

    // Parameter adjustments
    dint_dim1 = nint;

    ndm = iap->ndm;
    ndm2 = ndm*2;

    f[0] = 0.;
    f[1] = 0.;
    f[2] = -par[12];

    for (int64 i = 0; i < ndm; ++i) {
        if (global_rotations.nrtn[i] == 0) {
            f[0] += u[i]*upold[i];
        }
        f[1] = f[1] + u[ndm + i]*u[ndm2 + i] - u[ndm2 + i]*u[ndm + i];
        f[2] = f[2] + u[ndm + i]*u[ndm + i] + u[ndm2 + i]*u[ndm2 + i];
    }

    if (ijac == 0) {
        return 0;
    }

    jtmp = NPARX;
    nn = ndim + jtmp;
    for (int64 i = 0; i < nint; ++i) {
        for (int64 j = 0; j < nn; ++j) {
            ARRAY2D(dint, i, j) = 0.;
        }
    }

    for (int64 i = 0; i < ndm; ++i) {
        if (global_rotations.nrtn[i] == 0) {
            ARRAY2D(dint, 0, i) = upold[i];
        } else {
            ARRAY2D(dint, 0, i) = 0.;
        }
        ARRAY2D(dint, 1, ndm + i) = u[ndm2 + i];
        ARRAY2D(dint, 1, ndm2 + i) = -u[ndm + i];
        ARRAY2D(dint, 2, ndm + i) = u[ndm + i]*2;
        ARRAY2D(dint, 2, ndm2 + i) = u[ndm2 + i]*2;
    }

    ARRAY2D(dint, 2, ndim + 12) = -1.;

    return 0;
}

int32
stpntr(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr, int64 *ncolrs,
       double *rlcur, double *rldot, int64 *ndxloc, double *ups, double *udotps, double *upoldp,
       double *tm, double *dtm, int64 *nodir, double *thl, double *thu) {
    int64 ups_dim1;
    int64 udotps_dim1;

    int64 ndim;
    double temp[7];
    int64 nfpr;
    int64 nfpr1;
    int64 ntpl1;
    int64 nrsp1;
    int64 ntot1;
    int64 found;
    int64 icprs[NPARX];
    int64 nparr;
    int64 k1;
    int64 k2;
    int64 k3;
    int64 nskip1;

    int64 ibr;
    int64 ndm;
    int64 k2p1;
    int64 irs;
    int64 lab1;
    int64 nar1;
    int64 itp1;
    int64 isw1;

    (void)thu;
    (void)thl;
    (void)dtm;
    (void)upoldp;

    // Generates starting data for the 2-parameter continuation of torus
    // bifurcations.

    // Local

    // Parameter adjustments
    udotps_dim1 = *ndxloc;
    ups_dim1 = *ndxloc;

    ndim = iap->ndim;
    irs = iap->irs;
    ndm = iap->ndm;
    nfpr = iap->nfpr;
    ibr = iap->ibr;

    findlb(iap, rap, irs, &nfpr1, &found);
    fscanf(fp3, "%ld", &ibr);
    fscanf(fp3, "%ld", &ntot1);
    fscanf(fp3, "%ld", &itp1);
    fscanf(fp3, "%ld", &lab1);
    fscanf(fp3, "%ld", &nfpr1);
    fscanf(fp3, "%ld", &isw1);
    fscanf(fp3, "%ld", &ntpl1);
    fscanf(fp3, "%ld", &nar1);
    fscanf(fp3, "%ld", &nskip1);
    fscanf(fp3, "%ld", &(*ntsr));
    fscanf(fp3, "%ld", &(*ncolrs));
    fscanf(fp3, "%ld", &nparr);
    iap->ibr = ibr;
    nrsp1 = *ntsr + 1;

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim;
            k2 = k1 + ndm - 1;
            fscanf(fp3, "%lf", &temp[i]);
            for (int64 k = k1; k <= k2; ++k) {
                fscanf(fp3, "%lf", &ARRAY2D(ups, j, k));
            }
            k2p1 = k2 + 1;
            k3 = k2 + ndm;
            for (int64 k = k2p1; k <= k3; ++k) {
                ARRAY2D(ups, j, k) = sin(temp[i])*(double)1e-4;
                ARRAY2D(ups, j, (k + ndm)) = cos(temp[i])*(double)1e-4;
            }
        }
        tm[j] = temp[0];
    }

    fscanf(fp3, "%lf", &tm[*ntsr]);
    for (int64 k = 0; k < ndm; ++k) {
        fscanf(fp3, "%lf", &ARRAY2D(ups, *ntsr, k));
    }
    for (int64 i = 0; i < ndm; ++i) {
        ARRAY2D(ups, *ntsr, (ndm + i)) = 0.;
        ARRAY2D(ups, *ntsr, ((ndm*2) + i)) = 0.;
    }

    fscanf(fp3, "%ld", icprs);
    fscanf(fp3, "%ld", &icprs[1]);
    fscanf(fp3, "%lf", rldot);
    fscanf(fp3, "%lf", &rldot[1]);
    rldot[2] = 0.;
    rldot[3] = 0.;

    // Read the direction vector and initialize the starting direction.

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim;
            k2 = k1 + ndm - 1;
            for (int64 k = k1; k <= k2; ++k) {
                fscanf(fp3, "%lf", &ARRAY2D(udotps, j, k));
            }
            k2p1 = k2 + 1;
            k3 = k2 + ndm;
            for (int64 k = k2p1; k <= k2; ++k) {
                ARRAY2D(udotps, j, k) = 0.;
                ARRAY2D(udotps, j, (k + ndm)) = 0.;
            }
        }
    }

    for (int64 k = 0; k < ndm; ++k) {
        fscanf(fp3, "%lf", &ARRAY2D(udotps, *ntsr, k));
    }
    for (int64 i = 0; i < ndm; ++i) {
        ARRAY2D(udotps, *ntsr, (ndm + i)) = 0.;
        ARRAY2D(udotps, *ntsr, ((ndm*2) + i)) = 0.;
    }

    // Read the parameter values.

    if (nparr > NPARX) {
        nparr = NPARX;
        printf("Warning : NPARX too small for restart data\n");
        printf("PAR(i) set to zero, fot i > %3ld\n", nparr);
    }
    for (int64 i = 0; i < nparr; ++i) {
        fscanf(fp3, "%lf", &par[i]);
    }

    par[12] = 0.;

    for (int64 i = 0; i < nfpr; ++i) {
        rlcur[i] = par[icp[i]];
    }

    *nodir = 0;

    return 0;
}

/* ----------------------------------------------------------------------- */
/*        Subroutines for Optimization of Periodic Solutions */
/* ----------------------------------------------------------------------- */

int32
fnpo(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 nfpr;
    double rtmp;
    double *upold, ep, period;
    int64 ndm;
    double umx;

    upold = xmalloc(sizeof(*upold)*(usize)(iap->ndim));

    // Generates the equations for periodic optimization problems.

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;
    nfpr = iap->nfpr;

    // Generate F(UOLD)

    funi(iap, rap, ndm, uold, uold, icp, par, 0, upold, dfdu, dfdp);
    period = par[10];
    for (int64 i = 0; i < ndm; ++i) {
        upold[i] = period*upold[i];
    }

    // Generate the function.

    ffpo(iap, rap, ndim, u, uold, upold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        free(upold);
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        ffpo(iap, rap, ndim, global_scratch.uu1, uold, upold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        ffpo(iap, rap, ndim, global_scratch.uu2, uold, upold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        ffpo(iap, rap, ndim, u, uold, upold, icp, par, global_scratch.ff1, ndm, global_scratch.dfu,
             global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdp, j, icp[i]) = (global_scratch.ff1[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }

    free(upold);
    return 0;
}

int32
ffpo(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, double *upold, int64 *icp,
     double *par, double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double gamma, rkappa, period, dfp[NPARX], *dfu, fop;

    (void)ndim;

    dfu = xmalloc(sizeof(*dfu)*(usize)(iap->ndim));
    // Local

    // Parameter adjustments
    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    period = par[10];
    rkappa = par[12];
    gamma = par[13];

    for (int64 i = 0; i < ndm; ++i) {
        for (int64 j = 0; j < NPARX; ++j) {
            ARRAY2D(dfdp, i, j) = 0.;
        }
    }
    funi(iap, rap, ndm, u, uold, icp, par, 1, f, dfdu, dfdp);
    for (int64 i = 0; i < NPARX; ++i) {
        dfp[i] = 0.;
    }
    fopi(iap, rap, ndm, u, icp, par, 1, &fop, dfu, dfp);

    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = 0.;
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm + i] -= ARRAY2D(dfdu, j, i)*u[ndm + j];
        }
        f[i] = period*f[i];
        f[ndm + i] = period*f[ndm + i] + rkappa*upold[i] + gamma*dfu[i];
    }

    free(dfu);
    return 0;
}

int32
bcpo(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nbc, double *u0,
     double *u1, double *f, int64 ijac, double *dbc) {
    int64 dbc_dim1;

    int64 nfpr;
    int64 nbc0;

    (void)rap;

    /* Generates the boundary conditions for periodic optimization problems.
     */

    // Parameter adjustments
    dbc_dim1 = nbc;

    nfpr = iap->nfpr;

    for (int64 i = 0; i < nbc; ++i) {
        f[i] = u0[i] - u1[i];
    }

    // Rotations
    if (global_rotations.irtn != 0) {
        nbc0 = iap->nbc0;
        for (int64 i = 0; i < nbc0; ++i) {
            if (global_rotations.nrtn[i] != 0) {
                f[i] += par[18]*(double)global_rotations.nrtn[i];
            }
        }
    }

    if (ijac == 0) {
        return 0;
    }

    for (int64 i = 0; i < nbc; ++i) {
        for (int64 j = 0; j <= (ndim*2); ++j) {
            ARRAY2D(dbc, i, j) = 0.;
        }
        ARRAY2D(dbc, i, i) = 1.;
        ARRAY2D(dbc, i, (ndim + i)) = -1.;
        for (int64 j = 0; j < nfpr; ++j) {
            ARRAY2D(dbc, i, (ndim*2) + icp[j]) = 0.;
        }
    }
    return 0;
}

int32
icpo(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, double *u,
     double *uold, double *udot, double *upold, double *f, int64 ijac, double *dint) {
    int64 dint_dim1;

    int64 nfpr;
    double rtmp;
    double *f1, *f2, ep;
    int64 ndm;
    double *dnt;
    double umx;
    int64 nnt0;

    f1 = xmalloc(sizeof(*(f1))*(usize)(iap->nint));
    f2 = xmalloc(sizeof(*(f2))*(usize)(iap->nint));
    dnt = xmalloc(sizeof(*dnt)*(usize)((iap->nint)*(iap->ndim + NPARX)));

    // Generates integral conditions for periodic optimization problems.

    // Local

    // Parameter adjustments
    dint_dim1 = nint;

    ndm = iap->ndm;
    nnt0 = iap->nnt0;
    nfpr = iap->nfpr;

    // Generate the function.

    fipo(iap, rap, ndim, par, icp, nint, nnt0, u, uold, udot, upold, f, dnt, ndm,
         global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        free(f1);
        free(f2);
        free(dnt);
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        fipo(iap, rap, ndim, par, icp, nint, nnt0, global_scratch.uu1, uold, udot, upold, f1, dnt,
             ndm, global_scratch.dfu, global_scratch.dfp);
        fipo(iap, rap, ndim, par, icp, nint, nnt0, global_scratch.uu2, uold, udot, upold, f2, dnt,
             ndm, global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < nint; ++j) {
            ARRAY2D(dint, j, i) = (f2[j] - f1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        fipo(iap, rap, ndim, par, icp, nint, nnt0, u, uold, udot, upold, f1, dnt, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < nint; ++j) {
            ARRAY2D(dint, j, ndim + icp[i]) = (f1[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }
    free(f1);
    free(f2);
    free(dnt);

    return 0;
}

int32
fipo(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, int64 nnt0,
     double *u, double *uold, double *udot, double *upold, double *fi, double *dint, int64 ndmt,
     double *dfdu, double *dfdp) {
    int64 dint_dim1;
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 nfpr;
    int64 indx;
    double *f;
    int64 l;
    double dfp[NPARX], *dfu;
    int64 ndm;
    double fop;

    (void)dint;
    (void)udot;
    (void)ndim;

    f = xmalloc(sizeof(*f)*(usize)(iap->ndim));
    dfu = xmalloc(sizeof(*dfu)*(usize)(iap->ndim));
    // Local

    // Parameter adjustments
    dint_dim1 = nnt0;
    dfdp_dim1 = ndmt;
    dfdu_dim1 = ndmt;

    ndm = iap->ndm;
    nfpr = iap->nfpr;

    fi[0] = 0.;
    for (int64 i = 0; i < ndm; ++i) {
        if (global_rotations.nrtn[i] == 0) {
            fi[0] += u[i]*upold[i];
        }
    }

    for (int64 i = 0; i < NPARX; ++i) {
        dfp[i] = 0.;
    }
    fopi(iap, rap, ndm, u, icp, par, 2, &fop, dfu, dfp);
    fi[1] = par[9] - fop;

    // Computing 2nd power
    fi[2] = par[12]*par[12] + par[13]*par[13] - par[11];
    for (int64 i = 0; i < ndm; ++i) {
        // Computing 2nd power
        fi[2] += u[ndm + i]*u[ndm + i];
    }

    for (int64 i = 0; i < ndm; ++i) {
        for (int64 j = 0; j < NPARX; ++j) {
            ARRAY2D(dfdp, i, j) = 0.;
        }
    }
    funi(iap, rap, ndm, u, uold, icp, par, 2, f, dfdu, dfdp);

    for (l = 3; l < nint; ++l) {
        indx = icp[nfpr + l - 3];
        if (indx == 10) {
            fi[l] = -par[13]*dfp[indx] - par[indx + 20];
            for (int64 i = 0; i < ndm; ++i) {
                fi[l] += f[i]*u[ndm + i];
            }
        } else {
            fi[l] = -par[13]*dfp[indx] - par[indx + 20];
            for (int64 i = 0; i < ndm; ++i) {
                fi[l] += par[10]*ARRAY2D(dfdp, i, (indx))*u[ndm + i];
            }
        }
    }
    free(f);
    free(dfu);

    return 0;
}

int32
stpnpo(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr, int64 *ncolrs,
       double *rlcur, double *rldot, int64 *ndxloc, double *ups, double *udotps, double *upoldp,
       double *tm, double *dtm, int64 *nodir, double *thl, double *thu) {
    int64 ups_dim1;
    int64 udotps_dim1;
    int64 upoldp_dim1;

    int64 ndim;
    double temp[7];
    int64 nfpr;
    double dump;

    double dumu;
    int64 nfpr1;
    int64 ntpl1;
    int64 nrsp1;
    int64 ntot1;
    double *u;
    int64 found;
    int64 icprs[NPARX];
    int64 nparr;

    int64 k1;
    int64 k2;
    int64 nskip1;
    double fs;

    int64 ibr;
    int64 ndm;
    int64 irs;
    int64 lab1;
    int64 nar1;
    double rld1;
    double rld2;
    int64 itp1;
    int64 isw1;

    double *temporary_storage;
    int64 temporary_storage_dim1;
    usize size;

    (void)thu;
    (void)thl;
    (void)upoldp;
    (void)rldot;

    /* This is a little funky.  In the older version, upoldp was used for some
       temporary storage in a loop later on.  I wanted to get rid of that
       my adding a local varialbe.  Unfortunately, things are never that easy.
       The size of this has the same problems as computing the sizes in
       rsptbv.  The are various places the sizes are defined (fort.2 and fort.8)
       and you have to pick the maximum, multiplied by a constant (something
       like 4 to take into account the increase in size for certain
       calculations). So, that is why I use ndxloc here.  Also, iap->ncol MAY BE
       tool small, but I am not sure how to get value from the fort.8 file into
       here. */
    size = sizeof(*temporary_storage)*(usize)((*ndxloc)*(iap->ndim*iap->ncol));
    temporary_storage = xmalloc(size);
    u = xmalloc(sizeof(*u)*(usize)(iap->ndim));

    // Generates starting data for optimization of periodic solutions.

    // Local

    // Parameter adjustments
    upoldp_dim1 = *ndxloc;
    udotps_dim1 = *ndxloc;
    ups_dim1 = *ndxloc;
    temporary_storage_dim1 = *ndxloc;

    ndim = iap->ndim;
    irs = iap->irs;
    ndm = iap->ndm;
    nfpr = iap->nfpr;
    ibr = iap->ibr;

    findlb(iap, rap, irs, &nfpr1, &found);
    fscanf(fp3, "%ld", &ibr);
    fscanf(fp3, "%ld", &ntot1);
    fscanf(fp3, "%ld", &itp1);
    fscanf(fp3, "%ld", &lab1);
    fscanf(fp3, "%ld", &nfpr1);
    fscanf(fp3, "%ld", &isw1);
    fscanf(fp3, "%ld", &ntpl1);
    fscanf(fp3, "%ld", &nar1);
    fscanf(fp3, "%ld", &nskip1);
    fscanf(fp3, "%ld", &(*ntsr));
    fscanf(fp3, "%ld", &(*ncolrs));
    fscanf(fp3, "%ld", &nparr);
    iap->ibr = ibr;
    nrsp1 = *ntsr + 1;

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim;
            k2 = k1 + ndm - 1;
            fscanf(fp3, "%lf", &temp[i]);
            for (int64 k = k1; k <= k2; ++k) {
                fscanf(fp3, "%lf", &ARRAY2D(ups, j, k));
            }
        }
        tm[j] = temp[0];
    }
    fscanf(fp3, "%lf", &tm[*ntsr]);
    for (int64 k = 0; k < ndm; ++k) {
        fscanf(fp3, "%lf", &ARRAY2D(ups, *ntsr, k));
    }
    for (int64 j = 0; j < *ntsr; ++j) {
        dtm[j] = tm[j + 1] - tm[j];
    }

    fscanf(fp3, "%ld", icprs);
    fscanf(fp3, "%ld", &icprs[1]);
    fscanf(fp3, "%lf", &rld1);
    fscanf(fp3, "%lf", &rld2);

    // Read U-dot (derivative with respect to arclength).
    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim;
            k2 = k1 + ndm - 1;
            for (int64 k = k1; k <= k2; ++k) {
                fscanf(fp3, "%lf", &ARRAY2D(udotps, j, k));
            }
        }
    }
    for (int64 k = 0; k < ndm; ++k) {
        fscanf(fp3, "%lf", &ARRAY2D(udotps, *ntsr, k));
    }

    // Read the parameter values.
    if (nparr > NPARX) {
        nparr = NPARX;
        printf("Warning : NPARX too small for restart data\n");
        printf("PAR(i) set to zero, fot i > %3ld\n", nparr);
    }
    for (int64 i = 0; i < nparr; ++i) {
        fscanf(fp3, "%lf", &par[i]);
    }

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim;
            k2 = k1 + ndm - 1;
            for (int64 k = k1; k <= k2; ++k) {
                u[k - k1] = ARRAY2D(ups, j, k);
            }
            fopt(ndm, u, icp, par, 0, &fs, &dumu, &dump);
#define TEMPORARY_STORAGE
#ifdef TEMPORARY_STORAGE
            ARRAY2D(temporary_storage, j, k1) = fs;
#else
            ARRAY2D(upoldp, j, k1) = fs;
#endif
        }
    }
    for (int64 k = 0; k < ndm; ++k) {
        u[k] = ARRAY2D(ups, *ntsr, k);
    }
    fopt(ndm, u, icp, par, 0, &fs, &dumu, &dump);
#ifdef TEMPORARY_STORAGE
    temporary_storage[*ntsr] = fs;
    par[9] = rintg(iap, ndxloc, 1, temporary_storage, dtm);
#else
    upoldp[*ntsr] = fs;
    par[9] = rintg(iap, ndxloc, 1, upoldp, dtm);
#endif

    // Complement starting data

    for (int64 i = 11; i < NPARX; ++i) {
        par[i] = 0.;
    }

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim + ndm;
            k2 = (i + 1)*ndim - 1;
            for (int64 k = k1; k <= k2; ++k) {
                ARRAY2D(ups, j, k) = 0.;
            }
        }
    }
    for (int64 k = ndm; k < ndim; ++k) {
        ARRAY2D(ups, *ntsr, k) = 0.;
    }

    for (int64 i = 0; i < nfpr; ++i) {
        rlcur[i] = par[icp[i]];
    }

    *nodir = 1;

    free(u);
    free(temporary_storage);
    return 0;
}

/* ----------------------------------------------------------------------- */
/*        Subroutines for the Continuation of Folds for BVP. */
/* ----------------------------------------------------------------------- */

int32
fnbl(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 nfpr;
    double rtmp;
    double ep;
    int64 ndm;
    double umx;

    // Generates the equations for the 2-parameter continuation
    // of folds (BVP).

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;
    nfpr = iap->nfpr;

    // Generate the function.

    ffbl(iap, rap, ndim, u, uold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

    if (ijac == 0) {
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        ffbl(iap, rap, ndim, global_scratch.uu1, uold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        ffbl(iap, rap, ndim, global_scratch.uu2, uold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        ffbl(iap, rap, ndim, u, uold, icp, par, global_scratch.ff1, ndm, global_scratch.dfu,
             global_scratch.dfp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdp, j, icp[i]) = (global_scratch.ff1[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }

    return 0;
}

int32
ffbl(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 nfpr;
    int64 nfpx;

    (void)ndim;

    // Parameter adjustments
    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    nfpr = iap->nfpr;

    funi(iap, rap, ndm, u, uold, icp, par, 2, f, dfdu, dfdp);

    nfpx = nfpr / 2 - 1;
    for (int64 i = 0; i < ndm; ++i) {
        f[ndm + i] = 0.;
        for (int64 j = 0; j < ndm; ++j) {
            f[ndm + i] += ARRAY2D(dfdu, i, j)*u[ndm + j];
        }
        if (nfpx > 0) {
            for (int64 j = 0; j < nfpx; ++j) {
                f[ndm + i] += ARRAY2D(dfdp, i, icp[j + 1])*par[icp[nfpr - nfpx + j]];
            }
        }
    }

    return 0;
}

int32
bcbl(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nbc, double *u0,
     double *u1, double *f, int64 ijac, double *dbc) {
    int64 dbc_dim1;

    int64 nfpr;
    double rtmp;
    double ep, *ff1, *ff2, *uu1, *uu2, *dfu, umx;
    int64 nbc0;

    ff1 = xmalloc(sizeof(*(ff1))*(usize)(iap->nbc));
    ff2 = xmalloc(sizeof(*(ff2))*(usize)(iap->nbc));
    uu1 = xmalloc(sizeof(*(uu1))*(usize)(iap->ndim));
    uu2 = xmalloc(sizeof(*(uu2))*(usize)(iap->ndim));
    dfu = xmalloc(sizeof(*dfu)*(usize)((iap->nbc)*(2*iap->ndim + NPARX)));

    // Generates the boundary conditions for the 2-parameter continuation
    // of folds (BVP).

    // Local

    // Parameter adjustments
    dbc_dim1 = nbc;

    nbc0 = iap->nbc0;
    nfpr = iap->nfpr;

    // Generate the function.

    fbbl(iap, rap, ndim, par, icp, nbc, nbc0, u0, u1, f, dfu);

    if (ijac == 0) {
        free(ff1);
        free(ff2);
        free(uu1);
        free(uu2);
        free(dfu);
        return 0;
    }

    // Derivatives with respect to U0.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u0[i]) > umx) {
            umx = fabs(u0[i]);
        }
    }
    rtmp = HMACH;
    ep = rtmp*(umx + 1);
    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            uu1[j] = u0[j];
            uu2[j] = u0[j];
        }
        uu1[i] -= ep;
        uu2[i] += ep;
        fbbl(iap, rap, ndim, par, icp, nbc, nbc0, uu1, u1, ff1, dfu);
        fbbl(iap, rap, ndim, par, icp, nbc, nbc0, uu2, u1, ff2, dfu);
        for (int64 j = 0; j < nbc; ++j) {
            ARRAY2D(dbc, j, i) = (ff2[j] - ff1[j]) / (ep*2);
        }
    }

    // Derivatives with respect to U1.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u1[i]) > umx) {
            umx = fabs(u1[i]);
        }
    }
    rtmp = HMACH;
    ep = rtmp*(umx + 1);
    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            uu1[j] = u1[j];
            uu2[j] = u1[j];
        }
        uu1[i] -= ep;
        uu2[i] += ep;
        fbbl(iap, rap, ndim, par, icp, nbc, nbc0, u0, uu1, ff1, dfu);
        fbbl(iap, rap, ndim, par, icp, nbc, nbc0, u0, uu2, ff2, dfu);
        for (int64 j = 0; j < nbc; ++j) {
            ARRAY2D(dbc, j, (ndim + i)) = (ff2[j] - ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        fbbl(iap, rap, ndim, par, icp, nbc, nbc0, u0, u1, ff2, dfu);
        for (int64 j = 0; j < nbc; ++j) {
            ARRAY2D(dbc, j, (ndim*2) + icp[i]) = (ff2[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }
    free(ff1);
    free(ff2);
    free(uu1);
    free(uu2);
    free(dfu);

    return 0;
}

int32
fbbl(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nbc, int64 nbc0,
     double *u0, double *u1, double *f, double *dbc) {
    int64 dbc_dim1;

    int64 nfpr;
    int64 nfpx;
    int64 ndm;

    (void)nbc;

    // Parameter adjustments
    dbc_dim1 = nbc0;

    ndm = iap->ndm;
    nfpr = iap->nfpr;

    nfpx = nfpr / 2 - 1;
    bcni(iap, rap, ndm, par, icp, nbc0, u0, u1, f, 2, dbc);
    for (int64 i = 0; i < nbc0; ++i) {
        f[nbc0 + i] = 0.;
        for (int64 j = 0; j < ndm; ++j) {
            f[nbc0 + i] += ARRAY2D(dbc, i, j)*u0[ndm + j];
            f[nbc0 + i] += ARRAY2D(dbc, i, (ndm + j))*u1[ndm + j];
        }
        if (nfpx != 0) {
            for (int64 j = 0; j < nfpx; ++j) {
                f[nbc0 + i] += ARRAY2D(dbc, i, ndim + icp[j + 1])*par[icp[nfpr - nfpx + j]];
            }
        }
    }

    return 0;
}

int32
icbl(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, double *u,
     double *uold, double *udot, double *upold, double *f, int64 ijac, double *dint) {
    int64 dint_dim1;

    int64 nfpr;
    double rtmp;
    double ep, *ff1, *ff2, *uu1, *uu2, *dfu, umx;
    int64 nnt0;

    ff1 = xmalloc(sizeof(*(ff1))*(usize)(iap->nint));
    ff2 = xmalloc(sizeof(*(ff2))*(usize)(iap->nint));
    uu1 = xmalloc(sizeof(*(uu1))*(usize)(iap->ndim));
    uu2 = xmalloc(sizeof(*(uu2))*(usize)(iap->ndim));
    dfu = xmalloc(sizeof(*dfu)*(usize)((iap->ndim)*(iap->ndim + NPARX)));

    // Generates integral conditions for the 2-parameter continuation of
    // folds (BVP).

    // Local

    // Parameter adjustments
    dint_dim1 = nint;

    nnt0 = iap->nnt0;
    nfpr = iap->nfpr;

    // Generate the function.

    fibl(iap, rap, ndim, par, icp, nint, nnt0, u, uold, udot, upold, f, dfu);

    if (ijac == 0) {
        free(ff1);
        free(ff2);
        free(uu1);
        free(uu2);
        free(dfu);
        return 0;
    }

    // Generate the Jacobian.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            uu1[j] = u[j];
            uu2[j] = u[j];
        }
        uu1[i] -= ep;
        uu2[i] += ep;
        fibl(iap, rap, ndim, par, icp, nint, nnt0, uu1, uold, udot, upold, ff1, dfu);
        fibl(iap, rap, ndim, par, icp, nint, nnt0, uu2, uold, udot, upold, ff2, dfu);
        for (int64 j = 0; j < nint; ++j) {
            ARRAY2D(dint, j, i) = (ff2[j] - ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        fibl(iap, rap, ndim, par, icp, nint, nnt0, u, uold, udot, upold, ff1, dfu);
        for (int64 j = 0; j < nint; ++j) {
            ARRAY2D(dint, j, ndim + icp[i]) = (ff1[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }

    free(ff1);
    free(ff2);
    free(uu1);
    free(uu2);
    free(dfu);
    return 0;
}

int32
fibl(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, int64 nnt0,
     double *u, double *uold, double *udot, double *upold, double *f, double *dint) {
    int64 dint_dim1;

    int64 nfpr;
    int64 nfpx = 0;
    int64 ndm;

    (void)ndim;

    // Parameter adjustments
    dint_dim1 = nnt0;

    ndm = iap->ndm;
    nfpr = iap->nfpr;

    if (nnt0 > 0) {
        nfpx = nfpr / 2 - 1;
        icni(iap, rap, ndm, par, icp, nnt0, u, uold, udot, upold, f, 2, dint);
        for (int64 i = 0; i < nnt0; ++i) {
            f[nnt0 + i] = 0.;
            for (int64 j = 0; j < ndm; ++j) {
                f[nnt0 + i] += ARRAY2D(dint, i, j)*u[ndm + j];
            }
            if (nfpx != 0) {
                for (int64 j = 0; j < nfpx; ++j) {
                    f[nnt0 + i] += ARRAY2D(dint, i, ndm + icp[j + 1])*par[icp[nfpr - nfpx + j]];
                }
            }
        }
    }

    // Note that PAR(11+NFPR/2) is used to keep the norm of the null vector
    f[-1 + nint] = -par[-1 + nfpr / 2 + 11];
    for (int64 i = 0; i < ndm; ++i) {
        f[-1 + nint] += u[ndm + i]*u[ndm + i];
    }
    if (nfpx != 0) {
        for (int64 i = 0; i < nfpx; ++i) {
            // Computing 2nd power
            f[-1 + nint] += par[icp[nfpr - nfpx + i]]*par[icp[nfpr - nfpx + i]];
        }
    }

    return 0;
}

int32
stpnbl(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr, int64 *ncolrs,
       double *rlcur, double *rldot, int64 *ndxloc, double *ups, double *udotps, double *upoldp,
       double *tm, double *dtm, int64 *nodir, double *thl, double *thu) {
    int64 ups_dim1;
    int64 udotps_dim1;

    int64 ndim;
    double temp[7];
    int64 nfpr;
    int64 nfpx;
    int64 nfpr0;
    int64 nfpr1;
    int64 ntpl1;
    int64 nrsp1;
    int64 ntot1;
    int64 found;
    int64 icprs[NPARX];
    int64 nparr;
    int64 k1;
    int64 k2;
    int64 nskip1;

    int64 ibr;
    int64 ndm;
    int64 irs;
    int64 lab1;
    int64 nar1;
    int64 itp1;
    int64 isw1;

    (void)thu;
    (void)thl;
    (void)dtm;
    (void)upoldp;
    (void)udotps;

    // Generates starting data for the 2-parameter continuation of folds.
    // (BVP).

    // Local

    // Parameter adjustments
    udotps_dim1 = *ndxloc;
    ups_dim1 = *ndxloc;

    ndim = iap->ndim;
    irs = iap->irs;
    ndm = iap->ndm;
    nfpr = iap->nfpr;
    ibr = iap->ibr;

    findlb(iap, rap, irs, &nfpr1, &found);
    fscanf(fp3, "%ld", &ibr);
    fscanf(fp3, "%ld", &ntot1);
    fscanf(fp3, "%ld", &itp1);
    fscanf(fp3, "%ld", &lab1);
    fscanf(fp3, "%ld", &nfpr1);
    fscanf(fp3, "%ld", &isw1);
    fscanf(fp3, "%ld", &ntpl1);
    fscanf(fp3, "%ld", &nar1);
    fscanf(fp3, "%ld", &nskip1);
    fscanf(fp3, "%ld", &(*ntsr));
    fscanf(fp3, "%ld", &(*ncolrs));
    fscanf(fp3, "%ld", &nparr);
    iap->ibr = ibr;
    nrsp1 = *ntsr + 1;

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim;
            k2 = k1 + ndm - 1;
            fscanf(fp3, "%lf", &temp[i]);
            for (int64 k = k1; k <= k2; ++k) {
                fscanf(fp3, "%lf", &ARRAY2D(ups, j, k));
            }
        }
        tm[j] = temp[0];
    }
    fscanf(fp3, "%lf", &tm[*ntsr]);
    for (int64 k = 0; k < ndm; ++k) {
        fscanf(fp3, "%lf", &ARRAY2D(ups, *ntsr, k));
    }

    nfpr0 = nfpr / 2;
    fscanf(fp3, "%ld", icprs);
    for (int64 i = 0; i < nfpr0; ++i) {
        fscanf(fp3, "%lf", &rldot[i]);
    }

    // Read U-dot (Derivative with respect to arclength).

    for (int64 j = 0; j < *ntsr; ++j) {
        for (int64 i = 0; i < *ncolrs; ++i) {
            k1 = i*ndim + ndm;
            k2 = (i + 1)*ndim - 1;
            for (int64 k = k1; k <= k2; ++k) {
                fscanf(fp3, "%lf", &ARRAY2D(ups, j, k));
            }
        }
    }
    for (int64 k = ndm; k < ndim; ++k) {
        fscanf(fp3, "%lf", &ARRAY2D(ups, *ntsr, k));
    }

    // Read the parameter values.

    if (nparr > NPARX) {
        nparr = NPARX;
        printf("Warning : NPARX too small for restart data\n");
        printf("PAR(i) set to zero, for i > %3ld\n", nparr);
    }
    for (int64 i = 0; i < nparr; ++i) {
        fscanf(fp3, "%lf", &par[i]);
    }

    nfpx = nfpr / 2 - 1;
    if (nfpx > 0) {
        for (int64 i = 0; i < nfpx; ++i) {
            par[icp[nfpr0 + 1 + i]] = rldot[i + 1];
        }
    }
    // Initialize the norm of the null vector
    par[-1 + nfpr / 2 + 11] = (double)0.;

    for (int64 i = 0; i < nfpr; ++i) {
        rlcur[i] = par[icp[i]];
    }

    *nodir = 1;

    return 0;
}

/* ----------------------------------------------------------------------- */
/*          Routines for Interface with User Supplied Routines */
/*  (To generate Jacobian by differencing, if not supplied analytically) */
/* ----------------------------------------------------------------------- */

int32
funi(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double *u1zz;
    double *u2zz;

    int64 nfpr;
    double rtmp;
    double ep;
    int64 jac;
    int64 ijc;
    double umx, *f1zz, *f2zz;

    (void)uold;
    (void)rap;

    u1zz = xmalloc(sizeof(*(u1zz))*(usize)(iap->ndim));
    u2zz = xmalloc(sizeof(*(u2zz))*(usize)(iap->ndim));
    f1zz = xmalloc(sizeof(*(f1zz))*(usize)(iap->ndim));
    f2zz = xmalloc(sizeof(*(f2zz))*(usize)(iap->ndim));

    // Interface subroutine to user supplied FUNC.

    // Local

    // Parameter adjustments
    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    jac = iap->jac;
    nfpr = iap->nfpr;

    // Generate the function.

    if (jac == 0) {
        ijc = 0;
    } else {
        ijc = ijac;
    }
    func(ndim, u, icp, par, ijc, f, dfdu, dfdp);

    if (jac == 1 || ijac == 0) {
        free(u1zz);
        free(u2zz);
        free(f1zz);
        free(f2zz);
        return 0;
    }

    // Generate the Jacobian by differencing.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            u1zz[j] = u[j];
            u2zz[j] = u[j];
        }
        u1zz[i] -= ep;
        u2zz[i] += ep;
        func(ndim, u1zz, icp, par, 0, f1zz, dfdu, dfdp);
        func(ndim, u2zz, icp, par, 0, f2zz, dfdu, dfdp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (f2zz[j] - f1zz[j]) / (ep*2);
        }
    }

    if (ijac == 1) {
        free(u1zz);
        free(u2zz);
        free(f1zz);
        free(f2zz);
        return 0;
    }

    for (int64 i = 0; i < nfpr; ++i) {
        rtmp = HMACH;
        ep = rtmp*(fabs(par[icp[i]]) + 1);
        par[icp[i]] += ep;
        func(ndim, u, icp, par, 0, f1zz, dfdu, dfdp);
        for (int64 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdp, j, icp[i]) = (f1zz[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }

    free(u1zz);
    free(u2zz);
    free(f1zz);
    free(f2zz);
    return 0;
}

int32
bcni(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nbc, double *u0,
     double *u1, double *f, int64 ijac, double *dbc) {
    int64 dbc_dim1;

    double *u1zz;
    double *u2zz;
    int64 nfpr;
    double rtmp;
    double ep;
    int64 jac;
    int64 ijc;
    double umx, *f1zz, *f2zz;

    (void)rap;

    u1zz = xmalloc(sizeof(*(u1zz))*(usize)(iap->ndim));
    u2zz = xmalloc(sizeof(*(u2zz))*(usize)(iap->ndim));
    f1zz = xmalloc(sizeof(*(f1zz))*(usize)(iap->nbc));
    f2zz = xmalloc(sizeof(*(f2zz))*(usize)(iap->nbc));

    // Interface subroutine to the user supplied BCND.

    // Local

    // Parameter adjustments
    dbc_dim1 = nbc;

    jac = iap->jac;
    nfpr = iap->nfpr;

    // Generate the function.

    if (jac == 0) {
        ijc = 0;
    } else {
        ijc = ijac;
    }
    bcnd(ndim, par, icp, nbc, u0, u1, ijc, f, dbc);

    if (jac == 1 || ijac == 0) {
        return 0;
    }

    // Generate the Jacobian by differencing.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u0[i]) > umx) {
            umx = fabs(u0[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            u1zz[j] = u0[j];
            u2zz[j] = u0[j];
        }
        u1zz[i] -= ep;
        u2zz[i] += ep;
        bcnd(ndim, par, icp, nbc, u1zz, u1, 0, f1zz, dbc);
        bcnd(ndim, par, icp, nbc, u2zz, u1, 0, f2zz, dbc);
        for (int64 j = 0; j < nbc; ++j) {
            ARRAY2D(dbc, j, i) = (f2zz[j] - f1zz[j]) / (ep*2);
        }
    }

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u1[i]) > umx) {
            umx = fabs(u1[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            u1zz[j] = u1[j];
            u2zz[j] = u1[j];
        }
        u1zz[i] -= ep;
        u2zz[i] += ep;
        bcnd(ndim, par, icp, nbc, u0, u1zz, 0, f1zz, dbc);
        bcnd(ndim, par, icp, nbc, u0, u2zz, 0, f2zz, dbc);
        for (int64 j = 0; j < nbc; ++j) {
            ARRAY2D(dbc, j, (ndim + i)) = (f2zz[j] - f1zz[j]) / (ep*2);
        }
    }

    if (ijac == 1) {
        free(u1zz);
        free(u2zz);
        free(f1zz);
        free(f2zz);
        return 0;
    }

    for (int64 i = 0; i < nfpr; ++i) {
        rtmp = HMACH;
        ep = rtmp*(fabs(par[icp[i]]) + 1);
        par[icp[i]] += ep;
        bcnd(ndim, par, icp, nbc, u0, u1, 0, f1zz, dbc);
        for (int64 j = 0; j < nbc; ++j) {
            ARRAY2D(dbc, j, (ndim*2) + icp[i]) = (f1zz[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }
    free(u1zz);
    free(u2zz);
    free(f1zz);
    free(f2zz);

    return 0;
}

int32
icni(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, double *u,
     double *uold, double *udot, double *upold, double *f, int64 ijac, double *dint) {
    int64 dint_dim1;

    double *u1zz;
    double *u2zz;

    int64 nfpr;
    double rtmp;
    double ep;
    int64 jac;
    int64 ijc;
    double umx, *f1zz, *f2zz;

    (void)rap;

    f1zz = xmalloc(sizeof(*(f1zz))*(usize)(iap->nint));
    f2zz = xmalloc(sizeof(*(f2zz))*(usize)(iap->nint));
    u1zz = xmalloc(sizeof(*(u1zz))*(usize)(iap->ndim));
    u2zz = xmalloc(sizeof(*(u2zz))*(usize)(iap->ndim));
    // Interface subroutine to user supplied ICND.

    // Local

    // Parameter adjustments

    dint_dim1 = nint;

    jac = iap->jac;
    nfpr = iap->nfpr;

    // Generate the integrand.

    if (jac == 0) {
        ijc = 0;
    } else {
        ijc = ijac;
    }
    icnd(ndim, par, icp, nint, u, uold, udot, upold, ijc, f, dint);

    if (jac == 1 || ijac == 0) {
        return 0;
    }

    // Generate the Jacobian by differencing.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            u1zz[j] = u[j];
            u2zz[j] = u[j];
        }
        u1zz[i] -= ep;
        u2zz[i] += ep;
        icnd(ndim, par, icp, nint, u1zz, uold, udot, upold, 0, f1zz, dint);
        icnd(ndim, par, icp, nint, u2zz, uold, udot, upold, 0, f2zz, dint);
        for (int64 j = 0; j < nint; ++j) {
            ARRAY2D(dint, j, i) = (f2zz[j] - f1zz[j]) / (ep*2);
        }
    }

    if (ijac == 1) {
        free(f1zz);
        free(f2zz);
        free(u1zz);
        free(u2zz);
        return 0;
    }

    for (int64 i = 0; i < nfpr; ++i) {
        rtmp = HMACH;
        ep = rtmp*(fabs(par[icp[i]]) + 1);
        par[icp[i]] += ep;
        icnd(ndim, par, icp, nint, u, uold, udot, upold, 0, f1zz, dint);
        for (int64 j = 0; j < nint; ++j) {
            ARRAY2D(dint, j, ndim + icp[i]) = (f1zz[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }
    free(f1zz);
    free(f2zz);
    free(u1zz);
    free(u2zz);

    return 0;
}

int32
fopi(iap_type *iap, rap_type *rap, int64 ndim, double *u, int64 *icp, double *par, int64 ijac,
     double *f, double *dfdu, double *dfdp) {
    double *u1zz;
    double *u2zz;
    int64 nfpr;

    double rtmp;
    double f1;
    double f2;
    double ep;
    int64 jac;
    int64 ijc;
    double umx;

    (void)rap;

    u1zz = xmalloc(sizeof(*(u1zz))*(usize)(iap->ndim));
    u2zz = xmalloc(sizeof(*(u2zz))*(usize)(iap->ndim));

    // Interface subroutine to user supplied FOPT.

    // Local

    // Parameter adjustments

    jac = iap->jac;
    nfpr = iap->nfpr;

    // Generate the objective function.

    if (jac == 0) {
        ijc = 0;
    } else {
        ijc = ijac;
    }
    fopt(ndim, u, icp, par, ijc, f, dfdu, dfdp);

    if (jac == 1 || ijac == 0) {
        free(u1zz);
        free(u2zz);
        return 0;
    }

    // Generate the Jacobian by differencing.

    umx = 0.;
    for (int64 i = 0; i < ndim; ++i) {
        if (fabs(u[i]) > umx) {
            umx = fabs(u[i]);
        }
    }

    rtmp = HMACH;
    ep = rtmp*(umx + 1);

    for (int64 i = 0; i < ndim; ++i) {
        for (int64 j = 0; j < ndim; ++j) {
            u1zz[j] = u[j];
            u2zz[j] = u[j];
        }
        u1zz[i] -= ep;
        u2zz[i] += ep;
        fopt(ndim, u1zz, icp, par, 0, &f1, dfdu, dfdp);
        fopt(ndim, u2zz, icp, par, 0, &f2, dfdu, dfdp);
        dfdu[i] = (f2 - f1) / (ep*2);
    }

    if (ijac == 1) {
        free(u1zz);
        free(u2zz);
        return 0;
    }

    for (int64 i = 0; i < nfpr; ++i) {
        rtmp = HMACH;
        ep = rtmp*(fabs(par[icp[i]]) + 1);
        par[icp[i]] += ep;
        fopt(ndim, u, icp, par, 0, &f1, dfdu, dfdp);
        dfdp[icp[i]] = (f1 - *f) / ep;
        par[icp[i]] -= ep;
    }

    free(u1zz);
    free(u2zz);
    return 0;
}
