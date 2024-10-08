/* autlib5.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

#include "auto_f2c.h"
#include "auto_c.h"
#include "x_auto.h"
#include "integers.h"
#include "autlib.h"
#include "xmalloc.h"

/* All of these global structures correspond to common
   blocks in the original code.  They are ONLY used within
   the Homcont code.
*/
static struct {
    int64 itwist;
    int64 istart;
    int64 iequib;
    int64 nfixed;
    int64 npsi;
    int64 nunstab;
    int64 nstab;
    int64 nrev;
} blhom_1;

static struct {
    int64 *ipsi, *ifixed, *irev;
} blhmp_1 = {NULL, NULL, NULL};

static struct {
    double *pu0;
    double *pu1;
} blhmu_1 = {NULL, NULL};

static struct {
    int64 nbcn;
} bcnn_1;

static struct {
    double compzero;
} blhma_1;

static struct {
    double *rr, *ri, *v, *vt, *xequib;
    int64 ineig;
} bleig_1 = {NULL, NULL, NULL, NULL, NULL, 0};

static struct {
    double *vrprev;
    int64 *ieigc;
} blhme_1 = {NULL, NULL};

static struct {
    double *cprev;
    int64 *iflag;
} beyn_1 = {NULL, NULL};

/* ----------------------------------------------------------------------- */
/*        Subroutines for Homoclinic Bifurcation Analysis */
/*       (A. R. Champneys, Yu. A. Kuznetsov, B. Sandstede) */
/* ----------------------------------------------------------------------- */

int32
fnho(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     int64 ijac, double *f, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 nfpr;
    double rtmp;
    double ep;
    int64 ndm;
    double umx;

    // Generates the equations for homoclinic bifurcation analysis

    // Local

    dfdp_dim1 = ndim;
    dfdu_dim1 = ndim;

    ndm = iap->ndm;
    nfpr = iap->nfpr;

    // Generate the function.

    ffho(iap, rap, ndim, u, uold, icp, par, f, ndm, global_scratch.dfu, global_scratch.dfp);

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
        for (int32 j = 0; j < ndim; ++j) {
            global_scratch.uu1[j] = u[j];
            global_scratch.uu2[j] = u[j];
        }
        global_scratch.uu1[i] -= ep;
        global_scratch.uu2[i] += ep;
        ffho(iap, rap, ndim, global_scratch.uu1, uold, icp, par, global_scratch.ff1, ndm,
             global_scratch.dfu, global_scratch.dfp);
        ffho(iap, rap, ndim, global_scratch.uu2, uold, icp, par, global_scratch.ff2, ndm,
             global_scratch.dfu, global_scratch.dfp);
        for (int32 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdu, j, i) = (global_scratch.ff2[j] - global_scratch.ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        ffho(iap, rap, ndim, u, uold, icp, par, global_scratch.ff1, ndm, global_scratch.dfu,
             global_scratch.dfp);
        for (int32 j = 0; j < ndim; ++j) {
            ARRAY2D(dfdp, j, icp[i]) = (global_scratch.ff1[j] - f[j]) / ep;
        }
        par[icp[i]] -= ep;
    }

    return 0;
}

int32
ffho(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, int64 *icp, double *par,
     double *f, int64 ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    double dum1;

    dfdp_dim1 = ndm;
    dfdu_dim1 = ndm;

    ndm = iap->ndm;

    if (blhom_1.itwist == 0) {
        //        *Evaluate the R.-H. sides
        funi(iap, rap, ndm, u, uold, icp, par, 0, f, dfdu, dfdp);
    } else {
        //        *Adjoint variational equations for normal vector
        funi(iap, rap, ndm, u, uold, icp, par, 1, f, dfdu, dfdp);
        //        *Set F = - (Df)^T u
        for (int32 j = 0; j < ndm; ++j) {
            dum1 = 0.;
            for (int64 i = 0; i < ndm; ++i) {
                dum1 += ARRAY2D(dfdu, i, j)*u[ndm + i];
            }
            f[ndm + j] = -dum1;
        }
        //        *Set F =  F + PAR(10)*f
        for (int32 j = 0; j < ndm; ++j) {
            f[ndm + j] += par[9]*f[j];
        }
    }

    // Scale by truncation interval T=PAR(11)

    for (int64 i = 0; i < ndim; ++i) {
        f[i] = par[10]*f[i];
    }

    return 0;
}

int32
bcho(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nbc, double *u0,
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

    /* Generates the boundary conditions for homoclinic bifurcation analysis
     */

    // Local

    dbc_dim1 = nbc;

    nbc0 = iap->nbc0;
    nfpr = iap->nfpr;

    // Generate the function.

    fbho(iap, rap, ndim, par, icp, nbc, nbc0, &u0[0], &u1[0], &f[0], dfu);

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
        for (int32 j = 0; j < ndim; ++j) {
            uu1[j] = u0[j];
            uu2[j] = u0[j];
        }
        uu1[i] -= ep;
        uu2[i] += ep;
        fbho(iap, rap, ndim, par, icp, nbc, nbc0, uu1, u1, ff1, dfu);
        fbho(iap, rap, ndim, par, icp, nbc, nbc0, uu2, u1, ff2, dfu);
        for (int32 j = 0; j < nbc; ++j) {
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
        for (int32 j = 0; j < ndim; ++j) {
            uu1[j] = u1[j];
            uu2[j] = u1[j];
        }
        uu1[i] -= ep;
        uu2[i] += ep;
        fbho(iap, rap, ndim, par, icp, nbc, nbc0, u0, uu1, ff1, dfu);
        fbho(iap, rap, ndim, par, icp, nbc, nbc0, u0, uu2, ff2, dfu);
        for (int32 j = 0; j < nbc; ++j) {
            ARRAY2D(dbc, j, (ndim + i)) = (ff2[j] - ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        fbho(iap, rap, ndim, par, icp, nbc, nbc0, u0, u1, ff2, dfu);
        for (int32 j = 0; j < nbc; ++j) {
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
fbho(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nbc, int64 nbc0,
     double *u0, double *u1, double *fb, double *dbc) {
    int64 dbc_dim1;

    int64 ieig;

    int64 k;

    int64 ineig;

    int64 jb;
    int64 ip;
    int64 kp;

    int64 ijc = 0;
    int64 ndm;
    double dum;
    double dum1;
    double dum2;

    double *f;
    double *bound2;
    double *fj;
    double *ri;
    double *rr, *vr, *vt;
    double *xequib1;
    double *xequib2;

    (void)nbc0;
    (void)rap;

    /* I am not 100% sure if this is supposed to be iap->ndm or iap->ndim,
       but it appears from looking at the code that it should be iap->ndm.
       Also, note that I have replaced the occurances of N X in the algorithms
       with (iap->ndm), so if you change it to iap->ndim here you will
       need to make the similiar changes in the algorithms.
       Finally, the routines called from here prjcti and eighi
       also depend on these arrays, and more importantly the algorithm,
       having N X.  So, they all need to be changed at once.
    */
    f = xmalloc(sizeof(*f)*(usize)(iap->ndm));
    bound2 = xmalloc(sizeof(*bound2)*(usize)((iap->ndm)*(iap->ndm)));
    fj = xmalloc(sizeof(*fj)*(usize)(iap->ndm));
    ri = xmalloc(sizeof(*ri)*(usize)(iap->ndm));
    rr = xmalloc(sizeof(*rr)*(usize)(iap->ndm));
    vr = xmalloc(sizeof(*vr)*(usize)((iap->ndm)*(iap->ndm)));
    vt = xmalloc(sizeof(*vt)*(usize)((iap->ndm)*(iap->ndm)));
    xequib1 = xmalloc(sizeof(*(xequib1))*(usize)(iap->ndm));
    xequib2 = xmalloc(sizeof(*(xequib2))*(usize)(iap->ndm));

    // Generates the boundary conditions for homoclinic orbits.

    // Local

    dbc_dim1 = nbc;

    ndm = iap->ndm;

    //     *Initialization
    for (int64 i = 1; i <= nbc; ++i) {
        fb[-1 + i] = 0.;
    }
    jb = 1;

    //     *Update pu0,pu1
    for (int64 i = 0; i < ndim; ++i) {
        blhmu_1.pu0[i] = u0[i];
        blhmu_1.pu1[i] = u1[i];
    }

    if (blhom_1.iequib == 0 || blhom_1.iequib == -1) {
        pvls(ndm, u0, par);
    }
    //              write(9,*) 'Xequib:'
    for (int64 i = 0; i < ndm; ++i) {
        xequib1[i] = par[i + 11];
        //              write(9,*) I,XEQUIB1(I)
    }
    if (blhom_1.iequib >= 0) {
        for (int64 i = 0; i < ndm; ++i) {
            xequib2[i] = par[i + 11];
        }
    } else {
        for (int64 i = 0; i < ndm; ++i) {
            xequib2[i] = par[ndm + 11 + i];
        }
    }

    //     **Regular Continuation**
    if (blhom_1.istart != 3) {
        //        *Projection boundary conditions for the homoclinic orbit
        //        *NSTAB boundary conditions at t=0
        prjcti(bound2, xequib1, icp, par, -1, 1, 1, &ndm);
        for (int64 i = 0; i < blhom_1.nstab; ++i) {
            for (k = 0; k < ndm; ++k) {
                fb[-1 + jb] += (u0[k] - xequib1[k])*bound2[i + k*(iap->ndm)];
            }
            //         write(9,*) 'fb',jb,fb(jb)
            ++jb;
        }
        //        *NUNSTAB boundary conditions at t=1
        if (blhom_1.nrev == 0) {
            prjcti(bound2, xequib2, icp, par, 1, 2, 1, &ndm);
            for (int64 i = ndm - blhom_1.nunstab; i < ndm; ++i) {
                for (k = 0; k < ndm; ++k) {
                    fb[-1 + jb] += (u1[k] - xequib2[k])*bound2[i + k*(iap->ndm)];
                }
                ++jb;
            }
        } else {
            /*         *NUNSTAB symmetric boundary conditions at t=1 if NREV=1
             */

            for (int64 i = 0; i < ndim; ++i) {
                if (blhmp_1.irev[i] > 0) {
                    if (blhmp_1.irev[i] == 1) {
                        // *****NOTE MODIFICATION FROM GENERAL CASE
                        fb[-1 + jb] = sin(u1[i]);
                        //                        FB(JB)=U1(I)
                    } else {
                        fb[-1 + jb] = u1[i];
                    }
                    ++jb;
                }
            }
        }
        ieig = 0;
        ineig = 0;
        /*        *NFIXED extra boundary conditions for the fixed conditions
         */
        if (blhom_1.nfixed > 0) {
            if (ieig == 0) {
                eighi(1, 2, rr, ri, vr, xequib1, icp, par, &ndm);
                ieig = 1;
            }
            for (int64 i = 0; i < blhom_1.nfixed; ++i) {
                if (blhmp_1.ifixed[i] > 10 && ineig == 0) {
                    eighi(1, 1, rr, ri, vt, xequib1, icp, par, &ndm);
                    ineig = 1;
                }
                fb[-1 + jb] = psiho(iap, blhmp_1.ifixed[i], rr, ri, vr, vt, icp, par);
                ++jb;
            }
        }
        /*        *NDM initial conditions for the equilibrium if IEQUIB=1,2,-2
         */
        if (blhom_1.iequib != 0 && blhom_1.iequib != -1) {
            func(ndm, xequib1, icp, par, 0, f, &dum1, &dum2);
            for (int64 i = 0; i < ndm; ++i) {
                fb[-1 + jb] = f[i];
                ++jb;
            }
            /*        *NDM extra initial conditions for the equilibrium if IEQ
                      UIB=-2 */
            if (blhom_1.iequib == -2) {
                func(ndm, xequib2, icp, par, 0, f, &dum1, &dum2);
                for (int64 i = 0; i < ndm; ++i) {
                    fb[-1 + jb] = f[i];
                    ++jb;
                }
            }
        }
        /*       *extra boundary condition in the case of a saddle-node homocl
                 inic*/
        if (blhom_1.iequib == 2) {
            if (ineig == 0) {
                eighi(1, 1, rr, ri, vt, xequib1, icp, par, &ndm);
                ineig = 1;
            }
            fb[-1 + jb] = rr[blhom_1.nstab];
            ++jb;
        }
        //        *boundary conditions for normal vector
        if (blhom_1.itwist == 1) {
            /*           *-orthogonal to the unstable directions of A  at t=0
             */
            prjcti(bound2, xequib1, icp, par, 1, 1, 2, &ndm);
            for (int64 i = ndm - blhom_1.nunstab; i < ndm; ++i) {
                dum = 0.;
                for (k = 0; k < ndm; ++k) {
                    dum += u0[ndm + k]*bound2[i + k*(iap->ndm)];
                }
                fb[-1 + jb] = dum;
                ++jb;
            }
            //           *-orthogonal to the stable directions of A  at t=1
            prjcti(bound2, xequib2, icp, par, -1, 2, 2, &ndm);
            for (int64 i = 0; i < blhom_1.nstab; ++i) {
                dum = 0.;
                for (k = 0; k < ndm; ++k) {
                    dum += u1[ndm + k]*bound2[i + k*(iap->ndm)];
                }
                fb[-1 + jb] = dum;
                ++jb;
            }
            free(f);
            free(bound2);
            free(fj);
            free(ri);
            free(rr);
            free(vr);
            free(vt);
            free(xequib1);
            free(xequib2);
            return 0;
        }
    } else {
        //     **Starting Solutions using Homotopy**
        jb = 0;
        for (int64 i = 1; i <= nbc; ++i) {
            fb[-1 + i] = 0.;
        }
        ineig = 0;
        ip = 12;
        if (blhom_1.iequib >= 0) {
            ip += ndm;
        } else {
            ip += ndm << 1;
        }
        kp = ip;
        //        *Explicit boundary conditions for homoclinic orbit at t=0
        eighi(1, 2, rr, ri, vr, xequib1, icp, par, &ndm);
        ieig = 1;
        if (blhom_1.nunstab > 1) {
            dum = 0.;
            kp = ip + blhom_1.nunstab;
            jb = ndm + 1;
            for (int32 j = 0; j < blhom_1.nunstab; ++j) {
                for (int64 i = 0; i < ndm; ++i) {
                    fb[i] =
                        u0[i] - xequib1[i] - par[ip + j]*vr[blhom_1.nstab + j + i*(iap->ndm)];
                }
                // Computing 2nd power
                dum += par[ip + j]*par[ip + j];
            }
            jb = ndm + 1;
            fb[-1 + jb] = dum - par[-1 + ip];
            ++jb;
        } else {
            kp = ip + 1;
            jb = ndm;
            for (int64 i = 0; i < ndm; ++i) {
                fb[i] = u0[i] - xequib1[i] -
                        par[-1 + ip]*par[ip]*vr[blhom_1.nstab + i*(iap->ndm)];
            }
            jb = ndm + 1;
        }
        /*        *Projection boundary conditions for the homoclinic orbit at
    t=1 */
        if (ineig == 0) {
            eighi(1, 1, rr, ri, vt, xequib2, icp, par, &ndm);
            ineig = 1;
        }
        for (int64 i = 0; i < blhom_1.nunstab; ++i) {
            k = i + blhom_1.nstab;
            dum = 0.;
            for (int32 j = 0; j < ndm; ++j) {
                dum += (u1[j] - xequib2[j])*vt[k + j*(iap->ndm)];
            }
            ++kp;
            fb[-1 + jb] = dum - par[-1 + kp];
            ++jb;
        }
        /*        *NDM initial conditions for the equilibrium if IEQUIB=1,2,-2
         */
        if (blhom_1.iequib != 0 && blhom_1.iequib != -1) {
            func(ndm, xequib1, icp, par, 0, f, &dum1, &dum2);
            for (int64 i = 0; i < ndm; ++i) {
                fb[-1 + jb] = f[i];
                ++jb;
            }
            /*        *NDM extra initial conditions for the equilibrium if IEQ
                      UIB=-2 */
            if (blhom_1.iequib == -2) {
                func(ndm, xequib2, icp, par, 0, f, &dum1, &dum2);
                for (int64 i = 0; i < ndm; ++i) {
                    fb[-1 + jb] = f[i];
                    ++jb;
                }
            }
        }
    }

    //      write(9,*) NBCN,NBC
    // *user defined extra boundary conditions
    if (bcnn_1.nbcn > 0) {
        bcnd(ndim, par, icp, bcnn_1.nbcn, u0, u1, ijc, fj, dbc);
        for (k = 0; k < bcnn_1.nbcn; ++k) {
            fb[-1 + jb] = fj[k];
            //            write(9,*),fb(jb),par(30)
            ++jb;
        }
    }

    free(f);
    free(bound2);
    free(fj);
    free(ri);
    free(rr);
    free(vr);
    free(vt);
    free(xequib1);
    free(xequib2);

    return 0;
}

int32
icho(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, double *u,
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

    // Generates integral conditions for homoclinic bifurcation analysis

    // Local

    dint_dim1 = nint;

    nnt0 = iap->nnt0;
    nfpr = iap->nfpr;

    // Generate the function.

    fiho(iap, rap, ndim, par, icp, nint, nnt0, u, uold, udot, upold, f, dfu);

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
        for (int32 j = 0; j < ndim; ++j) {
            uu1[j] = u[j];
            uu2[j] = u[j];
        }
        uu1[i] -= ep;
        uu2[i] += ep;
        fiho(iap, rap, ndim, par, icp, nint, nnt0, uu1, uold, udot, upold, ff1, dfu);
        fiho(iap, rap, ndim, par, icp, nint, nnt0, uu2, uold, udot, upold, ff2, dfu);
        for (int32 j = 0; j < nint; ++j) {
            ARRAY2D(dint, j, i) = (ff2[j] - ff1[j]) / (ep*2);
        }
    }

    for (int64 i = 0; i < nfpr; ++i) {
        par[icp[i]] += ep;
        fiho(iap, rap, ndim, par, icp, nint, nnt0, u, uold, udot, upold, ff1, dfu);
        for (int32 j = 0; j < nint; ++j) {
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
fiho(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, int64 nint, int64 nnt0,
     double *u, double *uold, double *udot, double *upold, double *fi, double *dint) {
    int64 dint_dim1;

    int64 ijac = 0;

    int64 jb;
    double *fj;
    int64 ndm;
    double dum;

    (void)nnt0;
    (void)ndim;
    (void)rap;

    fj = xmalloc(sizeof(*fj)*(usize)(iap->ndim));
    // Generates the integral conditions for homoclinic orbits.

    dint_dim1 = nint;

    ndm = iap->ndm;
    jb = 0;

    // Integral phase condition for homoclinic orbit

    if (blhom_1.nrev == 0) {
        dum = 0.;
        for (int64 i = 0; i < ndm; ++i) {
            dum += upold[i]*(u[i] - uold[i]);
        }
        ++jb;
        fi[-1 + jb] = dum;
    }

    // Integral phase condition for adjoint equation

    if (blhom_1.itwist == 1) {
        dum = 0.;
        for (int64 i = 0; i < ndm; ++i) {
            dum += uold[ndm + i]*(u[ndm + i] - uold[ndm + i]);
        }
        ++jb;
        fi[1] = dum;
    }

    // User-defined integral constraints

    if (jb < nint) {
        icnd(ndm, par, icp, nint, u, uold, udot, upold, ijac, fj, dint);
        for (int64 i = 0; i < nint - jb; ++i) {
            fi[i + jb] = fj[i];
        }
    }
    free(fj);
    return 0;
}

int32
inho(iap_type *iap, int64 *icp, double *par) {
    int64 ndim;
    int64 nint;
    int64 nuzr;
    int64 nfree;
    int64 icorr;
    int64 nbc;
    int64 ndm;
    int64 irs;
    int64 isw;

    // Allocate memory for global structures.
    free(blhmp_1.ipsi);
    free(blhmp_1.ifixed);
    free(blhmp_1.irev);

    free(blhmu_1.pu0);
    free(blhmu_1.pu1);

    free(bleig_1.rr);
    free(bleig_1.ri);
    free(bleig_1.v);
    free(bleig_1.vt);
    free(bleig_1.xequib);

    free(blhme_1.vrprev);
    free(blhme_1.ieigc);

    free(beyn_1.cprev);
    free(beyn_1.iflag);

    blhmp_1.ipsi = xmalloc(sizeof(*(blhmp_1.ipsi))*NPARX);
    blhmp_1.ifixed = xmalloc(sizeof(*(blhmp_1.ifixed))*NPARX);
    blhmp_1.irev = xmalloc(sizeof(*(blhmp_1.irev))*(usize)(iap->ndim));

    blhme_1.ieigc = xmalloc(sizeof(*(blhme_1.ieigc))*2);

    beyn_1.iflag = xmalloc(sizeof(*(beyn_1.iflag))*4);

    /* the prjctn_ function uses this array to test if this is
         the first time the prjctn_ function has been called.
         Accordingly, I initialize it to zero here after I
         have created it. */
    for (int64 i = 0; i < 4; i++) {
        beyn_1.iflag[i] = 0;
    }

    // Reads from fort.11 specific constants for homoclinic continuation.
    // Sets up re-defined constants in IAP.
    // Sets other constants in the following common blocks.

    // set various constants

    ndim = iap->ndim;
    irs = iap->irs;
    isw = iap->isw;
    nbc = iap->nbc;
    nint = iap->nint;
    nuzr = iap->nuzr;
    ndm = ndim;
    blhma_1.compzero = HMACHHO;
    blhom_1.nunstab = x_auto.nunstab;
    blhom_1.nstab = x_auto.nstab;
    blhom_1.iequib = x_auto.iequib;
    blhom_1.itwist = 0;
    blhom_1.istart = 2;
    blhom_1.nrev = 0;
    blhom_1.nfixed = 0;
    blhom_1.npsi = 0;
    /* printf("%d %d %d %d
     * %d\n",blhom_1.nunstab,blhom_1.nstab,blhom_1.iequib,blhom_1.itwist,blhom_1.istart);
     */
    // updated reading in of constants for reversible equations
    // replaces location in datafile of compzero

    ndim = ndm*(blhom_1.itwist + 1);
    /* Allocate memory for global structures.  We didn't know the
       size for these until ndim was computed. */

    blhmu_1.pu0 = xmalloc(sizeof(*(blhmu_1.pu0))*(usize)ndim);
    blhmu_1.pu1 = xmalloc(sizeof(*(blhmu_1.pu1))*(usize)ndim);

    bleig_1.rr = xmalloc(sizeof(*(bleig_1.rr))*(usize)ndim);
    bleig_1.ri = xmalloc(sizeof(*(bleig_1.ri))*(usize)ndim);
    bleig_1.v = xmalloc(sizeof(*(bleig_1.v))*(usize)(ndim*ndim));
    bleig_1.vt = xmalloc(sizeof(*(bleig_1.vt))*(usize)(ndim*ndim));
    bleig_1.xequib = xmalloc(sizeof(*(bleig_1.xequib))*(usize)ndim);

    blhme_1.vrprev = xmalloc(sizeof(*(blhme_1.vrprev))*(usize)(2*ndim*ndim));

    beyn_1.cprev = xmalloc(sizeof(*(beyn_1.cprev))*(usize)(2*2*ndim*ndim));
    nfree = blhom_1.nfixed + 2 - blhom_1.nrev + nint + nbc;
    bcnn_1.nbcn = nbc;

    // Free parameter (artificial parameter for psi)
    // nondegeneracy parameter of the adjoint

    if (blhom_1.itwist == 1) {
        ++nfree;
        icp[-1 + nfree] = 9;
        par[9] = 0.;
    }

    // Extra free parameters for equilibrium if iequib=1,2,-2

    if (blhom_1.iequib != 0 && blhom_1.iequib != -1) {
        for (int64 i = 0; i < ndm; ++i) {
            icp[nfree + i] = i + 11;
        }
    }

    if (blhom_1.iequib == -2) {
        for (int64 i = 0; i < ndm; ++i) {
            icp[nfree + ndm + i] = ndm + 11 + i;
        }
    }

    if (blhom_1.istart != 3) {
        //     *regular continuation

        nint = nint + blhom_1.itwist + 1 - blhom_1.nrev;

        if (isw == 2) {
            icorr = 2;
        } else {
            icorr = 1;
        }
        nbc = blhom_1.nstab + blhom_1.nunstab + (blhom_1.itwist + blhom_1.iequib)*ndm + nfree -
              nint - icorr;
        if (blhom_1.iequib == 2) {
            nbc = nbc - ndm + 1;
        }
        if (blhom_1.iequib < 0) {
            nbc -= (blhom_1.iequib*3 + 2)*ndm;
        }
    } else {
        //     *starting solutions using homotopy
        if (blhom_1.nunstab == 1) {
            nbc = ndm*(blhom_1.iequib + 1) + 1;
        } else {
            nbc = ndm*(blhom_1.iequib + 1) + blhom_1.nunstab + 1;
        }
        if (blhom_1.iequib == 2) {
            fprintf(fp9, "WARNING: IEQUIB=2 NOT ALLOWED WITH ISTART=3\n");
        }
        if (blhom_1.iequib < 0) {
            nbc -= ndm*(blhom_1.iequib*3 + 2);
        }
        nint = 0;
    }

    // write new constants into IAP

    iap->ndim = ndim;
    iap->nbc = nbc;
    iap->nint = nint;
    iap->nuzr = nuzr;
    iap->ndm = ndm;

    return 0;
}

int32
preho(int64 *ndx, int64 *ntsr, int64 *nar, int64 *ndim, int64 *ncolrs, double *ups, double *udotps,
      double *tm, double *par) {
    int64 ups_dim1;
    int64 udotps_dim1;

    int64 jmin;
    double upsi;
    int64 i;
    int64 j;
    double tmmin;
    int64 k1;
    int64 k2;
    int64 ii;
    double upsmin;
    int64 ist;

    // Preprocesses (perturbs) restart data to enable
    // initial computation of the adjoint variable

    udotps_dim1 = *ndx;
    ups_dim1 = *ndx;

    if (*nar < *ndim) {
        for (j = 0; j < *ntsr; ++j) {
            for (i = 0; i < *ncolrs; ++i) {
                k1 = i**ndim;
                k2 = (i + 1)**ndim - 1;
                for (int64 k = k1 + *nar; k <= k2; ++k) {
                    ARRAY2D(ups, j, k) = .1;
                }
            }
        }
        for (int64 k = *nar; k < *ndim; ++k) {
            ARRAY2D(ups, *ntsr, k) = .1;
        }
    }

    // Shift phase if necessary if continu(e)ing from
    // a periodic orbit into a homoclinic one

    if (blhom_1.istart == 1) {
        // First find smallest value in norm

        upsmin = 1e20;
        jmin = 1;
        for (j = 0; j < *ntsr + 1; ++j) {
            upsi = 0.;
            for (i = 0; i < *nar; ++i) {
                upsi += (ARRAY2D(ups, j, i) - par[i + 11])*(ARRAY2D(ups, j, i) - par[i + 11]);
            }
            if (upsi <= upsmin) {
                upsmin = upsi;
                jmin = j + 1;
            }
        }
        tmmin = tm[-1 + jmin];

        // And then do the actual shift

        if (jmin != 1) {
            ist = 0;
            j = *ntsr + 1;
            for (ii = 0; ii < *ntsr; ++ii) {
                if (j == *ntsr + 1) {
                    ++ist;
                    tm[-1 + j] = tm[-1 + ist];
                    for (int64 k = 0; k < *ncolrs**ndim; ++k) {
                        ARRAY2D(ups, (j - 1), k) = ARRAY2D(ups, (ist - 1), k);
                        ARRAY2D(udotps, (j - 1), k) = ARRAY2D(udotps, (ist - 1), k);
                    }
                    j = ist;
                }
                i = j;
                j = j + jmin - 1;
                if (j > *ntsr) {
                    j -= *ntsr;
                }
                if (j == ist) {
                    j = *ntsr + 1;
                }
                tm[-1 + i] = tm[-1 + j] - tmmin;
                if (tm[-1 + i] < 0.) {
                    tm[-1 + i] += 1.;
                }
                for (int64 k = 0; k < *ncolrs**ndim; ++k) {
                    ARRAY2D(ups, (i - 1), k) = ARRAY2D(ups, (j - 1), k);
                    ARRAY2D(udotps, (i - 1), k) = ARRAY2D(udotps, (j - 1), k);
                }
            }

            // Last equal to first

            tm[*ntsr] = 1.;
            for (int64 k = 0; k < *ncolrs**ndim; ++k) {
                ARRAY2D(ups, *ntsr, k) = ARRAY2D(ups, 0, k);
                ARRAY2D(udotps, *ntsr, k) = ARRAY2D(udotps, 0, k);
            }
        }
    }

    return 0;
}

int32
stpnho(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr, int64 *ncolrs,
       double *rlcur, double *rldot, int64 *ndxloc, double *ups, double *udotps, double *upoldp,
       double *tm, double *dtm, int64 *nodir, double *thl, double *thu) {
    int64 ups_dim1;
    int64 udotps_dim1;

    int64 ndim;
    int64 ncol;
    int64 nfpr;
    int64 ntst;
    int64 ncol1;
    double t;
    double *u;
    int64 k1;
    int64 k2;

    double dt;
    int64 lab;
    int64 ibr;

    (void)thu;
    (void)thl;
    (void)dtm;
    (void)upoldp;
    (void)udotps;
    (void)rldot;
    (void)rap;

    u = xmalloc(sizeof(*u)*(usize)(iap->ndim));
    // Generates a starting point for the continuation of a branch of
    // of solutions to general boundary value problems by calling the user
    // supplied subroutine STPNT where an analytical solution is given.

    // Local

    udotps_dim1 = *ndxloc;
    ups_dim1 = *ndxloc;

    ndim = iap->ndim;
    ntst = iap->ntst;
    ncol = iap->ncol;
    nfpr = iap->nfpr;

    // Generate the (initially uniform) mesh.

    msh(iap, tm);
    dt = 1. / (double)(ntst*ncol);

    for (int32 j = 0; j < ntst + 1; ++j) {
        if (j == ntst) {
            ncol1 = 1;
        } else {
            ncol1 = ncol;
        }
        for (int64 i = 0; i < ncol1; ++i) {
            t = tm[j] + (double)i*dt;
            k1 = i*ndim;
            k2 = (i + 1)*ndim - 1;
            stpho(iap, icp, u, par, &t);
            for (int64 k = k1; k <= k2; ++k) {
                ARRAY2D(ups, j, k) = u[k - k1];
            }
        }
    }

    *ntsr = ntst;
    *ncolrs = ncol;
    ibr = 1;
    iap->ibr = ibr;
    lab = 0;
    iap->lab = lab;

    for (int64 i = 0; i < nfpr; ++i) {
        rlcur[i] = par[icp[i]];
    }

    *nodir = 1;
    free(u);
    return 0;
}

int32
stpho(iap_type *iap, int64 *icp, double *u, double *par, double *t) {
    int64 ip;
    int64 kp;
    int64 ndm;

    double *ri;
    double *rr, *vr, *vt, *xequib;

    ri = xmalloc(sizeof(*ri)*(usize)(iap->ndm));
    rr = xmalloc(sizeof(*rr)*(usize)(iap->ndm));
    vr = xmalloc(sizeof(*vr)*(usize)((iap->ndm)*(iap->ndm)));
    vt = xmalloc(sizeof(*vt)*(usize)((iap->ndm)*(iap->ndm)));
    xequib = xmalloc(sizeof(*xequib)*(usize)(iap->ndm));

    // Generates a starting point for homoclinic continuation
    // If ISTART=2 it calls STPNHO.
    // If ISTART=3 it sets up the homotopy method.

    // Local

    ndm = iap->ndm;

    // Initialize parameters

    stpnt(ndm, *t, u, par);

    // Initialize solution and additional parameters

    switch ((int32)blhom_1.istart) {
    case 1:
        goto L1;
    case 2:
        goto L2;
    case 3:
        goto L3;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }

L1:
    // Obsolete option

    free(ri);
    free(rr);
    free(vr);
    free(vt);
    free(xequib);
    return 0;

L2:
    //     *Regular continuation (explicit solution in STHO)

    free(ri);
    free(rr);
    free(vr);
    free(vt);
    free(xequib);
    return 0;

L3:
    //     *Starting solutions using homotopy

    pvls(ndm, u, par);
    for (int64 i = 0; i < ndm; ++i) {
        xequib[i] = par[i + 11];
    }
    eighi(1, 1, rr, ri, vt, xequib, icp, par, &ndm);
    eighi(1, 2, rr, ri, vr, xequib, icp, par, &ndm);

    // Set up artificial parameters at the left-hand end point of orbit

    ip = 12;
    if (blhom_1.iequib >= 0) {
        ip += ndm;
    } else {
        ip += ndm*2;
    }
    kp = ip;

    // Parameters xi 1=1, xi i=0, i=2,NSTAB

    par[ip] = 1.;
    if (blhom_1.nunstab > 1) {
        for (int64 i = 1; i < blhom_1.nunstab; ++i) {
            par[ip + i] = 0.;
        }
    }
    ip += blhom_1.nunstab;

    /*Starting guess for homoclinic orbit in real principal unstable direction
     */

    for (int64 i = 0; i < ndm; ++i) {
        u[i] = xequib[i] + vr[blhom_1.nstab + i*(iap->ndm)]*par[-1 + kp]*par[kp] *
                               exp(rr[blhom_1.nstab]**t*par[10]);
    }
    for (int64 i = 0; i < ndm; ++i) {
        fprintf(fp9, "stpho %20.10f\n", u[i]);
    }
    fprintf(fp9, "\n");

    // Artificial parameters at the right-hand end point of the orbit
    // omega_i=<x(1)-x_o,w_i^*>

    for (int64 i = 0; i < blhom_1.nunstab; ++i) {
        par[ip + i] = 0.;
        for (int32 j = 0; j < ndm; ++j) {
            par[ip + i] += vr[blhom_1.nstab + j*(iap->ndm)]*par[-1 + kp]*par[kp] *
                           exp(rr[blhom_1.nstab]*par[10]) *
                           vt[blhom_1.nstab + i + j*(iap->ndm)];
        }
    }
    ip += blhom_1.nunstab;
    free(ri);
    free(rr);
    free(vr);
    free(vt);
    free(xequib);
    return 0;
}

int32
pvlsho(iap_type *iap, rap_type *rap, int64 *icp, double *dtm, int64 *ndxloc, double *ups,
       int64 *ndim, double *p0, double *p1, double *par) {
    int64 ups_dim1;
    int64 p0_dim1;
    int64 p1_dim1;

    double orient;

    int64 iid;
    int64 ndm;

    ups_dim1 = *ndxloc;
    p1_dim1 = *ndim;
    p0_dim1 = *ndim;

    iid = iap->iid;
    ndm = iap->ndm;

    pvlsbv(iap, rap, icp, dtm, ndxloc, ups, ndim, p0, p1, par);

    //      *Compute eigenvalues
    bleig_1.ineig = 0;
    for (int64 i = 0; i < ndm; ++i) {
        bleig_1.xequib[i] = par[i + 11];
    }
    eighi(1, 2, bleig_1.rr, bleig_1.ri, bleig_1.v, bleig_1.xequib, icp, par, &ndm);
    if (iid >= 3) {
        fprintf(fp9, "EIGENVALUES\n");
        for (int32 j = 0; j < ndm; ++j) {
            fprintf(fp9, " (%12.7f %12.7f)\n", bleig_1.rr[j], bleig_1.ri[j]);
        }
    }
    if (blhom_1.itwist == 1) {
        eighi(1, 1, bleig_1.rr, bleig_1.ri, bleig_1.vt, bleig_1.xequib, icp, par, &ndm);
        bleig_1.ineig = 1;
        orient = psiho(iap, 0, bleig_1.rr, bleig_1.ri, bleig_1.v, bleig_1.vt, icp, par);
        if (iid >= 3) {
            if (orient < 0.) {
                fprintf(fp9, " Non-orientable, (%20.10f)\n", orient);
            } else {
                fprintf(fp9, " Orientable (%20.10f)\n", orient);
            }
        }
    }

    for (int64 i = 0; i < blhom_1.npsi; ++i) {
        if (blhmp_1.ipsi[i] > 10 && bleig_1.ineig == 0) {
            eighi(1, 1, bleig_1.rr, bleig_1.ri, bleig_1.vt, bleig_1.xequib, icp, par, &ndm);
            bleig_1.ineig = 1;
        }
        par[blhmp_1.ipsi[i] + 19] =
            psiho(iap, blhmp_1.ipsi[i], bleig_1.rr, bleig_1.ri, bleig_1.v, bleig_1.vt, icp, par);
        if (iid >= 3) {
            fprintf(fp9, " PSI(%2ld)=%20.10f\n", blhmp_1.ipsi[i], par[blhmp_1.ipsi[i] + 19]);
        }
    }

    return 0;
}

double
psiho(iap_type *iap, int64 is, double *rr, double *ri, double *v, double *vt, int64 *icp,
      double *par) {
    double ret_val;

    double *f0, *f1, droot, s1, s2, f0norm, f1norm, u0norm, u1norm;
    int64 ndm;
    double dum1;
    double dum2;

    f0 = xmalloc(sizeof(*(f0))*(usize)(iap->ndm));
    f1 = xmalloc(sizeof(*(f1))*(usize)(iap->ndm));

    /* The conditions for degenerate homoclinic orbits are given by PSI(IS)=0.
     */

    /* RR and RI contain the real and imaginary parts of eigenvalues which are
     */
    // ordered with respect to their real parts (smallest first).
    // The (generalised) real eigenvectors are stored as the ROWS of V.
    // The (generalised) real left eigenvectors are in the ROWS of VT.
    // In the block ENDPTS are stored the co-ordinates of the left (PU0)
    /*and right (PU1) endpoints of the solution (+  vector if that is computed
    )*/

    // Local

    vt -= ((iap->ndm) + 1);
    v -= ((iap->ndm) + 1);

    ndm = iap->ndm;

    func(ndm, blhmu_1.pu0, icp, par, 0, f0, &dum1, &dum2);
    func(ndm, blhmu_1.pu1, icp, par, 0, f1, &dum1, &dum2);

    ret_val = 0.;

    //  Compute orientation

    if (is == 0) {
        s1 = 0.;
        s2 = 0.;
        f0norm = 0.;
        f1norm = 0.;
        u0norm = 0.;
        u1norm = 0.;
        for (int32 j = 0; j < ndm; ++j) {
            s1 += f1[j]*blhmu_1.pu0[ndm + j];
            s2 += f0[j]*blhmu_1.pu1[ndm + j];
            // Computing 2nd power
            f0norm += f0[j]*f0[j];
            // Computing 2nd power
            f1norm += f1[j]*f1[j];
            // Computing 2nd power
            u0norm += blhmu_1.pu0[j + ndm]*blhmu_1.pu0[j + ndm];
            // Computing 2nd power
            u1norm += blhmu_1.pu1[j + ndm]*blhmu_1.pu1[j + ndm];
        }
        droot = sqrt(f0norm*f1norm*u0norm*u1norm);
        if (droot != 0.) {
            ret_val = -s1*s2 / droot;
        } else {
            ret_val = 0.;
        }
        free(f0);
        free(f1);
        return ret_val;
    }

    switch ((int32)is) {
    case 1:
        goto L1;
    case 2:
        goto L2;
    case 3:
        goto L3;
    case 4:
        goto L4;
    case 5:
        goto L5;
    case 6:
        goto L6;
    case 7:
        goto L7;
    case 8:
        goto L8;
    case 9:
        goto L9;
    case 10:
        goto L10;
    case 11:
        goto L11;
    case 12:
        goto L12;
    case 13:
        goto L13;
    case 14:
        goto L14;
    case 15:
        goto L15;
    case 16:
        goto L16;
    default:
        fprintf(stderr, "Unexpected case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }

    // Resonant eigenvalues (neutral saddle)

L1:
    ret_val =
        rr[-1 + blhom_1.nstab] + rr[blhom_1.nstab] + ri[-1 + blhom_1.nstab] + ri[blhom_1.nstab];
    goto free_f0f1;

    // Double real leading eigenvalues (stable)
    //   (saddle, saddle-focus transition)

L2:
    if (fabs(ri[-1 + blhom_1.nstab]) > blhma_1.compzero) {
        // Computing 2nd power
        double tmp = ri[-1 + blhom_1.nstab] - ri[-1 + blhom_1.nstab - 1];
        ret_val = -(tmp*tmp);
    } else {
        // Computing 2nd power
        double tmp = rr[-1 + blhom_1.nstab] - rr[-1 + blhom_1.nstab - 1];
        ret_val = tmp*tmp;
    }
    goto free_f0f1;

    // Double real positive eigenvalues (unstable)
    //   (saddle, saddle-focus transition)

L3:
    if (fabs(ri[blhom_1.nstab]) > blhma_1.compzero) {
        // Computing 2nd power
        double tmp = ri[blhom_1.nstab] - ri[blhom_1.nstab + 1];
        ret_val = -(tmp*tmp);
    } else {
        // Computing 2nd power
        double tmp = rr[blhom_1.nstab] - rr[blhom_1.nstab + 1];
        ret_val = tmp*tmp;
    }
    goto free_f0f1;

    // Neutral saddle, saddle-focus or bi-focus (includes 1, above, also)

L4:
    ret_val = rr[-1 + blhom_1.nstab] + rr[blhom_1.nstab];
    goto free_f0f1;

    // Neutrally-divergent saddle-focus (stable eigenvalues floatcomplex)

L5:
    ret_val = rr[-1 + blhom_1.nstab] + rr[blhom_1.nstab] + rr[blhom_1.nstab - 2];
    goto free_f0f1;

    // Neutrally-divergent saddle-focus (unstable eigenvalues floatcomplex)

L6:
    ret_val = rr[-1 + blhom_1.nstab] + rr[blhom_1.nstab] + rr[blhom_1.nstab + 1];
    goto free_f0f1;

    // Three leading eigenvalues (stable)

L7:
    ret_val = rr[-1 + blhom_1.nstab] - rr[blhom_1.nstab - 3];
    goto free_f0f1;

    // Three leading eigenvalues (ustable)

L8:
    ret_val = rr[blhom_1.nstab] - rr[blhom_1.nunstab + 2];
    goto free_f0f1;

    // Local bifurcation (zero eigenvalue or Hopf): NSTAB decreases
    //  (nb. the problem becomes ill-posed after a zero of 9 or 10)

L9:
    ret_val = rr[-1 + blhom_1.nstab];
    goto free_f0f1;

    // Local bifurcation (zero eigenvalue or Hopf): NSTAB increases

L10:
    ret_val = rr[blhom_1.nstab];
    goto free_f0f1;

    // Orbit flip (with respect to leading stable direction)
    //     e.g. 1D unstable manifold

L11:
    for (int32 j = 0; j < ndm; ++j) {
        ret_val += f1[j]*vt[blhom_1.nstab + (j + 1)*(iap->ndm)];
    }
    ret_val *= exp(-par[10]*rr[-1 + blhom_1.nstab] / 2.);
    goto free_f0f1;

    // Orbit flip (with respect to leading unstable direction)
    //     e.g. 1D stable manifold

L12:
    for (int32 j = 0; j < ndm; ++j) {
        ret_val += f0[j]*vt[blhom_1.nstab + 1 + (j + 1)*(iap->ndm)];
    }
    ret_val *= exp(par[10]*rr[blhom_1.nstab] / 2.);
    goto free_f0f1;

    /* Inclination flip (critically twisted) with respect to stable manifold
     */
    //   e.g. 1D unstable manifold

L13:
    for (int64 i = 0; i < ndm; ++i) {
        ret_val += blhmu_1.pu0[ndm + i]*v[blhom_1.nstab + (i + 1)*(iap->ndm)];
    }
    ret_val *= exp(-par[10]*rr[-1 + blhom_1.nstab] / 2.);
    goto free_f0f1;

    /* Inclination flip (critically twisted) with respect to unstable manifold
     */
    //   e.g. 1D stable manifold

L14:
    for (int64 i = 0; i < ndm; ++i) {
        ret_val += blhmu_1.pu1[ndm + i]*v[blhom_1.nstab + 1 + (i + 1)*(iap->ndm)];
    }
    ret_val *= exp(par[10]*rr[blhom_1.nstab] / 2.);
    goto free_f0f1;

    // Non-central homoclinic to saddle-node (in stable manifold)

L15:
    for (int64 i = 0; i < ndm; ++i) {
        ret_val += (par[i + 11] - blhmu_1.pu1[i])*v[blhom_1.nstab + 1 + (i + 1)*(iap->ndm)];
    }
    goto free_f0f1;

    // Non-central homoclinic to saddle-node (in unstable manifold)

L16:
    for (int64 i = 0; i < ndm; ++i) {
        ret_val += (par[i + 11] - blhmu_1.pu0[i])*v[blhom_1.nstab + 1 + (i + 1)*(iap->ndm)];
    }
    goto free_f0f1;

free_f0f1:
    free(f0);
    free(f1);
    return ret_val;
}

int32
eighi(int64 isign, int64 itrans, double *rr, double *ri, double *vret, double *xequib, int64 *icp,
      double *par, int64 *ndm) {
    double *dfdp;
    double *dfdu;
    double *zz;

    dfdp = xmalloc(sizeof(*dfdp)*(usize)((*ndm)*NPARX));
    dfdu = xmalloc(sizeof(*dfdu)*(usize)((*ndm)*(*ndm)));
    zz = xmalloc(sizeof(*zz)*(usize)((*ndm)*(*ndm)));

    eigho(&isign, &itrans, rr, ri, vret, xequib, icp, par, ndm, dfdu, dfdp, zz);

    free(dfdp);
    free(dfdu);
    free(zz);

    return 0;
}

int32
eigho(int64 *isign, int64 *itrans, double *rr, double *ri, double *vret, double *xequib, int64 *icp,
      double *par, int64 *ndm, double *dfdu, double *dfdp, double *zz) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;
    int64 zz_dim1;

    int64 ifail;
    double vdot;

    double *f;
    double *ridum, *vidum, *rrdum, *vrdum;

    double *vi, *vr, *fv1;
    int64 *iv1;

    (void)isign;

    f = xmalloc(sizeof(*f)*(usize)(*ndm));
    ridum = xmalloc(sizeof(*ridum)*(usize)(*ndm));
    vidum = xmalloc(sizeof(*vidum)*(usize)((*ndm)*(*ndm)));
    rrdum = xmalloc(sizeof(*rrdum)*(usize)(*ndm));
    vrdum = xmalloc(sizeof(*vrdum)*(usize)((*ndm)*(*ndm)));
    vi = xmalloc(sizeof(*vi)*(usize)((*ndm)*(*ndm)));
    vr = xmalloc(sizeof(*vr)*(usize)((*ndm)*(*ndm)));
    fv1 = xmalloc(sizeof(*(fv1))*(usize)(*ndm));
    iv1 = xmalloc(sizeof(*(iv1))*(usize)(*ndm));

    // Uses EISPACK routine RG to calculate the eigenvalues/eigenvectors
    // of the linearization matrix a (obtained from DFHO) and orders them
    // according to their real parts. Simple continuity with respect
    // previous call with same value of ITRANS.

    // 	input variables
    // 		ISIGN  = 1 => left-hand endpoint
    //       	       = 2 => right-hand endpoint
    //               ITRANS = 1 use transpose of A
    //                      = 2 otherwise

    //       output variables
    // 		RR,RI real and imaginary parts of eigenvalues, ordered w.r.t
    // 	           real parts (largest first)
    // 	        VRET the rows of which are real parts of corresponding
    //                  eigenvectors

    // Local

    vret -= ((*ndm) + 1);
    zz_dim1 = *ndm;
    dfdp_dim1 = *ndm;
    dfdu_dim1 = *ndm;

    ifail = 0;

    func(*ndm, xequib, icp, par, 1, f, dfdu, dfdp);

    if (*itrans == 1) {
        for (int64 i = 0; i < *ndm; ++i) {
            for (int32 j = 0; j < *ndm; ++j) {
                vrdum[i + j*(*ndm)] = ARRAY2D(dfdu, j, i);
            }
        }
        for (int64 i = 0; i < *ndm; ++i) {
            for (int32 j = 0; j < *ndm; ++j) {
                ARRAY2D(dfdu, i, j) = vrdum[i + j*(*ndm)];
            }
        }
    }

    // EISPACK call for eigenvalues and eigenvectors
    rg(*ndm, *ndm, dfdu, rr, ri, 1, zz, iv1, fv1, &ifail);

    if (ifail != 0) {
        fprintf(fp9, "EISPACK EIGENVALUE ROUTINE FAILED !\n");
    }

    for (int32 j = 0; j < *ndm; ++j) {
        if (ri[j] > 0.) {
            for (int64 i = 0; i < *ndm; ++i) {
                vr[i + j*(*ndm)] = ARRAY2D(zz, i, j);
                vi[i + j*(*ndm)] = ARRAY2D(zz, i, (j + 1));
            }
        } else if (ri[j] < 0.) {
            for (int64 i = 0; i < *ndm; ++i) {
                vr[i + j*(*ndm)] = ARRAY2D(zz, i, (j - 1));
                vi[i + j*(*ndm)] = -ARRAY2D(zz, i, j);
            }
        } else {
            for (int64 i = 0; i < *ndm; ++i) {
                vr[i + j*(*ndm)] = ARRAY2D(zz, i, j);
                vi[i + j*(*ndm)] = 0.;
            }
        }
    }
    /*Order the eigenvectors/values according size of real part of eigenvalue.
     */
    //     (smallest first)

    for (int64 i = 0; i < *ndm - 1; ++i) {
        for (int64 j = i + 1; j < *ndm; ++j) {
            if (rr[i] > rr[j]) {
                rrdum[i] = rr[i];
                ridum[i] = ri[i];
                rr[i] = rr[j];
                rr[j] = rrdum[i];
                ri[i] = ri[j];
                ri[j] = ridum[i];
                for (int64 k = 0; k < *ndm; ++k) {
                    vrdum[k + i*(*ndm)] = vr[k + i*(*ndm)];
                    vr[k + i*(*ndm)] = vr[k + j*(*ndm)];
                    vr[k + j*(*ndm)] = vrdum[k + i*(*ndm)];
                    vidum[k + i*(*ndm)] = vi[k + i*(*ndm)];
                    vi[k + i*(*ndm)] = vi[k + j*(*ndm)];
                    vi[k + j*(*ndm)] = vidum[k + i*(*ndm)];
                }
            }
        }
    }

    // Choose sign of real part of eigenvectors to be
    // commensurate with that of the corresponding eigenvector
    // from the previous call with the same value of ISIGN

    if (blhme_1.ieigc[*itrans - 1] == 0) {
        for (int32 j = 0; j < *ndm; ++j) {
            for (int64 i = 0; i < *ndm; ++i) {
                blhme_1.vrprev[*itrans + (i*2 + j*(*ndm)*2) - 1] = vr[i + j*(*ndm)];
            }
        }
        blhme_1.ieigc[*itrans - 1] = 1;
    }
    for (int64 i = 0; i < *ndm; ++i) {
        vdot = 0.;
#define GCC_2_96_FIX
#ifdef GCC_2_96_FIX
        {
            int64 tmp;
            tmp = *ndm;
            for (int32 j = 0; j < tmp; ++j) {
                vdot += vr[j + i*tmp]*blhme_1.vrprev[*itrans + (j*2 + i*tmp*2) - 1];
            }
        }
#else
        for (int32 j = 0; j < *ndm; ++j) {
            vdot += vr[j + i*(*ndm)]*blhme_1.vrprev[*itrans + (j*2 + i*(*ndm)*2) - 1];
        }
#endif
        if (vdot < 0.) {
            for (int32 j = 0; j < *ndm; ++j) {
                vr[j + i*(*ndm)] = -vr[j + i*(*ndm)];
                //               VI(J,I)=-VI(J,I)
            }
        }
        for (int32 j = 0; j < *ndm; ++j) {
            blhme_1.vrprev[*itrans + (j*2 + i*(*ndm)*2) - 1] = vr[j + i*(*ndm)];
        }
    }

    // Send back the transpose of the matrix of real parts of eigenvectors
    for (int64 i = 0; i < *ndm; ++i) {
        for (int32 j = 0; j < *ndm; ++j) {
            vret[(i + 1) + (j + 1)*(*ndm)] = vr[j + i*(*ndm)];
        }
    }

    free(f);
    free(ridum);
    free(vidum);
    free(rrdum);
    free(vrdum);
    free(vi);
    free(vr);
    free(fv1);
    free(iv1);
    return 0;
}

int32
prjcti(double *bound2, double *xequib, int64 *icp, double *par, int64 imfd, int64 is, int64 itrans,
       int64 *ndm) {
    double *dfdp;
    double *dfdu;

    dfdp = xmalloc(sizeof(*dfdp)*(usize)((*ndm)*NPARX));
    dfdu = xmalloc(sizeof(*dfdu)*(usize)((*ndm)*(*ndm)));

    prjctn(bound2, xequib, icp, par, &imfd, &is, &itrans, ndm, dfdu, dfdp);

    free(dfdp);
    free(dfdu);
    return 0;
}

int32
prjctn(double *bound2, double *xequib, int64 *icp, double *par, int64 *imfd, int64 *is,
       int64 *itrans, int64 *ndm, double *dfdu, double *dfdp) {
    int64 dfdu_dim1;
    int64 dfdp_dim1;

    int64 mcond;
    int64 k1;
    int64 k2;
    int64 m0;

    double det;
    double eps;

    double *fdum;
    double *cnow;
    int64 *type__;
    double *a;
    double *d;
    double *v;
    double *ei;
    double *er;
    double *ort, *dum1, *dum2;

    fdum = xmalloc(sizeof(*fdum)*(usize)(*ndm));
    cnow = xmalloc(sizeof(*cnow)*(usize)((*ndm)*(*ndm)));
    type__ = xmalloc(sizeof(*type__)*(usize)(*ndm));
    a = xmalloc(sizeof(*a)*(usize)((*ndm)*(*ndm)));
    d = xmalloc(sizeof(*d)*(usize)((*ndm)*(*ndm)));
    v = xmalloc(sizeof(*v)*(usize)((*ndm)*(*ndm)));
    ei = xmalloc(sizeof(*ei)*(usize)(*ndm));
    er = xmalloc(sizeof(*er)*(usize)(*ndm));
    ort = xmalloc(sizeof(*ort)*(usize)(*ndm));
    dum1 = xmalloc(sizeof(*(dum1))*(usize)((*ndm)*(*ndm)));
    dum2 = xmalloc(sizeof(*(dum2))*(usize)((*ndm)*(*ndm)));

    // Compute NUNSTAB (or NSTAB) projection boundary condition functions
    /*onto to the UNSTABLE (or STABLE) manifold of the appropriate equilibrium
     */

    //    IMFD   = -1 stable eigenspace
    //           =  1 unstable eigenspace
    //    ITRANS =  1 use transpose of A
    //           =  2 otherwise
    //    IS     =  I (1 or 2) implies use the ith equilibrium in XEQUIB

    // Use the normalization in Beyn 1990 (4.4) to ensure continuity
    // w.r.t parameters.
    // For the purposes of this routine the "previous point on the
    // branch" is at the values of PAR at which the routine was last
    // called with the same values of IS and ITRANS.

    // Local

    bound2 -= ((*ndm) + 1);
    dfdp_dim1 = *ndm;
    dfdu_dim1 = *ndm;

    func(*ndm, xequib, icp, par, 1, fdum, dfdu, dfdp);

    // Compute transpose of A if ITRANS=1
    if (*itrans == 1) {
        for (int64 i = 0; i < *ndm; ++i) {
            for (int32 j = 0; j < *ndm; ++j) {
                a[i + j*(*ndm)] = ARRAY2D(dfdu, j, i);
            }
        }
    } else {
        for (int64 i = 0; i < *ndm; ++i) {
            for (int32 j = 0; j < *ndm; ++j) {
                a[i + j*(*ndm)] = ARRAY2D(dfdu, i, j);
            }
        }
    }

    // Compute basis V to put A in upper Hessenberg form
    {
        /* This is here since I don't want to change the calling sequence of the
           BLAS routines. */
        int64 tmp = 1;
        orthes((ndm), ndm, &tmp, ndm, a, ort);
        ortran((ndm), ndm, &tmp, ndm, a, ort, v);
    }

    // Force A to be upper Hessenberg
    if (*ndm > 2) {
        for (int64 i = 2; i < *ndm; ++i) {
            for (int32 j = 0; j < i - 1; ++j) {
                a[i + j*(*ndm)] = 0.;
            }
        }
    }

    // Computes basis to put A in "Quasi Upper-Triangular form"
    // with the positive (negative) eigenvalues first if IMFD =-1 (=1)
    eps = blhma_1.compzero;
    {
        /* This is here since I don't want to change the calling sequence of the
           BLAS routines. */
        int64 tmp = 1;
        hqr3lc(a, v, ndm, &tmp, ndm, &eps, er, ei, type__, (ndm), (ndm), imfd);
    }
    // Put the basis in the appropriate part of the matrix CNOW
    if (*imfd == 1) {
        k1 = *ndm - blhom_1.nunstab + 1;
        k2 = *ndm;
    } else {
        k1 = 1;
        k2 = blhom_1.nstab;
    }
    mcond = k2 - k1 + 1;
    m0 = k1 - 1;

    for (int64 i = k1 - 1; i < k2; ++i) {
        for (int32 j = 0; j < *ndm; ++j) {
            cnow[i + j*(*ndm)] = v[j + (i - k1 + 1)*(*ndm)];
        }
    }

    // Set previous matrix to be the present one if this is the first call

    /* Note by Randy Paffenroth:  There is a slight problem here
       in that this array is used before it is assigned to,
       hence its value is, in general, undefined.  It has
       worked because the just happened to be filled
       with zeros, even though this is not guaranteed.*/
    if (beyn_1.iflag[*is + (*itrans*2) - 3] == 0) {
        for (int64 i = k1 - 1; i < k2; ++i) {
            for (int32 j = 0; j < *ndm; ++j) {
                beyn_1.cprev[i + (j + ((*is - 1) + ((*itrans - 1)*2))*(*ndm))*(*ndm)] =
                    cnow[i + j*(*ndm)];
                bound2[(i + 1) + (j + 1)*(*ndm)] = cnow[i + j*(*ndm)];
            }
        }
        beyn_1.iflag[*is + (*itrans*2) - 3] = 1;
        free(fdum);
        free(cnow);
        free(type__);
        free(a);
        free(d);
        free(v);
        free(ei);
        free(er);
        free(ort);
        free(dum1);
        free(dum2);
        return 0;
    }

    // Calculate the (transpose of the) BEYN matrix D and hence bound2
    for (int64 i = 0; i < mcond; ++i) {
        for (int32 j = 0; j < mcond; ++j) {
            dum1[i + j*(*ndm)] = 0.;
            dum2[i + j*(*ndm)] = 0.;
#define GCC_2_96_FIX
#ifdef GCC_2_96_FIX
            {
                int64 tmp;
                tmp = *ndm;
                for (int64 k = 0; k < tmp; ++k) {
                    dum1[i + j*(tmp)] +=
                        beyn_1.cprev[i + m0 +
                                     (k + ((*is - 1) + ((*itrans - 1)*2))*(tmp))*(tmp)] *
                        cnow[j + m0 + k*(tmp)];
                    dum2[i + j*(tmp)] +=
                        beyn_1.cprev[i + m0 +
                                     (k + ((*is - 1) + ((*itrans - 1)*2))*(tmp))*(tmp)] *
                        beyn_1.cprev[j + m0 +
                                     (k + ((*is - 1) + ((*itrans - 1)*2))*(tmp))*(tmp)];
                }
            }
#else
            for (int64 k = 0; k < *ndm; ++k) {
                dum1[i + j*(*ndm)] +=
                    beyn_1
                        .cprev[i + m0 + (k + ((*is - 1) + ((*itrans - 1)*2))*(*ndm))*(*ndm)] *
                    cnow[j + m0 + k*(*ndm)];
                dum2[i + j*(*ndm)] +=
                    beyn_1
                        .cprev[i + m0 + (k + ((*is - 1) + ((*itrans - 1)*2))*(*ndm))*(*ndm)] *
                    beyn_1
                        .cprev[j + m0 + (k + ((*is - 1) + ((*itrans - 1)*2))*(*ndm))*(*ndm)];
            }
#endif
        }
    }

    if (mcond > 0) {
        ge(mcond, *ndm, dum1, mcond, *ndm, d, *ndm, dum2, &det);
    }

    for (int64 i = 0; i < mcond; ++i) {
        for (int32 j = 0; j < *ndm; ++j) {
            bound2[(i + 1) + m0 + (j + 1)*(*ndm)] = 0.;
            for (int64 k = 0; k < mcond; ++k) {
                bound2[(i + 1) + m0 + (j + 1)*(*ndm)] +=
                    d[k + i*(*ndm)]*cnow[k + m0 + j*(*ndm)];
            }
        }
    }

    for (int64 i = k1 - 1; i < k2; ++i) {
        for (int32 j = 0; j < *ndm; ++j) {
            beyn_1.cprev[i + (j + ((*is - 1) + ((*itrans - 1)*2))*(*ndm))*(*ndm)] =
                bound2[(i + 1) + (j + 1)*(*ndm)];
        }
    }

    free(fdum);
    free(cnow);
    free(type__);
    free(a);
    free(d);
    free(v);
    free(ei);
    free(er);
    free(ort);
    free(dum1);
    free(dum2);

    return 0;
}
