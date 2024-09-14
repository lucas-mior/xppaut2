#include "auto_f2c.h"
#include "integers.h"
#include "auto_c.h"
#include "auto_types.h"
#include <string.h>

#ifdef TIME
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

static double
time_start(void) {
    struct timeval time;
    double seconds, microseconds;
    gettimeofday(&time, NULL);
    seconds = (double)time.tv_sec;
    microseconds = (double)time.tv_usec;
    return seconds + microseconds / 1e6;
}

static double
time_end(double start) {
    struct timeval time;
    double seconds, microseconds;
    gettimeofday(&time, NULL);
    seconds = (double)time.tv_sec;
    microseconds = (double)time.tv_usec;
    return (seconds + microseconds / 1e6) - start;
}
#endif

void *
setubv_make_aa_bb_cc(void *arg) {
    /* System generated locals */
    int64 aa_dim1, aa_dim2, bb_dim1, bb_dim2, cc_dim1, cc_dim2, ups_dim1,
        uoldps_dim1, udotps_dim1, upoldp_dim1, dbc_dim1, dicd_dim1, wploc_dim1,
        dfdu_dim1, dfdp_dim1, wp_dim1, wt_dim1;

    /* Local variables */
    int64 i, j, k, l, m;
    int64 k1, l1;
    int64 i1, j1;

    int64 ib, ic, jj;
    double dt;
    int64 ib1, ic1;
    int64 jp1;
    double ddt;

    setubv_parallel_arglist *larg = (setubv_parallel_arglist *)arg;

    double *dicd, *ficd, *dfdp, *dfdu, *uold;
    double *f;
    double *u, *wploc;
    double *dbc, *fbc, *uic, *uio, *prm, *uid, *uip, *ubc0, *ubc1;

    double *ups = larg->ups;
    double *upoldp = larg->upoldp;
    double *udotps = larg->udotps;
    double *uoldps = larg->uoldps;

    double *aa = larg->aa;
    double *bb = larg->bb;
    double *cc = larg->cc;

    double *wp = larg->wp;
    double *wt = larg->wt;

    dicd = malloc(sizeof(*dicd)*(larg->nint)*(larg->ndim + NPARX));
    ficd = malloc(sizeof(*ficd)*(larg->nint));
    dfdp = malloc(sizeof(*dfdp)*(larg->ndim)*NPARX);
    dfdu = malloc(sizeof(*dfdu)*(larg->ndim)*(larg->ndim));
    uold = malloc(sizeof(*uold)*(larg->ndim));
    f = malloc(sizeof(*f)*(larg->ndim));
    u = malloc(sizeof(*u)*(larg->ndim));
    wploc = malloc(sizeof(*wploc)*(larg->ncol)*(larg->ncol + 1));
    dbc = malloc(sizeof(*dbc)*(larg->nbc)*(2*larg->ndim + NPARX));
    fbc = malloc(sizeof(*fbc)*(larg->nbc));
    uic = malloc(sizeof(*uic)*(larg->ndim));
    uio = malloc(sizeof(*uio)*(larg->ndim));
    prm = malloc(sizeof(*prm)*NPARX);
    uid = malloc(sizeof(*uid)*(larg->ndim));
    uip = malloc(sizeof(*uip)*(larg->ndim));
    ubc0 = malloc(sizeof(*(ubc0))*(larg->ndim));
    ubc1 = malloc(sizeof(*(ubc1))*(larg->ndim));

    upoldp_dim1 = larg->ndxloc;
    udotps_dim1 = larg->ndxloc;
    uoldps_dim1 = larg->ndxloc;
    ups_dim1 = larg->ndxloc;
    dicd_dim1 = larg->nint;
    dbc_dim1 = larg->nbc;
    dfdu_dim1 = larg->ndim;
    dfdp_dim1 = larg->ndim;

    bb_dim1 = larg->ncb;
    bb_dim2 = larg->nra;

    cc_dim1 = larg->nca;
    cc_dim2 = larg->nrc;

    aa_dim1 = larg->nca;
    aa_dim2 = larg->nra;

    wploc_dim1 = larg->ncol + 1;
    wp_dim1 = larg->ncol + 1;
    wt_dim1 = larg->ncol + 1;

    /* Generate AA and BB: */

    /*      Partition the mesh intervals */
    /*jj will be replaced with loop_start and loop_end*/
    for (jj = larg->loop_start; jj < larg->loop_end; ++jj) {
        j = jj;
        jp1 = j + 1;
        dt = larg->dtm[j];
        ddt = 1. / dt;
        for (ic = 0; ic < larg->ncol; ++ic) {
            for (ib = 0; ib < larg->ncol + 1; ++ib) {
                ARRAY2D(wploc, ib, ic) = ddt*ARRAY2D(wp, ib, ic);
            }
        }
        /*this loop uses the loop_offset variable since up and uoldps
          and sent by the MPI version in their entirety, but
          loop_start and loop_end have been shifted.  The loop_offset
          variable contains the original value of loop_start and removes
          the shift*/
        for (ic = 0; ic < larg->ncol; ++ic) {
            for (k = 0; k < larg->ndim; ++k) {
                u[k] = ARRAY2D(wt, larg->ncol, ic) *
                       ARRAY2D(ups, jp1 + larg->loop_offset, k);
                uold[k] = ARRAY2D(wt, larg->ncol, ic) *
                          ARRAY2D(uoldps, jp1 + larg->loop_offset, k);
                for (l = 0; l < larg->ncol; ++l) {
                    l1 = l*larg->ndim + k;
                    u[k] += ARRAY2D(wt, l, ic) *
                            ARRAY2D(ups, j + larg->loop_offset, l1);
                    uold[k] += ARRAY2D(wt, l, ic) *
                               ARRAY2D(uoldps, j + larg->loop_offset, l1);
                }
            }

            for (i = 0; i < NPARX; ++i) {
                prm[i] = larg->par[i];
            }
            /*
                Ok this is a little wierd, so hold tight.  This function
                is actually a pointer to a wrapper function, which eventually
                calls the user defined func_.  Which wrapper is used
                depends on what kind of problem it is.  The need for
                the mutex is because some of these wrappers use a common
                block for temporary storage
                NOTE!!!:  The icni and bcni wrappers do the same thing,
                so if they ever get parallelized they need to be
                checked as well.
            */
            (*(larg->funi))(larg->iap, larg->rap, larg->ndim, u, uold,
                            larg->icp, prm, 2, f, dfdu, dfdp);

            ic1 = ic*(larg->ndim);
            for (ib = 0; ib < larg->ncol + 1; ++ib) {
                double wt_tmp = ARRAY2D(wt, ib, ic);
                double wploc_tmp = ARRAY2D(wploc, ib, ic);
                ib1 = ib*larg->ndim;
                for (i = 0; i < larg->ndim; ++i) {
                    ARRAY3D(aa, ib1 + i, ic1 + i, jj) = wploc_tmp;
                    for (k = 0; k < larg->ndim; ++k) {
                        ARRAY3D(aa, ib1 + k, ic1 + i, jj) -=
                            wt_tmp*ARRAY2D(dfdu, i, k);
                    }
                }
            }
            for (i = 0; i < larg->ndim; ++i) {
                for (k = 0; k < larg->ncb; ++k) {
                    ARRAY3D(bb, k, ic1 + i, jj) =
                        -ARRAY2D(dfdp, i, larg->icp[k]);
                }
            }
        }
    }

    /*     Generate CC : */

    /*     Boundary conditions : */
    if (larg->nbc > 0) {
        for (i = 0; i < larg->ndim; ++i) {
            ubc0[i] = ARRAY2D(ups, 0, i);
            ubc1[i] = ARRAY2D(ups, larg->na, i);
        }

        (*(larg->bcni))(larg->iap, larg->rap, larg->ndim, larg->par, larg->icp,
                        larg->nbc, ubc0, ubc1, fbc, 2, dbc);
        for (i = 0; i < larg->nbc; ++i) {
            for (k = 0; k < larg->ndim; ++k) {
                /*NOTE!!
                  This needs to split up.  Only the first processor does the
                  first part and only the last processors does the last part.*/
                if (larg->loop_offset + larg->loop_start == 0) {
                    ARRAY3D(cc, k, i, 0) = ARRAY2D(dbc, i, k);
                }
                if (larg->loop_offset + larg->loop_end == larg->na) {
                    ARRAY3D(cc, larg->nra + k, i,
                            larg->na - 1 - larg->loop_offset) =
                        ARRAY2D(dbc, i, larg->ndim + k);
                }
            }
        }
    }

    /*     Integral constraints : */
    if (larg->nint > 0) {
        for (jj = larg->loop_start; jj < larg->loop_end; ++jj) {
            j = jj;
            jp1 = j + 1;
            for (k = 0; k < (larg->ncol + 1); ++k) {
                for (i = 0; i < larg->ndim; ++i) {
                    i1 = k*larg->ndim + i;
                    j1 = j;
                    if (k + 1 == (larg->ncol + 1)) {
                        i1 = i;
                    }
                    if (k + 1 == (larg->ncol + 1)) {
                        j1 = jp1;
                    }
                    uic[i] = ARRAY2D(ups, j1 + larg->loop_offset, i1);
                    uio[i] = ARRAY2D(uoldps, j1 + larg->loop_offset, i1);
                    uid[i] = ARRAY2D(udotps, j1 + larg->loop_offset, i1);
                    uip[i] = ARRAY2D(upoldp, j1 + larg->loop_offset, i1);
                }

                (*(larg->icni))(larg->iap, larg->rap, larg->ndim, larg->par,
                                larg->icp, larg->nint, uic, uio, uid, uip, ficd,
                                2, dicd);

                for (m = 0; m < larg->nint; ++m) {
                    for (i = 0; i < larg->ndim; ++i) {
                        k1 = k*larg->ndim + i;
                        ARRAY3D(cc, k1, larg->nbc + m, jj) =
                            larg->dtm[j]*larg->wi[k]*ARRAY2D(dicd, m, i);
                    }
                }
            }
        }
    }
    /*     Pseudo-arclength equation : */
    for (jj = larg->loop_start; jj < larg->loop_end; ++jj) {
        for (i = 0; i < larg->ndim; ++i) {
            for (k = 0; k < larg->ncol; ++k) {
                k1 = k*larg->ndim + i;
                ARRAY3D(cc, k1, larg->nrc - 1, jj) =
                    larg->dtm[jj]*larg->thu[i]*larg->wi[k] *
                    ARRAY2D(udotps, jj + larg->loop_offset, k1);
            }
            ARRAY3D(cc, larg->nra + i, larg->nrc - 1, jj) =
                larg->dtm[jj]*larg->thu[i]*larg->wi[larg->ncol] *
                ARRAY2D(udotps, jj + 1 + larg->loop_offset, i);
        }
    }

    free(dicd);
    free(ficd);
    free(dfdp);
    free(dfdu);
    free(uold);
    free(f);
    free(u);
    free(wploc);
    free(dbc);
    free(fbc);
    free(uic);
    free(uio);
    free(prm);
    free(uid);
    free(uip);
    free(ubc0);
    free(ubc1);

    return NULL;
}

int32
setubv_default_wrapper(setubv_parallel_arglist data) {
    setubv_make_aa_bb_cc((void *)&data);
    return 0;
}

int32
setubv(int64 ndim, int64 ips, int64 na, int64 ncol, int64 nbc, int64 nint,
       int64 ncb, int64 nrc, int64 nra, int64 nca, FUNI_TYPE((*funi)),
       BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), int64 ndxloc, iap_type *iap,
       rap_type *rap, double *par, int64 *icp, double rds, double *aa,
       double *bb, double *cc, double *dd, double *fa, double *fc,
       double *rlcur, double *rlold, double *rldot, double *ups, double *uoldps,
       double *udotps, double *upoldp, double *dups, double *dtm, double *thl,
       double *thu, double *p0, double *p1) {
    (void) p0;
    (void) p1;
    /* System generated locals */
    int64 aa_dim1, aa_dim2, bb_dim1, bb_dim2, cc_dim1, cc_dim2, dd_dim1;

    /* Local variables */
    int64 i, j, k;

    double *wi, *wp, *wt;

    wi = malloc(sizeof(*wi)*(ncol + 1));
    wp = malloc(sizeof(*wp)*(ncol)*(ncol + 1));
    wt = malloc(sizeof(*wt)*(ncol)*(ncol + 1));

    dd_dim1 = ncb;

    bb_dim1 = ncb;
    bb_dim2 = nra;

    cc_dim1 = nca;
    cc_dim2 = nrc;

    aa_dim1 = nca;
    aa_dim2 = nra;

    wint(ncol + 1, wi);
    genwts(ncol, ncol + 1, wt, wp);

    /* Initialize to zero. */
    for (i = 0; i < nrc; ++i) {
        fc[i] = 0.;
        for (k = 0; k < ncb; ++k) {
            ARRAY2D(dd, k, i) = 0.;
        }
    }

    /* Set constants. */
    for (i = 0; i < ncb; ++i) {
        par[icp[i]] = rlcur[i];
    }

    /*  NA is the local node's mesh interval number. */

    for (i = 0; i < na; ++i) {
        for (j = 0; j < nra; ++j) {
            for (k = 0; k < nca; ++k) {
                ARRAY3D(aa, k, j, i) = 0.;
            }
        }
        for (j = 0; j < nra; ++j) {
            for (k = 0; k < ncb; ++k) {
                ARRAY3D(bb, k, j, i) = 0.;
            }
        }
        for (j = 0; j < nca; ++j) {
            for (k = 0; k < nrc; ++k) {
                ARRAY3D(cc, j, k, i) = 0.;
            }
        }
    }

    /*     ** Time evolution computations (parabolic systems) */
    if (ips == 14 || ips == 16) {
        rap->tivp = rlold[0];
    }

    {
        setubv_parallel_arglist arglist;
        setubv_parallel_arglist_constructor(
            ndim, ips, na, ncol, nbc, nint, ncb, nrc, nra, nca, funi, icni,
            ndxloc, iap, rap, par, icp, aa, bb, cc, dd, fa, fc, ups, uoldps,
            udotps, upoldp, dtm, wp, wt, wi, thu, thl, rldot, bcni, &arglist);

        switch (global_setubv_type) {
        default:
            setubv_default_wrapper(arglist);
            break;
        }
        setubv_make_fa(arglist);
        setubv_make_fc_dd(arglist, dups, rlcur, rlold, rds);
    }

    free(wi);
    free(wp);
    free(wt);
    return 0;
}

void
setubv_make_fa(setubv_parallel_arglist larg) {
    int64 i, j, k, l;
    int64 ic, k1, ib;
    int64 jj, jp1, l1, ic1;
    double dt, ddt;

    double *ups = larg.ups;
    int64 ups_dim1 = larg.ndxloc;

    double *uoldps = larg.uoldps;
    int64 uoldps_dim1 = larg.ndxloc;

    double *wp = larg.wp;
    int64 wp_dim1 = larg.ncol + 1;

    double *wt = larg.wt;
    int64 wt_dim1 = larg.ncol + 1;

    double *fa = larg.fa;
    int64 fa_dim1 = larg.nra;

    double *wploc = malloc(sizeof(double)*(larg.ncol)*(larg.ncol + 1));
    int64 wploc_dim1 = larg.ncol + 1;

    double *dfdp = malloc(sizeof(double)*(larg.ndim)*NPARX);
    double *dfdu = malloc(sizeof(double)*(larg.ndim)*(larg.ndim));
    double *u = malloc(sizeof(double)*(larg.ndim));
    double *uold = malloc(sizeof(double)*(larg.ndim));
    double *f = malloc(sizeof(double)*(larg.ndim));
    double *prm = malloc(sizeof(double)*NPARX);

    for (jj = 0; jj < larg.na; ++jj) {
        j = jj;
        jp1 = j + 1;
        dt = larg.dtm[j];
        ddt = 1. / dt;
        for (ic = 0; ic < larg.ncol; ++ic) {
            for (ib = 0; ib < larg.ncol + 1; ++ib) {
                ARRAY2D(wploc, ib, ic) = ddt*ARRAY2D(wp, ib, ic);
            }
        }
        for (ic = 0; ic < larg.ncol; ++ic) {
            for (k = 0; k < larg.ndim; ++k) {
                u[k] = ARRAY2D(wt, larg.ncol, ic)*ARRAY2D(ups, jp1, k);
                uold[k] = ARRAY2D(wt, larg.ncol, ic)*ARRAY2D(uoldps, jp1, k);
                for (l = 0; l < larg.ncol; ++l) {
                    l1 = l*larg.ndim + k;
                    u[k] += ARRAY2D(wt, l, ic) *
                            ARRAY2D(ups, j + larg.loop_offset, l1);
                    uold[k] += ARRAY2D(wt, l, ic) *
                               ARRAY2D(uoldps, j + larg.loop_offset, l1);
                }
            }

            for (i = 0; i < NPARX; ++i) {
                prm[i] = larg.par[i];
            }
            (*(larg.funi))(larg.iap, larg.rap, larg.ndim, u, uold, larg.icp,
                           prm, 2, f, dfdu, dfdp);

            ic1 = ic*(larg.ndim);
            for (i = 0; i < larg.ndim; ++i) {
                ARRAY2D(fa, ic1 + i, jj) =
                    f[i] - ARRAY2D(wploc, larg.ncol, ic) *
                               ARRAY2D(ups, jp1 + larg.loop_offset, i);
                for (k = 0; k < larg.ncol; ++k) {
                    k1 = k*larg.ndim + i;
                    ARRAY2D(fa, ic1 + i, jj) -=
                        ARRAY2D(wploc, k, ic) *
                        ARRAY2D(ups, j + larg.loop_offset, k1);
                }
            }
        }
    }
    free(wploc);
    free(dfdp);
    free(dfdu);
    free(u);
    free(uold);
    free(f);
    free(prm);
    return;
}

void
setubv_make_fc_dd(setubv_parallel_arglist larg, double *dups, double *rlcur,
                  double *rlold, double rds) {
    int64 i, j, jj, jp1, k, i1, m, j1;
    double rlsum;

    int64 dups_dim1 = larg.ndxloc;

    double *dd = larg.dd;
    int64 dd_dim1 = larg.ncb;

    double *ups = larg.ups;
    int64 ups_dim1 = larg.ndxloc;

    double *uoldps = larg.uoldps;
    int64 uoldps_dim1 = larg.ndxloc;

    double *udotps = larg.udotps;
    int64 udotps_dim1 = larg.ndxloc;

    double *upoldp = larg.upoldp;
    int64 upoldp_dim1 = larg.ndxloc;

    int64 dbc_dim1 = larg.nbc;
    double *dbc = malloc(sizeof(double)*(larg.nbc)*(2*larg.ndim + NPARX));
    double *fbc = malloc(sizeof(double)*(larg.nbc));
    double *ubc0 = malloc(sizeof(double)*(larg.ndim));
    double *ubc1 = malloc(sizeof(double)*(larg.ndim));
    int64 dicd_dim1 = larg.nint;
    double *dicd = malloc(sizeof(double)*(larg.nint)*(larg.ndim + NPARX));
    double *ficd = malloc(sizeof(double)*(larg.nint));
    double *uic = malloc(sizeof(double)*(larg.ndim));
    double *uio = malloc(sizeof(double)*(larg.ndim));
    double *uid = malloc(sizeof(double)*(larg.ndim));
    double *uip = malloc(sizeof(double)*(larg.ndim));

    /* Boundary condition part of FC */
    if (larg.nbc > 0) {
        for (i = 0; i < larg.ndim; ++i) {
            ubc0[i] = ARRAY2D(ups, 0, i);
            ubc1[i] = ARRAY2D(ups, larg.na, i);
        }

        (*(larg.bcni))(larg.iap, larg.rap, larg.ndim, larg.par, larg.icp,
                       larg.nbc, ubc0, ubc1, fbc, 2, dbc);
        for (i = 0; i < larg.nbc; ++i) {
            larg.fc[i] = -fbc[i];
            for (k = 0; k < larg.ncb; ++k) {
                ARRAY2D(dd, k, i) =
                    ARRAY2D(dbc, i, (larg.ndim*2) + larg.icp[k]);
            }
        }
        /*       Save difference : */
        for (j = 0; j < larg.na + 1; ++j) {
            for (i = 0; i < larg.nra; ++i) {
                ARRAY2D(dups, j, i) =
                    ARRAY2D(ups, j, i) - ARRAY2D(uoldps, j, i);
            }
        }
    }

    /* Integral constraint part of FC */
    if (larg.nint > 0) {
        for (jj = larg.loop_start; jj < larg.loop_end; ++jj) {
            j = jj;
            jp1 = j + 1;
            for (k = 0; k < (larg.ncol + 1); ++k) {
                for (i = 0; i < larg.ndim; ++i) {
                    i1 = k*larg.ndim + i;
                    j1 = j;
                    if (k + 1 == (larg.ncol + 1)) {
                        i1 = i;
                    }
                    if (k + 1 == (larg.ncol + 1)) {
                        j1 = jp1;
                    }
                    uic[i] = ARRAY2D(ups, j1, i1);
                    uio[i] = ARRAY2D(uoldps, j1, i1);
                    uid[i] = ARRAY2D(udotps, j1, i1);
                    uip[i] = ARRAY2D(upoldp, j1, i1);
                }

                (*(larg.icni))(larg.iap, larg.rap, larg.ndim, larg.par,
                               larg.icp, larg.nint, uic, uio, uid, uip, ficd, 2,
                               dicd);

                for (m = 0; m < larg.nint; ++m) {
                    larg.fc[larg.nbc + m] -= larg.dtm[j]*larg.wi[k]*ficd[m];
                    for (i = 0; i < larg.ncb; ++i) {
                        ARRAY2D(dd, i, larg.nbc + m) +=
                            larg.dtm[j]*larg.wi[k] *
                            ARRAY2D(dicd, m, larg.ndim + larg.icp[i]);
                    }
                }
            }
        }
    }

    for (i = 0; i < larg.ncb; ++i) {
        ARRAY2D(dd, i, (larg.nrc - 1)) = larg.thl[larg.icp[i]]*larg.rldot[i];
    }

    rlsum = 0.;
    for (i = 0; i < larg.ncb; ++i) {
        rlsum += larg.thl[larg.icp[i]]*(rlcur[i] - rlold[i])*larg.rldot[i];
    }

    larg.fc[larg.nrc - 1] = rds -
                            rinpr(larg.iap, &(larg.ndim), &(larg.ndxloc),
                                  larg.udotps, dups, larg.dtm, larg.thu) -
                            rlsum;

    free(dbc);
    free(fbc);
    free(ubc0);
    free(ubc1);
    free(dicd);
    free(ficd);
    free(uic);
    free(uio);
    free(uid);
    free(uip);
    return;
}

/* Copy a setubv_parallel_arglist */
void
setubv_parallel_arglist_copy(setubv_parallel_arglist *output,
                             const setubv_parallel_arglist input) {
    memcpy(output, &input, sizeof(setubv_parallel_arglist));
    return;
}

/* Fill in a setubv_parallel_arglist for the individual variables */
void
setubv_parallel_arglist_constructor(
    int64 ndim, int64 ips, int64 na, int64 ncol, int64 nbc, int64 nint,
    int64 ncb, int64 nrc, int64 nra, int64 nca, FUNI_TYPE((*funi)),
    ICNI_TYPE((*icni)), int64 ndxloc, iap_type *iap, rap_type *rap, double *par,
    int64 *icp, double *aa, double *bb, double *cc, double *dd, double *fa,
    double *fc, double *ups, double *uoldps, double *udotps, double *upoldp,
    double *dtm, double *wp, double *wt, double *wi, double *thu, double *thl,
    double *rldot, BCNI_TYPE((*bcni)), setubv_parallel_arglist *data) {
    data->ndim = ndim;
    data->ips = ips;
    data->ncol = ncol;
    data->nbc = nbc;
    data->nint = nint;
    data->ncb = ncb;
    data->nrc = nrc;
    data->nra = nra;
    data->nca = nca;
    data->na = na;
    data->funi = funi;
    data->icni = icni;
    data->ndxloc = ndxloc;
    data->iap = iap;
    data->rap = rap;
    data->par = par;
    data->icp = icp;
    data->aa = aa;
    data->bb = bb;
    data->cc = cc;
    data->dd = dd;
    data->fa = fa;
    data->fc = fc;
    data->ups = ups;
    data->uoldps = uoldps;
    data->udotps = udotps;
    data->upoldp = upoldp;
    data->dtm = dtm;
    data->loop_start = 0;
    data->loop_end = na;
    data->loop_offset = 0;
    data->wp = wp;
    data->wt = wt;
    data->wi = wi;
    data->thu = thu;
    data->thl = thl;
    data->rldot = rldot;
    data->bcni = bcni;
    return;
}
