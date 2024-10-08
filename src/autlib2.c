/* Autlib2.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
        -lf2c -lm   (in that order)
*/

#include "auto_f2c.h"
#include "somemath.h"
#include "auto_c.h"
#include "integers.h"
#include "xmalloc.h"

#ifdef ACCES_TEST
struct {
    double *a, *b, *c;
} test;
#endif

typedef struct {
    double *a;
    double *b;
    double *c;
    double *d;
    double *a1;
    double *a2;
    double *s1;
    double *s2;
    double *bb;
    double *cc;
    double *faa;
    double *ca1;

    int64 *icf;
    int64 *irf;
    int64 *ipr;
    int64 *icf11;
    int64 *icf1;
    int64 *icf2;
    int64 *np;
} MainAutoStorage;

static MainAutoStorage mas = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                              NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

static int32 gdsum(void);
static int32 gsendx(void);
static int32 infpar(int64 *iam, int64 *par, double *a, double *b, double *fa, double *sol1,
                    double *sol2, double *fc, int64 *na, int64 *nov, int64 *nra, int64 *nca,
                    int64 *ncb, int64 *irf, int64 *icf);
int32 cpyrhs(int64 *na, int64 *nov, int64 *nra, double *faa, double *fa, int64 *irf);

/* ----------------------------------------------------------------------- */
/*           Setting up of the Jacobian and right hand side */
/* ----------------------------------------------------------------------- */

int32
solvbv(int64 *ifst, iap_type *iap, rap_type *rap, double *par, int64 *icp, FUNI_TYPE((*funi)),
       BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), double *rds, int64 *nllv, double *rlcur,
       double *rlold, double *rldot, int64 *ndxloc, double *ups, double *dups, double *uoldps,
       double *udotps, double *upoldp, double *dtm, double *fa, double *fc, double *p0, double *p1,
       double *thl, double *thu) {
    int64 ndim;
    int64 ipar;
    int64 ncol;
    int64 nclm;
    int64 nfpr;
    int64 nint;
    int64 nrow;
    int64 ntst;
    int64 ntst0;

    double *ff;
    double *ft;

    int64 nbc;
    int64 iid;
    int64 iam;
    double det;
    int64 ips;
    int64 nrc;

    int64 kwt;

    //     N AX is the local N TSTX, which is smaller than the global N TSTX.
    //     NODES is the total number of nodes.

    /* Sets up and solves the linear equations for one Newton/Chord iteration
     */

    // Most of the required memory is allocated below
    /* This is an interesting section of code.  The main point
       is that setubv and conpar only get called when ifst
       is 1.  This is a optimization since you can solve
       the system using the previously factored jacobian.
       One thing to watch out for is that two seperate calls
       of solvbv_ talk to each other through these arrays,
       so it is only safe to get rid of them when ifst is
       1 (since their entries will then be recreated in conpar
       and setubv).
    */

    ff = xmalloc(sizeof(*ff)*(usize)((iap->ndim*iap->ncol)*iap->ntst + 1));
    ft = xmalloc(sizeof(*ft)*(usize)((iap->ndim*iap->ncol)*(iap->ntst + 1)));

    ndim = iap->ndim;
    ips = iap->ips;
    ntst = iap->ntst;
    ncol = iap->ncol;
    nbc = iap->nbc;
    nint = iap->nint;
    iid = iap->iid;
    nfpr = iap->nfpr;
    nrc = nbc + nint + 1;
    nrow = ndim*ncol;
    nclm = nrow + ndim;

    if (*ifst == 1) {
        /* The formulas used for the allocation are somewhat floatcomplex, but
           they are based on following macros (the space after the first letter
           is for the scripts which detect these things automatically, the
           original name does not have the space:

           M 1AAR =  (((iap->ndim*iap->ncol ) + iap->ndim ) )
           M 2AA  =	((iap->ndim*iap->ncol ) )
           N AX   =	(iap->ntst /NODES+1)
           M 1BB  =	(NPARX)
           M 2BB  =	((iap->ndim*iap->ncol ) )
           M 1CC  =	((((iap->ndim*iap->ncol ) + iap->ndim ) ) )
           M 2CC  =	(((iap->ndim +3) +NINTX+1) )
           M 1DD  =	(((iap->ndim +3) +NINTX+1) )
           M 2DD  =	(NPARX)
           N RCX  =	((iap->ndim +3) +NINTX+1)
           N CLMX =	((iap->ndim*iap->ncol ) + iap->ndim )
           N ROWX =	(iap->ndim*iap->ncol )
        */

        // Free floating point arrays
        free(mas.a);
        free(mas.b);
        free(mas.c);
        free(mas.d);
        free(mas.a1);
        free(mas.a2);
        free(mas.s1);
        free(mas.s2);
        free(mas.bb);
        free(mas.cc);
        free(mas.faa);
        free(mas.ca1);

        // Free int64 arrays
        free(mas.icf);
        free(mas.irf);
        free(mas.ipr);
        free(mas.icf11);
        free(mas.icf1);
        free(mas.icf2);
        free(mas.np);

        mas.a =
            xmalloc(sizeof(*(mas.a))*(usize)((ndim*ncol + ndim)*(ndim*ncol)*(ntst + 1)));
        mas.b = xmalloc(sizeof(*(mas.b))*(usize)(NPARX*(ndim*ncol)*(ntst + 1)));
        mas.c = xmalloc(sizeof(*(mas.c)) *
                        (usize)((ndim*ncol + ndim)*(nbc + nint + 1)*(ntst + 1)));
        mas.d = xmalloc(sizeof(*(mas.d))*(usize)((nbc + nint + 1)*NPARX));
        mas.a1 = xmalloc(sizeof(*(mas.a1))*(usize)(ndim*ndim*(ntst + 1)));
        mas.a2 = xmalloc(sizeof(*(mas.a2))*(usize)(ndim*ndim*(ntst + 1)));
        mas.s1 = xmalloc(sizeof(*(mas.s1))*(usize)(ndim*ndim*(ntst + 1)));
        mas.s2 = xmalloc(sizeof(*(mas.s2))*(usize)(ndim*ndim*(ntst + 1)));
        mas.bb = xmalloc(sizeof(*(mas.bb))*(usize)(ndim*NPARX*(ntst + 1)));
        mas.cc = xmalloc(sizeof(*(mas.cc))*(usize)((nbc + nint + 1)*ndim*(ntst + 1) + 1));
        mas.faa = xmalloc(sizeof(*(mas.faa))*(usize)(ndim*(ntst + 1)));
        mas.ca1 = xmalloc(sizeof(*(mas.ca1))*(usize)(ndim*ndim*KREDO));
        mas.icf = xmalloc(sizeof(*(mas.icf))*(usize)((ndim*ncol + ndim)*(ntst + 1)));
        mas.irf = xmalloc(sizeof(*(mas.irf))*(usize)(ndim*ncol*(ntst + 1)));
        mas.ipr = xmalloc(sizeof(*(mas.ipr))*(usize)(ndim*(ntst + 1)));
        mas.icf11 = xmalloc(sizeof(*(mas.icf11))*(usize)(ndim*KREDO));
        mas.icf1 = xmalloc(sizeof(*(mas.icf1))*(usize)(ndim*(ntst + 1)));
        mas.icf2 = xmalloc(sizeof(*(mas.icf2))*(usize)(ndim*(ntst + 1)));
        mas.np = xmalloc(sizeof(*(mas.np))*(2));
    }

    iam = iap->mynode;
    kwt = iap->numnodes;
    if (kwt > 1) {
        ipar = true;
    } else {
        ipar = false;
    }

    if (kwt > ntst) {
        printf("NTST is less than the number of nodes\n");
        exit(0);
    } else {
        partition(&ntst, &kwt, mas.np);
    }

    //     NTST0 is the global one, NTST is the local one.
    //     The value of NTST may be different in different nodes.
    ntst0 = ntst;
    ntst = mas.np[iam];

    if (*ifst == 1) {
        setubv(ndim, ips, ntst, ncol, nbc, nint, nfpr, nrc, nrow, nclm, funi, bcni, icni, *ndxloc,
               iap, rap, par, icp, *rds, mas.a, mas.b, mas.c, mas.d, ft, fc, rlcur, rlold, rldot,
               ups, uoldps, udotps, upoldp, dups, dtm, thl, thu, p0, p1);
#ifdef ACCES_TEST
        test.a = mas.a;
        test.b = mas.b;
        test.c = mas.c;
        mas.a = NULL;
        mas.b = NULL;
        mas.c = NULL;
#endif
    } else {
        setrhs(&ndim, &ips, &ntst, &ntst0, mas.np, &ncol, &nbc, &nint, &nfpr, &nrc, &nrow, &nclm,
               &iam, &kwt, &ipar, funi, bcni, icni, ndxloc, iap, rap, par, icp, rds, ft, fc, rlcur,
               rlold, rldot, ups, uoldps, udotps, upoldp, dups, dtm, thl, thu, p0, p1);
    }
    /*     The matrix D and FC are set to zero for all nodes except the first.
     */
    if (iam > 0) {
        setfcdd(ifst, mas.d, fc, &nfpr, &nrc);
    }

#ifdef MATLAB_OUTPUT
    print_jacobian(*iap, mas);
    {
        static num_calls = 0;
        char filename[80];
        sprintf(filename, "before%03d", num_calls);
        num_calls++;
        print_fa_fc(*iap, ft, fc, filename);
    }
#endif
    brbd(mas.a, mas.b, mas.c, mas.d, ft, fc, p0, p1, ifst, &iid, nllv, &det, &ndim, &ntst, &nbc,
         &nrow, &nclm, &nfpr, &nrc, &iam, &kwt, &ipar, mas.a1, mas.a2, mas.bb, mas.cc, mas.faa,
         mas.ca1, mas.s1, mas.s2, mas.icf11, mas.ipr, mas.icf1, mas.icf2, mas.irf, mas.icf);
#ifdef ACCES_TEST
    mas.a = test.a;
    mas.b = test.b;
    mas.c = test.c;
#endif

    /*
      This is some stuff from the parallel version that isn't needed anymore
      ----------------------------------------------------------------------
      lenft = ntst*nrow << 3;
      lenff = ntst0*nrow << 3;
      jtmp1 = M 2AA;   I added spaces so these don't get flagged as header file
      macro dependancies jtmp2 = M 3AA;   I added spaces so these don't get
      flagged as header file macro dependancies lenff2 = jtmp1*(jtmp2 + 1) <<
      3;
    */
    if (ipar) {
        //        Global concatenation of the solution from each node.
        int64 tmp;
        gcol();
        tmp = iap->ntst + 1;
        faft(ff, fa, &ntst0, &nrow, ndxloc);
    } else {
        int64 tmp;
        tmp = iap->ntst + 1;
        faft(ft, fa, &ntst0, &nrow, ndxloc);
    }
#ifdef MATLAB_OUTPUT
    {
        static num_calls = 0;
        char filename[80];
        sprintf(filename, "after%03d", num_calls);
        num_calls++;
        print_fa_fc(*iap, ft, fc, filename);
    }
#endif

    rap->det = det;
    free(ff);
    free(ft);
    return 0;
}

int32
setfcdd(int64 *ifst, double *dd, double *fc, int64 *ncb, int64 *nrc) {
    int64 dd_dim1;

    dd_dim1 = *ncb;

    for (int64 i = 0; i < *nrc; ++i) {
        if (*ifst == 1) {
            for (int64 j = 0; j < *ncb; ++j) {
                ARRAY2D(dd, j, i) = 0.;
            }
        }
        fc[i] = 0.;
    }

    return 0;
}

int32
faft(double *ff, double *fa, int64 *ntst, int64 *nrow, int64 *ndxloc) {
    int64 fa_dim1;
    int64 ff_dim1;

    ff_dim1 = *nrow;
    fa_dim1 = *ndxloc;

    for (int64 i = 0; i < *ntst; ++i) {
        for (int64 j = 0; j < *nrow; ++j) {
            ARRAY2D(fa, i, j) = ARRAY2D(ff, j, i);
        }
    }

    return 0;
}

int32
partition(int64 *n, int64 *kwt, int64 *m) {
    int64 s;
    int64 t;

    //     Linear distribution of NTST over all nodes

    t = *n / *kwt;
    s = *n % *kwt;

    for (int64 i = 0; i < *kwt; ++i) {
        m[i] = t;
    }

    for (int64 i = 0; i < s; ++i) {
        ++m[i];
    }

    return 0;
}

int64
mypart(int64 *iam, int64 *np) {
    int64 ret_val;

    int64 k;

    //     Partition the mesh

    k = 0;
    for (int64 i = 0; i < *iam; ++i) {
        k += np[i];
    }
    ret_val = k;

    return ret_val;
}

int32
setrhs(int64 *ndim, int64 *ips, int64 *na, int64 *ntst, int64 *np, int64 *ncol, int64 *nbc,
       int64 *nint, int64 *ncb, int64 *nrc, int64 *nra, int64 *nca, int64 *iam, int64 *kwt,
       int64 *ipar, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), int64 *ndxloc,
       iap_type *iap, rap_type *rap, double *par, int64 *icp, double *rds, double *fa, double *fc,
       double *rlcur, double *rlold, double *rldot, double *ups, double *uoldps, double *udotps,
       double *upoldp, double *dups, double *dtm, double *thl, double *thu, double *p0,
       double *p1) {
    int64 ups_dim1;
    int64 dups_dim1;
    int64 uoldps_dim1;
    int64 udotps_dim1;
    int64 upoldp_dim1;
    int64 fa_dim1;
    int64 wt_dim1;
    int64 wp_dim1;
    int64 wploc_dim1;

    int64 j;
    int64 l;
    int64 m;
    int64 mpart;
    int64 i1;
    int64 j1;
    int64 k1;
    int64 l1;

    double rlsum;
    int64 ib;
    int64 ic;
    int64 jj;
    int64 ic1;

    int64 jp1;
    int64 ncp1;
    double dt;
    double ddt;

    double *dicd, *ficd, *dfdp, *dfdu, *uold;
    double *f;
    double *u;
    double *wploc;
    double *wi, *wp, *wt;
    double *dbc, *fbc, *uic, *uio, *prm, *uid, *uip, *ubc0, *ubc1;

    (void)nbc;
    (void)nca;
    (void)p0;
    (void)p1;

    dicd = xmalloc(sizeof(*dicd)*(usize)((iap->nint)*(iap->ndim + NPARX)));
    ficd = xmalloc(sizeof(*ficd)*(usize)((iap->nint)));
    dfdp = xmalloc(sizeof(*dfdp)*(usize)((iap->ndim)*NPARX));
    dfdu = xmalloc(sizeof(*dfdu)*(usize)((iap->ndim)*(iap->ndim)));
    uold = xmalloc(sizeof(*uold)*(usize)((iap->ndim)));
    f = xmalloc(sizeof(*f)*(usize)(iap->ndim));
    u = xmalloc(sizeof(*u)*(usize)(iap->ndim));
    wploc = xmalloc(sizeof(*wploc)*(usize)((iap->ncol)*(iap->ncol + 1)));
    wi = xmalloc(sizeof(*wi)*(usize)(iap->ncol + 1));
    wp = xmalloc(sizeof(*wp)*(usize)((iap->ncol)*(iap->ncol + 1)));
    wt = xmalloc(sizeof(*wt)*(usize)((iap->ncol)*(iap->ncol + 1)));
    dbc = xmalloc(sizeof(*dbc)*(usize)((iap->nbc)*(2*iap->ndim + NPARX)));
    fbc = xmalloc(sizeof(*fbc)*(usize)(iap->nbc));
    uic = xmalloc(sizeof(*uic)*(usize)(iap->ndim));
    uio = xmalloc(sizeof(*uio)*(usize)(iap->ndim));
    prm = xmalloc(sizeof(*prm)*NPARX);
    uid = xmalloc(sizeof(*uid)*(usize)(iap->ndim));
    uip = xmalloc(sizeof(*uip)*(usize)(iap->ndim));
    ubc0 = xmalloc(sizeof(*(ubc0))*(usize)(iap->ndim));
    ubc1 = xmalloc(sizeof(*(ubc1))*(usize)(iap->ndim));

    fa_dim1 = *nra;
    dups_dim1 = *ndxloc;
    upoldp_dim1 = *ndxloc;
    udotps_dim1 = *ndxloc;
    uoldps_dim1 = *ndxloc;
    ups_dim1 = *ndxloc;
    wt_dim1 = iap->ncol + 1;
    wp_dim1 = iap->ncol + 1;
    wploc_dim1 = iap->ncol + 1;

    *iam = iap->mynode;
    *kwt = iap->numnodes;
    if (*kwt > 1) {
        *ipar = true;
    } else {
        *ipar = false;
    }

    wint(*ncol + 1, wi);
    genwts(*ncol, iap->ncol + 1, wt, wp);
    // Initialize to zero.
    for (int64 i = 0; i < *nrc; ++i) {
        fc[i] = 0.;
    }

    // Set constants.
    ncp1 = *ncol + 1;
    for (int64 i = 0; i < *ncb; ++i) {
        par[icp[i]] = rlcur[i];
    }

    // Generate FA :

    //      Partition the mesh intervals.
    mpart = mypart(iam, np);

    for (jj = 0; jj < *na; ++jj) {
        j = jj + mpart;
        jp1 = j + 1;
        dt = dtm[j];
        ddt = 1. / dt;
        for (ic = 0; ic < *ncol; ++ic) {
            for (ib = 0; ib < ncp1; ++ib) {
                ARRAY2D(wploc, ib, ic) = ddt*ARRAY2D(wp, ib, ic);
            }
        }
        for (ic = 0; ic < *ncol; ++ic) {
            for (int32 k = 0; k < *ndim; ++k) {
                u[k] = ARRAY2D(wt, *ncol, ic)*ARRAY2D(ups, jp1, k);
                uold[k] = ARRAY2D(wt, *ncol, ic)*ARRAY2D(uoldps, jp1, k);
                for (l = 0; l < *ncol; ++l) {
                    l1 = l**ndim + k;
                    u[k] += ARRAY2D(wt, l, ic)*ARRAY2D(ups, j, l1);
                    uold[k] += ARRAY2D(wt, l, ic)*ARRAY2D(uoldps, j, l1);
                }
            }
            //     ** Time evolution computations (parabolic systems)
            if (*ips == 14 || *ips == 16) {
                rap->tivp = rlold[0];
            }
            for (int64 i = 0; i < NPARX; ++i) {
                prm[i] = par[i];
            }
            (*funi)(iap, rap, *ndim, u, uold, icp, prm, 2, f, dfdu, dfdp);
            ic1 = ic**ndim;
            for (int64 i = 0; i < *ndim; ++i) {
                ARRAY2D(fa, ic1 + i, jj) = f[i] - ARRAY2D(wploc, *ncol, ic)*ARRAY2D(ups, jp1, i);
                for (int32 k = 0; k < *ncol; ++k) {
                    k1 = k**ndim + i;
                    ARRAY2D(fa, ic1 + i, jj) -= ARRAY2D(wploc, k, ic)*ARRAY2D(ups, j, k1);
                }
            }
            // L1:
        }
        // L2:
    }

    //     Generate FC :

    //     Boundary conditions :

    if (*nbc > 0) {
        for (int64 i = 0; i < *ndim; ++i) {
            ubc0[i] = ARRAY2D(ups, 0, i);
            ubc1[i] = ARRAY2D(ups, *ntst, i);
        }
        (*bcni)(iap, rap, *ndim, par, icp, *nbc, ubc0, ubc1, fbc, 2, dbc);
        for (int64 i = 0; i < *nbc; ++i) {
            fc[i] = -fbc[i];
        }
        //       Save difference :
        for (j = 0; j < *ntst + 1; ++j) {
            for (int64 i = 0; i < *nra; ++i) {
                ARRAY2D(dups, j, i) = ARRAY2D(ups, j, i) - ARRAY2D(uoldps, j, i);
            }
        }
    }

    //     Integral constraints :
    if (*nint > 0) {
        for (jj = 0; jj < *na; ++jj) {
            j = jj + mpart;
            jp1 = j + 1;
            for (int32 k = 0; k < ncp1; ++k) {
                for (int64 i = 0; i < *ndim; ++i) {
                    i1 = k**ndim + i;
                    j1 = j;
                    if (k + 1 == ncp1) {
                        i1 = i;
                    }
                    if (k + 1 == ncp1) {
                        j1 = jp1;
                    }
                    uic[i] = ARRAY2D(ups, j1, i1);
                    uio[i] = ARRAY2D(uoldps, j1, i1);
                    uid[i] = ARRAY2D(udotps, j1, i1);
                    uip[i] = ARRAY2D(upoldp, j1, i1);
                }
                (*icni)(iap, rap, *ndim, par, icp, *nint, uic, uio, uid, uip, ficd, 2, dicd);
                for (m = 0; m < *nint; ++m) {
                    fc[*nbc + m] -= dtm[j]*wi[k]*ficd[m];
                }
            }
        }
    }

    //     Pseudo-arclength equation :
    rlsum = 0.;
    for (int64 i = 0; i < *ncb; ++i) {
        rlsum += thl[icp[i]]*(rlcur[i] - rlold[i])*rldot[i];
    }

    fc[-1 + *nrc] = *rds - rinpr(iap, ndim, ndxloc, udotps, dups, dtm, thu) - rlsum;

    free(dicd);
    free(ficd);
    free(dfdp);
    free(dfdu);
    free(uold);
    free(f);
    free(u);
    free(wploc);
    free(wi);
    free(wp);
    free(wt);
    free(dbc);
    free(fbc);
    free(uic);
    free(uio);
    free(prm);
    free(uid);
    free(uip);
    free(ubc0);
    free(ubc1);

    return 0;
}

int32
brbd(double *a, double *b, double *c, double *d, double *fa, double *fc, double *p0, double *p1,
     int64 *ifst, int64 *idb, int64 *nllv, double *det, int64 *nov, int64 *na, int64 *nbc,
     int64 *nra, int64 *nca, int64 *ncb, int64 *nrc, int64 *iam, int64 *kwt, int64 *par, double *a1,
     double *a2, double *bb, double *cc, double *faa, double *ca1, double *s1, double *s2,
     int64 *icf11, int64 *ipr, int64 *icf1, int64 *icf2, int64 *irf, int64 *icf) {
    double *e;
    double *fcc;
    double *sol1, *sol2, *sol3;

    e = xmalloc(sizeof(*e)*(usize)((*nov + *nrc)*(*nov + *nrc)));
    fcc = xmalloc(sizeof(*fcc)*(usize)((*nov + *nrc) + (2*(*nov)*(*nov)) + 1));

    sol1 = xmalloc(sizeof(*sol1)*(usize)((*nov)*(*na + 1)));
    sol2 = xmalloc(sizeof(*sol2)*(usize)((*nov)*(*na + 1)));
    sol3 = xmalloc(sizeof(*sol3)*(usize)((*nov)*(*na + 1)));

    // Local

    if (*idb > 4 && *iam == 0) {
#ifndef ACCES_TEST
        print1(nov, na, nra, nca, ncb, nrc, a, b, c, d, &fa[0], fc);
#endif
    }
    if (*ifst == 1) {
#ifdef ACCES_TEST
        a = test.a;
        b = test.b;
        c = test.c;
#endif
        conpar(nov, na, nra, nca, a, ncb, b, nbc, nrc, c, d, irf, icf);
        copycp(iam, kwt, na, nov, nra, nca, a, ncb, b, nrc, c, a1, a2, bb, cc, irf);
    }

    if (*nllv == 0) {
        conrhs(nov, na, nra, nca, a, nbc, nrc, c, fa, fc, irf, icf, iam);
        cpyrhs(na, nov, nra, faa, fa, irf);
    } else {
#ifdef RANDY_FIX
        /* The faa array needs to be intialized as well, since it
           it used in the dimrge_ rountine to print stuff out,
           and in the bcksub_ routine for actual computations! */
        for (int64 k = 0; k < ((*nov)*(*na + 1)); k++) {
            faa[k] = 0.0;
        }
        setzero(fa, fc, na, nra, nrc);
#else
        setzero(fa, fc, na, nra, nrc);
        cpyrhs(na, nov, nra, faa, fa, irf);
#endif
    }

    if (*ifst == 1) {
        reduce(iam, kwt, par, a1, a2, bb, cc, d, na, nov, ncb, nrc, s1, s2, ca1, icf1, icf2, icf11,
               ipr, nbc);
    }

    if (*nllv == 0) {
        redrhs(iam, kwt, par, a1, a2, cc, faa, fc, na, nov, ncb, nrc, ca1, icf1, icf2, icf11, ipr,
               nbc);
    }

    dimrge(iam, kwt, par, e, cc, d, fc, ifst, na, nrc, nov, ncb, idb, nllv, fcc, p0, p1, det, s1,
           a2, faa, bb);

    bcksub(iam, kwt, par, s1, s2, a2, bb, faa, fc, fcc, sol1, sol2, sol3, na, nov, ncb, icf2);

    infpar(iam, par, a, b, fa, sol1, sol2, fc, na, nov, nra, nca, ncb, irf, icf);

    free(e);
    free(fcc);
    free(sol1);
    free(sol2);
    free(sol3);
    return 0;
}

int32
setzero(double *fa, double *fc, int64 *na, int64 *nra, int64 *nrc) {
    int64 fa_dim1 = *nra;

    for (int64 i = 0; i < *na; ++i) {
        for (int64 j = 0; j < *nra; ++j) {
            ARRAY2D(fa, j, i) = 0.;
        }
    }

    for (int64 i = 0; i < *nrc; ++i) {
        fc[i] = 0.;
    }

    return 0;
}

int32
conrhs(int64 *nov, int64 *na, int64 *nra, int64 *nca, double *a, int64 *nbc, int64 *nrc, double *c,
       double *fa, double *fc, int64 *irf, int64 *icf, int64 *iam) {
    int64 icf_dim1;
    int64 irf_dim1;
    int64 a_dim1;
    int64 a_dim2;
    int64 c_dim1;
    int64 c_dim2;
    int64 fa_dim1;

    int64 nbcp1;
    int64 icfic;
    int64 irfir;
    int64 m1;
    int64 m2;
    int64 ic;
    int64 ir;
    int64 irfirp;
    int64 ir1;
    int64 nex;
    int64 irp;

    (void)iam;

    irf_dim1 = *nra;
    fa_dim1 = *nra;
    icf_dim1 = *nca;
    a_dim1 = *nca;
    a_dim2 = *nra;
    c_dim1 = *nca;
    c_dim2 = *nrc;

    nex = *nca - (*nov*2);
    if (nex == 0) {
        return 0;
    }

    // Condensation of right hand side.

    nbcp1 = *nbc + 1;
    m1 = *nov + 1;
    m2 = *nov + nex;

    for (int64 i = 0; i < *na; ++i) {
        for (ic = *nov; ic < m2; ++ic) {
            ir1 = ic - *nov + 1;
            irp = ir1 - 1;
            irfirp = ARRAY2D(irf, irp, i);
            icfic = ARRAY2D(icf, ic, i);
            for (ir = ir1; ir < *nra; ++ir) {
                irfir = ARRAY2D(irf, ir, i);
                if (ARRAY3D(a, (icfic - 1), (irfir - 1), i) != (double)0.) {
                    ARRAY2D(fa, (irfir - 1), i) -=
                        ARRAY3D(a, (icfic - 1), (irfir - 1), i)*ARRAY2D(fa, (irfirp - 1), i);
                }
            }
            for (ir = *nbc; ir < *nrc; ++ir) {
                if (ARRAY3D(c, (icfic - 1), ir, i) != (double)0.) {
                    fc[ir] -= ARRAY3D(c, (icfic - 1), ir, i)*ARRAY2D(fa, (irfirp - 1), i);
                }
            }
        }
    }

    return 0;
}

int32
copycp(int64 *iam, int64 *kwt, int64 *na, int64 *nov, int64 *nra, int64 *nca, double *a, int64 *ncb,
       double *b, int64 *nrc, double *c, double *a1, double *a2, double *bb, double *cc,
       int64 *irf) {
    int64 irf_dim1, a_dim1, a_dim2, b_dim1, b_dim2, c_dim1, c_dim2, a1_dim1, a1_dim2, a2_dim1,
        a2_dim2, bb_dim1, bb_dim2, cc_dim1, cc_dim2;

    int64 irfir;
    int64 ic;
    int64 ir;
    int64 ic1;
    int64 nap1;

    (void)iam;
    (void)kwt;

    // Local

    // Copies the condensed sytem generated by CONPAR into workspace.

    a2_dim1 = *nov;
    a2_dim2 = *nov;
    a1_dim1 = *nov;
    a1_dim2 = *nov;
    irf_dim1 = *nra;
    a_dim1 = *nca;
    a_dim2 = *nra;
    bb_dim1 = *nov;
    bb_dim2 = *ncb;
    b_dim1 = *ncb;
    b_dim2 = *nra;
    cc_dim1 = *nov;
    cc_dim2 = *nrc;
    c_dim1 = *nca;
    c_dim2 = *nrc;

    nap1 = *na + 1;
    for (int64 i = 0; i < *na; ++i) {
        for (ir = 0; ir < *nov; ++ir) {
            irfir = ARRAY2D(irf, *nra - *nov + ir, i);
            for (ic = 0; ic < *nov; ++ic) {
                ic1 = *nca - *nov + ic;
                ARRAY3D(a1, ir, ic, i) = ARRAY3D(a, ic, (irfir - 1), i);
                ARRAY3D(a2, ir, ic, i) = ARRAY3D(a, ic1, (irfir - 1), i);
            }
            for (ic = 0; ic < *ncb; ++ic) {
                ARRAY3D(bb, ir, ic, i) = ARRAY3D(b, ic, (irfir - 1), i);
            }
        }
    }

    for (int64 i = 0; i < nap1; ++i) {
        for (ir = 0; ir < *nrc; ++ir) {
            for (ic = 0; ic < *nov; ++ic) {
                if (i + 1 == 1) {
                    ARRAY3D(cc, ic, ir, i) = ARRAY3D(c, ic, ir, i);
                } else if (i + 1 == nap1) {
                    ARRAY3D(cc, ic, ir, i) = ARRAY3D(c, *nra + ic, ir, (i - 1));
                } else {
                    ARRAY3D(cc, ic, ir, i) =
                        ARRAY3D(c, ic, ir, i) + ARRAY3D(c, *nra + ic, ir, (i - 1));
                }
            }
        }
    }

    return 0;
}

int32
cpyrhs(int64 *na, int64 *nov, int64 *nra, double *faa, double *fa, int64 *irf) {
    int64 irf_dim1;
    int64 fa_dim1;
    int64 faa_dim1;

    int64 irfir;
    int64 ir;

    faa_dim1 = *nov;
    irf_dim1 = *nra;
    fa_dim1 = *nra;

    for (int64 i = 0; i < *na; ++i) {
        for (ir = 0; ir < *nov; ++ir) {
            irfir = ARRAY2D(irf, *nra - *nov + ir, i);
            ARRAY2D(faa, ir, i) = ARRAY2D(fa, (irfir - 1), i);
        }
    }

    return 0;
}

int32
reduce(int64 *iam, int64 *kwt, int64 *par, double *a1, double *a2, double *bb, double *cc,
       double *dd, int64 *na, int64 *nov, int64 *ncb, int64 *nrc, double *s1, double *s2,
       double *ca1, int64 *icf1, int64 *icf2, int64 *icf11, int64 *ipr, int64 *nbc) {
    int64 icf1_dim1, icf2_dim1, icf11_dim1, a1_dim1, a1_dim2, a2_dim1, a2_dim2, s1_dim1, s1_dim2,
        s2_dim1, s2_dim2, bb_dim1, bb_dim2, cc_dim1, cc_dim2, dd_dim1, ca1_dim1, ca1_dim2, ipr_dim1;

    int64 oddc[KREDO];
    int64 niam, ibuf, ismc[KREDO], irmc[KREDO], info, irmm[KREDO], ismm[KREDO], nlev, itmp;
    double zero;
    double tpiv;
    double xkwt;
    int64 nbcp1;
    int64 ibuf1;
    int64 ipiv1;
    int64 jpiv1;
    int64 ipiv2;
    int64 jpiv2;
    int64 k;
    int64 l;

    int64 evenc[KREDO];

    int64 i1;
    int64 i2;
    int64 k1;
    int64 k2;
    int64 i3;
    int64 l1;
    int64 iprow;
    int64 k3;
    int64 l2;
    int64 l3;
    int64 ic;
    int64 ir;
    double rm;
    int64 master[KREDO];
    int64 ib1;
    int64 ib2;
    int64 myleft[KREDO];

    int64 worker[KREDO];
    int64 ir1;
    int64 iprown;
    int64 iprown2;
    int64 ism[KREDO];
    int64 irm[KREDO];
    int64 nrcmnbc;
    double tmp;
    int64 myleftc[KREDO];
    int64 notsend;
    int64 nap1;
    int64 myright[KREDO];
    int64 nam1;
    int64 len1;
    int64 len2;
    int64 icp1;
    double piv1;
    double piv2;
    double *buf = NULL;

    ipr_dim1 = *nov;
    icf11_dim1 = *nov;
    icf2_dim1 = *nov;
    icf1_dim1 = *nov;
    ca1_dim1 = *nov;
    ca1_dim2 = *nov;
    s2_dim1 = *nov;
    s2_dim2 = *nov;
    s1_dim1 = *nov;
    s1_dim2 = *nov;
    a2_dim1 = *nov;
    a2_dim2 = *nov;
    a1_dim1 = *nov;
    a1_dim2 = *nov;
    dd_dim1 = *ncb;
    bb_dim1 = *nov;
    bb_dim2 = *ncb;
    cc_dim1 = *nov;
    cc_dim2 = *nrc;

    zero = 0.;
    nbcp1 = *nbc + 1;
    nap1 = *na + 1;
    nam1 = *na - 1;
    nrcmnbc = *nrc - *nbc;
    len1 = (*nov*(*nrc - *nbc))*8;
    len2 = (*nov + *nrc - *nbc + 1)*8;
    xkwt = (double)(*kwt);
    {
        double tmp2 = r_lg10(xkwt) / r_lg10(2.0);
        nlev = i_nint(&tmp2);
    }
    notsend = true;

    //     FOR EACH REURSIVE LEVEL, CALCULATE THE MASTER(HOLDING THE
    //     PIVOT ROW AFTER ROW SWAPPING) NODE WHICH WILL SEND THE
    //     PIVOT ROW TO THE CORRESPONDING WORKER NODE WHICH IS DISTANCED
    //     2**(K-1) FROM THE MASTER WHERE K IS THE RECURSIVE LEVEL NUMBER.
    //     THE CORRESPONDING MESSAGE TYPE IN EACH RECURSIVE LEVEL IS
    //     ALSO CALCULATED HERE.

    // For each level in the recursion, determine the master node
    // (holding the pivot row after row swapping), which will send the
    // pivot row to the corresponding worker node at distance 2**(K-1)
    // from the master. Here K is the level in the recursion.
    // The message type at each level in the recursion is also determined.

    if (*par) {
        for (int64 i = 0; i < nlev; ++i) {
            oddc[i] = false;
            evenc[i] = false;
            master[i] = false;
            worker[i] = false;

            k1 = pow_ii(2, i);
            k2 = k1*2;
            niam = *iam / k1;

            if (notsend) {
                if (niam % 2 == 0) {
                    master[i] = true;
                    notsend = false;
                    ism[i] = (i + 1) + *iam;
                    irm[i] = ism[i] + k1;
                    myright[i] = *iam + k1;
                    irmm[i] = (i + 1) + *iam + 1 + (*kwt*2);
                    ismc[i] = (i + 1) + *iam + *kwt;
                    myleftc[i] = *iam - (k1 - 1);

                } else {

                    worker[i] = true;
                    ism[i] = (i + 1) + *iam;
                    irm[i] = ism[i] - k1;
                    myleft[i] = *iam - k1;
                }
            }

            k = *iam % k2;
            if (k == k1) {
                evenc[i] = true;
                ismm[i] = (i + 1) + *iam + (*kwt*2);
            }

            if (*iam % k2 == 0) {
                oddc[i] = true;
                irmc[i] = (i + 1) + *iam + *kwt + (k1 - 1);
            }

            // L1:
        }
    }

    // Initialization

    for (int64 i = 0; i < *na; ++i) {
        for (k1 = 0; k1 < *nov; ++k1) {
            ARRAY2D(icf1, k1, i) = k1 + 1;
            ARRAY2D(icf2, k1, i) = k1 + 1;
            ARRAY2D(ipr, k1, i) = k1 + 1;
            for (k2 = 0; k2 < *nov; ++k2) {
                ARRAY3D(s2, k1, k2, i) = 0.;
                ARRAY3D(s1, k1, k2, i) = 0.;
            }
        }
    }

    for (ir = 0; ir < *nov; ++ir) {
        for (ic = 0; ic < *nov; ++ic) {
            ARRAY2D(s1, ir, ic) = ARRAY2D(a1, ir, ic);
        }
    }

    // The reduction process is done concurrently
    for (i1 = 0; i1 < nam1; ++i1) {
        i2 = i1 + 1;
        i3 = i2 + 1;

        for (ic = 0; ic < *nov; ++ic) {
            icp1 = ic + 1;

            /* Complete pivoting; rows are swapped physically, columns swap in
               dices */
            piv1 = zero;
            ipiv1 = ic + 1;
            jpiv1 = ic + 1;
            for (k1 = ic; k1 < *nov; ++k1) {
                for (k2 = ic; k2 < *nov; ++k2) {
                    tpiv = ARRAY3D(a2, k1, ARRAY2D(icf2, k2, i1) - 1, i1);
                    if (tpiv < zero) {
                        tpiv = -tpiv;
                    }
                    if (piv1 < tpiv) {
                        piv1 = tpiv;
                        ipiv1 = k1 + 1;
                        jpiv1 = k2 + 1;
                    }
                }
            }

            piv2 = zero;
            ipiv2 = 1;
            jpiv2 = ic + 1;
            for (k1 = 0; k1 < *nov; ++k1) {
                for (k2 = ic; k2 < *nov; ++k2) {
                    tpiv = ARRAY3D(a1, k1, ARRAY2D(icf1, k2, i2) - 1, i2);
                    if (tpiv < zero) {
                        tpiv = -tpiv;
                    }
                    if (piv2 < tpiv) {
                        piv2 = tpiv;
                        ipiv2 = k1 + 1;
                        jpiv2 = k2 + 1;
                    }
                }
            }

            if (piv1 >= piv2) {
                ARRAY2D(ipr, ic, i1) = ipiv1;
                itmp = ARRAY2D(icf2, ic, i1);
                ARRAY2D(icf2, ic, i1) = ARRAY2D(icf2, (jpiv1 - 1), i1);
                ARRAY2D(icf2, (jpiv1 - 1), i1) = itmp;
                itmp = ARRAY2D(icf1, ic, i2);
                ARRAY2D(icf1, ic, i2) = ARRAY2D(icf1, (jpiv1 - 1), i2);
                ARRAY2D(icf1, (jpiv1 - 1), i2) = itmp;
                // Swapping
                for (l = 0; l < *nov; ++l) {
                    tmp = ARRAY3D(s1, ic, l, i1);
                    ARRAY3D(s1, ic, l, i1) = ARRAY3D(s1, (ipiv1 - 1), l, i1);
                    ARRAY3D(s1, (ipiv1 - 1), l, i1) = tmp;
                    if (l >= ic) {
                        tmp = ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1);
                        ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1) =
                            ARRAY3D(a2, (ipiv1 - 1), ARRAY2D(icf2, l, i1) - 1, i1);
                        ARRAY3D(a2, (ipiv1 - 1), ARRAY2D(icf2, l, i1) - 1, i1) = tmp;
                    }
                    tmp = ARRAY3D(s2, ic, l, i1);
                    ARRAY3D(s2, ic, l, i1) = ARRAY3D(s2, (ipiv1 - 1), l, i1);
                    ARRAY3D(s2, (ipiv1 - 1), l, i1) = tmp;
                }

                for (l = 0; l < *ncb; ++l) {
                    tmp = ARRAY3D(bb, ic, l, i1);
                    ARRAY3D(bb, ic, l, i1) = ARRAY3D(bb, (ipiv1 - 1), l, i1);
                    ARRAY3D(bb, (ipiv1 - 1), l, i1) = tmp;
                }
            } else {
                ARRAY2D(ipr, ic, i1) = *nov + ipiv2;
                itmp = ARRAY2D(icf2, ic, i1);
                ARRAY2D(icf2, ic, i1) = ARRAY2D(icf2, (jpiv2 - 1), i1);
                ARRAY2D(icf2, (jpiv2 - 1), i1) = itmp;
                itmp = ARRAY2D(icf1, ic, i2);
                ARRAY2D(icf1, ic, i2) = ARRAY2D(icf1, (jpiv2 - 1), i2);
                ARRAY2D(icf1, (jpiv2 - 1), i2) = itmp;
                // Swapping
                for (l = 0; l < *nov; ++l) {
                    if (l >= ic) {
                        tmp = ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1);
                        ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1) =
                            ARRAY3D(a1, (ipiv2 - 1), ARRAY2D(icf2, l, i1) - 1, i2);
                        ARRAY3D(a1, (ipiv2 - 1), ARRAY2D(icf2, l, i1) - 1, i2) = tmp;
                    }
                    tmp = ARRAY3D(s2, ic, l, i1);
                    ARRAY3D(s2, ic, l, i1) = ARRAY3D(a2, (ipiv2 - 1), l, i2);
                    ARRAY3D(a2, (ipiv2 - 1), l, i2) = tmp;
                    tmp = ARRAY3D(s1, ic, l, i1);
                    ARRAY3D(s1, ic, l, i1) = ARRAY3D(s1, (ipiv2 - 1), l, i2);
                    ARRAY3D(s1, (ipiv2 - 1), l, i2) = tmp;
                }
                for (l = 0; l < *ncb; ++l) {
                    tmp = ARRAY3D(bb, ic, l, i1);
                    ARRAY3D(bb, ic, l, i1) = ARRAY3D(bb, (ipiv2 - 1), l, i2);
                    ARRAY3D(bb, (ipiv2 - 1), l, i2) = tmp;
                }
            }

            // End of pivoting; Elimination starts here

            for (ir = icp1; ir < *nov; ++ir) {
                rm = ARRAY3D(a2, ir, ARRAY2D(icf2, ic, i1) - 1, i1) /
                     ARRAY3D(a2, ic, ARRAY2D(icf2, ic, i1) - 1, i1);
                ARRAY3D(a2, ir, ARRAY2D(icf2, ic, i1) - 1, i1) = rm;

                if (rm != (double)0.) {
                    for (l = icp1; l < *nov; ++l) {
                        ARRAY3D(a2, ir, ARRAY2D(icf2, l, i1) - 1, i1) -=
                            rm*ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1);
                    }

                    for (l = 0; l < *nov; ++l) {
                        ARRAY3D(s1, ir, l, i1) -= rm*ARRAY3D(s1, ic, l, i1);
                        ARRAY3D(s2, ir, l, i1) -= rm*ARRAY3D(s2, ic, l, i1);
                    }

                    for (l = 0; l < *ncb; ++l) {
                        ARRAY3D(bb, ir, l, i1) -= rm*ARRAY3D(bb, ic, l, i1);
                    }
                }
            }

            for (ir = 0; ir < *nov; ++ir) {
                rm = ARRAY3D(a1, ir, ARRAY2D(icf1, ic, i2) - 1, i2) /
                     ARRAY3D(a2, ic, ARRAY2D(icf2, ic, i1) - 1, i1);
                ARRAY3D(a1, ir, ARRAY2D(icf1, ic, i2) - 1, i2) = rm;

                if (rm != (double)0.) {
                    for (l = icp1; l < *nov; ++l) {
                        ARRAY3D(a1, ir, ARRAY2D(icf1, l, i2) - 1, i2) -=
                            rm*ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1);
                    }
                    for (l = 0; l < *nov; ++l) {
                        ARRAY3D(s1, ir, l, i2) -= rm*ARRAY3D(s1, ic, l, i1);
                        ARRAY3D(a2, ir, l, i2) -= rm*ARRAY3D(s2, ic, l, i1);
                    }
                    for (l = 0; l < *ncb; ++l) {
                        ARRAY3D(bb, ir, l, i2) -= rm*ARRAY3D(bb, ic, l, i1);
                    }
                }
            }

            for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
                rm = ARRAY3D(cc, ARRAY2D(icf2, ic, i1) - 1, ir, i2) /
                     ARRAY3D(a2, ic, ARRAY2D(icf2, ic, i1) - 1, i1);
                ARRAY3D(cc, ARRAY2D(icf2, ic, i1) - 1, ir, i2) = rm;

                if (rm != (double)0.) {
                    for (l = icp1; l < *nov; ++l) {
                        ARRAY3D(cc, ARRAY2D(icf2, l, i1) - 1, ir, i2) -=
                            rm*ARRAY3D(a2, ic, ARRAY2D(icf2, l, i1) - 1, i1);
                    }
                    for (l = 0; l < *nov; ++l) {
                        ARRAY3D(cc, l, ir, 0) -= rm*ARRAY3D(s1, ic, l, i1);
                        ARRAY3D(cc, l, ir, i3) -= rm*ARRAY3D(s2, ic, l, i1);
                    }
                    for (l = 0; l < *ncb; ++l) {
                        ARRAY2D(dd, l, ir) -= rm*ARRAY3D(bb, ic, l, i1);
                    }
                }
            }

            // L2:
        }
        // L3:
    }

    // Initialization
    for (int64 i = 0; i < *nov; ++i) {
        ARRAY2D(icf2, i, (*na - 1)) = i + 1;
    }

    //     INTER NODES REDUCE IS DONE VIA COMMUNICATION
    //     BETWEEN MASTER NODES AND WORKER NODES.
    //     THE SUMMATION OF THE OVERLAPED PART C IN THE
    //     NEIGHBOR NODES IN THE CONDENSATION OF PARAMETE
    //     ROUTINE IS DELAYED TO HERE TO SUM.

    // Inter node reduction is done via communication between master node
    // and worker nodes. The summation over the overlapped part C of
    /*neighboring nodes in the condensation of parameters is delayed until her
    e.*/
    if (*par) {
        for (int64 i = 0; i < nlev; ++i) {
            if (master[i]) {
                crecv();
                for (ir = nbcp1; ir <= *nrc; ++ir) {
                    ir1 = ir - *nbc;
                    for (ic = 1; ic <= *nov; ++ic) {
                        l1 = ir1**nov + ic;
                        ARRAY3D(cc, ic, ir, (nap1 - 1)) += buf[l1 + 1];
                    }
                }
                for (ir = 0; ir < *nov; ++ir) {
                    for (ic = 0; ic < *nov; ++ic) {
                        ARRAY3D(s2, ir, ic, *na) = 0.;
                    }
                }
            }

            if (evenc[i]) {
                csend();
            }

            if (worker[i]) {
                for (ir = 0; ir < *nov; ++ir) {
                    for (ic = 0; ic < *nov; ++ic) {
                        ARRAY3D(ca1, ir, ic, i) = ARRAY3D(s1, ir, ic, (*na - 1));
                        ARRAY3D(s1, ir, ic, (*na - 1)) = 0.;
                    }
                }

                for (l = 0; l < *nov; ++l) {
                    ARRAY2D(icf11, l, i) = l + 1;
                }
            }

            for (ic = 0; ic < *nov; ++ic) {
                icp1 = ic + 1;
                iprow = *nov - ic + 1;
                iprown = iprow + *nov;
                iprown2 = iprown + *nov;
                ib1 = iprown2 + *ncb + 1;
                ib2 = ib1 + 1;
                ibuf = (ib2 + 1)*8;
                ibuf1 = (ib2 + *nrc - *nbc)*8;

                if (master[i]) {
                    // PIVOTING (COMPLETE PIVOTING)

                    piv1 = zero;
                    ipiv1 = ic + 1;
                    jpiv1 = ic + 1;
                    for (k1 = ic; k1 < *nov; ++k1) {
                        for (k2 = ic; k2 < *nov; ++k2) {
                            k3 = ARRAY2D(icf2, k2, (*na - 1));
                            tpiv = ARRAY3D(a2, k1, k3, (*na - 1));
                            if (tpiv < zero) {
                                tpiv = -tpiv;
                            }
                            if (piv1 < tpiv) {
                                piv1 = tpiv;
                                ipiv1 = k1 + 1;
                                jpiv1 = k2 + 1;
                            }
                        }
                    }

                    crecv();

                    jpiv2 = i_dnnt(&buf[ib1 + 1]);
                    ipiv2 = i_dnnt(&buf[ib2 + 1]);

                    piv2 = buf[1];
                    if (piv2 < 0.) {
                        piv2 = -piv2;
                    }

                    if (piv1 >= piv2) {
                        ARRAY2D(ipr, ic, (*na - 1)) = ipiv1;
                        itmp = ARRAY2D(icf2, ic, (*na - 1));
                        ARRAY2D(icf2, ic, (*na - 1)) = ARRAY2D(icf2, (jpiv1 - 1), (*na - 1));
                        ARRAY2D(icf2, (jpiv1 - 1), (*na - 1)) = itmp;

                        // Send pivot row to worker
                        for (l = 0; l < *nov; ++l) {
                            if (l >= ic) {
                                l1 = l - ic + 2;
                                l2 = ARRAY2D(icf2, l, (*na - 1)) - 1;
                                buf[l1 + 1] = ARRAY3D(a2, (ipiv1 - 1), l2, (*na - 1));
                            }
                            l1 = iprow + l;
                            l2 = iprown + l;
                            buf[l1 + 1] = ARRAY3D(s1, (ipiv1 - 1), l, (*na - 1));
                            buf[l2 + 1] = ARRAY3D(s2, (ipiv1 - 1), l, (*na - 1));
                        }

                        for (l = 0; l < *nbc; ++l) {
                            l1 = iprown2 + l;
                            buf[l1 + 1] = ARRAY3D(bb, (ipiv1 - 1), l, (*na - 1));
                        }

                        buf[ib1 + 1] = (double)jpiv1;

                        for (l = nbcp1 - 1; l < *nrc; ++l) {
                            l1 = l - *nbc;
                            l2 = ib1 + l1;
                            l3 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
                            buf[l2 + 1] = ARRAY3D(cc, l3, l, (nap1 - 1));
                        }

                        l1 = ib2 + nrcmnbc;
                        buf[l1 + 1] = 0.;
                        csend();

                        // Row swapping
                        for (l = 0; l < *nov; ++l) {
                            tmp = ARRAY3D(s1, ic, l, (*na - 1));
                            ARRAY3D(s1, ic, l, (*na - 1)) = ARRAY3D(s1, (ipiv1 - 1), l, (*na - 1));
                            ARRAY3D(s1, (ipiv1 - 1), l, (*na - 1)) = tmp;
                            if (l >= ic) {
                                l1 = ARRAY2D(icf2, l, (*na - 1)) - 1;
                                tmp = ARRAY3D(a2, ic, l1, (*na - 1));
                                ARRAY3D(a2, ic, l1, (*na - 1)) =
                                    ARRAY3D(a2, (ipiv1 - 1), l1, (*na - 1));
                                ARRAY3D(a2, (ipiv1 - 1), l1, (*na - 1)) = tmp;
                            }
                            tmp = ARRAY3D(s2, ic, l, (*na - 1));
                            ARRAY3D(s2, ic, l, (*na - 1)) = ARRAY3D(s2, (ipiv1 - 1), l, (*na - 1));
                            ARRAY3D(s2, (ipiv1 - 1), l, (*na - 1)) = tmp;
                        }

                        for (l = 0; l < *ncb; ++l) {
                            tmp = ARRAY3D(bb, ic, l, (*na - 1));
                            ARRAY3D(bb, ic, l, (*na - 1)) = ARRAY3D(bb, (ipiv1 - 1), l, (*na - 1));
                            ARRAY3D(bb, (ipiv1 - 1), l, (*na - 1)) = tmp;
                        }

                    } else {

                        ARRAY2D(ipr, ic, (*na - 1)) = *nov + ipiv2;
                        jpiv1 = jpiv2;
                        itmp = ARRAY2D(icf2, ic, (*na - 1));
                        ARRAY2D(icf2, ic, (*na - 1)) = ARRAY2D(icf2, (jpiv1 - 1), (*na - 1));
                        ARRAY2D(icf2, (jpiv1 - 1), (*na - 1)) = itmp;

                        for (l = 0; l < *nov; ++l) {
                            if (l >= ic) {
                                l1 = l - ic + 2;
                                l2 = ARRAY2D(icf2, l, (*na - 1)) - 1;
                                tmp = buf[l1 + 1];
                                buf[l1 + 1] = ARRAY3D(a2, ic, l2, (*na - 1));
                                ARRAY3D(a2, ic, l2, (*na - 1)) = tmp;
                            }
                            l1 = iprow + l;
                            l2 = iprown + l;
                            tmp = buf[l1 + 1];
                            buf[l1 + 1] = ARRAY3D(s1, ic, l, (*na - 1));
                            ARRAY3D(s1, ic, l, (*na - 1)) = tmp;
                            tmp = buf[l2 + 1];
                            buf[l2 + 1] = ARRAY3D(s2, ic, l, (*na - 1));
                            ARRAY3D(s2, ic, l, (*na - 1)) = tmp;
                        }

                        for (l = 0; l < *nbc; ++l) {
                            l1 = iprown2 + l;
                            tmp = buf[l1 + 1];
                            buf[l1 + 1] = ARRAY3D(bb, ic, l, (*na - 1));
                            ARRAY3D(bb, ic, l, (*na - 1)) = tmp;
                        }

                        buf[ib1 + 1] = (double)jpiv2;

                        for (l = nbcp1; l <= *nrc; ++l) {
                            l1 = l - *nbc;
                            l2 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
                            l3 = ib1 + l1 + 1;
                            buf[l3 + 1] = ARRAY3D(cc, l2, l, (nap1 - 1));
                        }
                        l1 = ib2 + nrcmnbc;
                        buf[l1 + 1] = 1.;

                        csend();
                    }
                    // End pivoting in master

                    // Send data to worker nodes
                    for (l = 0; l < *nov; ++l) {
                        buf[l + 1] = ARRAY3D(s1, ic, l, (*na - 1));
                    }

                    for (l = nbcp1 - 1; l < *nrc; ++l) {
                        l1 = l - *nbc;
                        l2 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
                        l3 = *nov + l1;
                        buf[l3 + 1] = ARRAY3D(cc, l2, l, (nap1 - 1));
                    }

                    l2 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
                    l1 = *nov + nrcmnbc;
                    buf[l1 + 1] = ARRAY3D(a2, ic, l2, (*na - 1));

                    csend();

                    // Elimination
                    for (ir = icp1; ir < *nov; ++ir) {
                        l2 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
                        rm = ARRAY3D(a2, ir, l2, (*na - 1)) / ARRAY3D(a2, ic, l2, (*na - 1));
                        ARRAY3D(a2, ir, l2, (*na - 1)) = rm;
                        if (rm != zero) {
                            for (l = icp1; l < *nov; ++l) {
                                l1 = ARRAY2D(icf2, l, (*na - 1)) - 1;
                                ARRAY3D(a2, ir, l1, (*na - 1)) -=
                                    rm*ARRAY3D(a2, ic, l1, (*na - 1));
                            }
                            for (l = 0; l < *nov; ++l) {
                                ARRAY3D(s1, ir, l, (*na - 1)) -= rm*ARRAY3D(s1, ic, l, (*na - 1));
                                ARRAY3D(s2, ir, l, (*na - 1)) -= rm*ARRAY3D(s2, ic, l, (*na - 1));
                            }
                            for (l = 0; l < *ncb; ++l) {
                                ARRAY3D(bb, ir, l, (*na - 1)) -= rm*ARRAY3D(bb, ic, l, (*na - 1));
                            }
                        }
                    }

                    for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
                        l2 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
                        rm = ARRAY3D(cc, l2, ir, (nap1 - 1)) / ARRAY3D(a2, ic, l2, (*na - 1));
                        ARRAY3D(cc, l2, ir, (nap1 - 1)) = rm;
                        if (rm != zero) {
                            for (l = icp1; l <= *nov; ++l) {
                                l1 = ARRAY2D(icf2, l, (*na - 1)) - 1;
                                ARRAY3D(cc, l1, ir, (nap1 - 1)) -=
                                    rm*ARRAY3D(a2, ic, l1, (*na - 1));
                            }
                            for (l = 0; l < *nbc; ++l) {
                                ARRAY2D(dd, l, ir) -= rm*ARRAY3D(bb, ic, l, (*na - 1));
                            }
                        }
                    }
                }

                if (worker[i]) {
                    // Pivoting
                    piv2 = zero;
                    ipiv2 = 1;
                    jpiv2 = ic + 1;
                    for (k1 = 0; k1 < *nov; ++k1) {
                        for (k2 = ic; k2 < *nov; ++k2) {
                            k3 = ARRAY2D(icf11, k2, i) - 1;
                            tpiv = ARRAY3D(ca1, k1, k3, i);
                            if (tpiv < zero) {
                                tpiv = -tpiv;
                            }
                            if (piv2 < tpiv) {
                                piv2 = tpiv;
                                ipiv2 = k1 + 1;
                                jpiv2 = k2 + 1;
                            }
                        }
                    }

                    itmp = ARRAY2D(icf11, ic, i);
                    ARRAY2D(icf11, ic, i) = ARRAY2D(icf11, (jpiv2 - 1), i);
                    ARRAY2D(icf11, (jpiv2 - 1), i) = itmp;

                    for (l = 0; l < *nov; ++l) {
                        if (l >= ic) {
                            l1 = l - ic + 2;
                            l2 = ARRAY2D(icf11, l, i) - 1;
                            buf[l1 + 1] = ARRAY3D(ca1, (ipiv2 - 1), l2, (*na - 1));
                        }
                        l1 = iprow + l;
                        l2 = l1 + *nov;
                        buf[l1 + 1] = ARRAY3D(s1, (ipiv2 - 1), l, i);
                        buf[l2 + 1] = ARRAY3D(a2, (ipiv2 - 1), l, (*na - 1));
                    }

                    for (l = 0; l < *ncb; ++l) {
                        l1 = iprown2 + l;
                        buf[l1 + 1] = ARRAY3D(bb, (ipiv2 - 1), l, (*na - 1));
                    }

                    buf[ib1 + 1] = (double)jpiv2;
                    buf[ib2 + 1] = (double)ipiv2;

                    csend();
                    crecv();

                    l1 = ib2 + nrcmnbc;
                    info = i_dnnt(&buf[l1 + 1]);

                    if (info == 1) {
                        // Send pivot row to master
                        for (l = 0; l < *nov; ++l) {
                            if (l >= ic) {
                                l1 = l - ic + 2;
                                l2 = ARRAY2D(icf11, l, i) - 1;
                                tmp = ARRAY3D(ca1, (ipiv2 - 1), l2, i);
                                ARRAY3D(ca1, (ipiv2 - 1), l2, i) = buf[l1 + 1];
                                buf[l1 + 1] = tmp;
                            }
                            l1 = iprow + l;
                            l2 = l1 + *nov;
                            tmp = ARRAY3D(s1, (ipiv2 - 1), l, (*na - 1));
                            ARRAY3D(s1, (ipiv2 - 1), l, (*na - 1)) = buf[l1 + 1];
                            buf[l1 + 1] = tmp;
                            tmp = ARRAY3D(a2, (ipiv2 - 1), l, (*na - 1));
                            ARRAY3D(a2, (ipiv2 - 1), l, (*na - 1)) = buf[l2 + 1];
                            buf[l2 + 1] = tmp;
                        }
                        for (l = 0; l < *nbc; ++l) {
                            l1 = iprown2 + l;
                            tmp = ARRAY3D(bb, (ipiv2 - 1), l, (*na - 1));
                            ARRAY3D(bb, (ipiv2 - 1), l, (*na - 1)) = buf[l1 + 1];
                            buf[l1 + 1] = tmp;
                        }
                    } else {

                        itmp = ARRAY2D(icf11, ic, i);
                        ARRAY2D(icf11, ic, i) = ARRAY2D(icf11, (jpiv2 - 1), i);
                        ARRAY2D(icf11, (jpiv2 - 1), i) = itmp;

                        jpiv2 = i_dnnt(&buf[ib1 + 1]);
                        itmp = ARRAY2D(icf11, ic, i);
                        ARRAY2D(icf11, ic, i) = ARRAY2D(icf11, (jpiv2 - 1), i);
                        ARRAY2D(icf11, (jpiv2 - 1), i) = itmp;
                    }

                    // Elimination
                    for (ir = 1; ir <= *nov; ++ir) {
                        l2 = ARRAY2D(icf11, ic, i) - 1;
                        rm = ARRAY3D(ca1, ir, l2, i) / buf[1];
                        ARRAY3D(ca1, ir, l2, i) = rm;

                        if (rm != zero) {
                            for (l = icp1; l < *nov; ++l) {
                                l1 = l - icp1 + 3;
                                l3 = ARRAY2D(icf11, l, i) - 1;
                                ARRAY3D(ca1, ir, l3, i) -= rm*buf[l1 + 1];
                            }
                            for (l = 0; l < *nov; ++l) {
                                l1 = iprow + l;
                                l2 = l1 + *nov;
                                ARRAY3D(s1, ir, l, (*na - 1)) -= rm*buf[l1 + 1];
                                ARRAY3D(a2, ir, l, (*na - 1)) -= rm*buf[l2 + 1];
                            }
                            for (l = 0; l < *ncb; ++l) {
                                l1 = iprown2 + l;
                                ARRAY3D(bb, ir, l, (*na - 1)) -= rm*buf[l1 + 1];
                            }
                        }
                    }

                    for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
                        l1 = ir - *nbc;
                        l2 = ib1 + l1;
                        rm = buf[l2 + 1] / buf[1];
                        if (rm != zero) {
                            for (l = 0; l < *nov; ++l) {
                                l3 = iprown + l;
                                ARRAY3D(cc, l, ir, (nap1 - 1)) -= rm*buf[l3 + 1];
                            }
                        }
                    }
                }

                if (oddc[i]) {
                    crecv();
                    for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
                        ir1 = ir - *nbc;
                        l1 = *nov + nrcmnbc;
                        l2 = *nov + ir1;
                        rm = buf[l2 + 1] / buf[l1 + 1];
                        if (rm != zero) {
                            for (l = 0; l < *nov; ++l) {
                                ARRAY3D(cc, l, ir, 0) -= rm*buf[l + 1];
                            }
                        }
                    }
                }
            }

            // L4:
        }

        // Global sum for D by recursive doubling
        {
            int64 tmp2 = (*nrc - *nbc)**ncb;
            rd0(iam, kwt, &ARRAY2D(dd, 0, (nbcp1 - 1)), &tmp2);
        }
    }

    return 0;
}

int32
redrhs(int64 *iam, int64 *kwt, int64 *par, double *a1, double *a2, double *cc, double *faa,
       double *fc, int64 *na, int64 *nov, int64 *ncb, int64 *nrc, double *ca1, int64 *icf1,
       int64 *icf2, int64 *icf11, int64 *ipr, int64 *nbc) {
    int64 icf1_dim1, icf2_dim1, icf11_dim1, a1_dim1, a1_dim2, a2_dim1, a2_dim2, cc_dim1, cc_dim2,
        faa_dim1, ca1_dim1, ca1_dim2, ipr_dim1;

    int64 niam;
    int64 nlev;
    double xkwt;
    int64 nbcp1;
    int64 ipiv1;
    int64 ipiv2;

    int64 i1;
    int64 i2;
    int64 k1;
    int64 l1;
    int64 ic;
    int64 ir;
    double rm;
    int64 master[KREDO];
    int64 myleft[KREDO];
    int64 worker[KREDO];
    double buf[2];
    int64 ism[KREDO];
    int64 irm[KREDO];
    double tmp;
    int64 notsend;
    int64 nap1;
    int64 nam1;
    int64 myright[KREDO];
    int64 icp1;

    (void)ncb;

    ipr_dim1 = *nov;
    icf11_dim1 = *nov;
    icf2_dim1 = *nov;
    icf1_dim1 = *nov;
    ca1_dim1 = *nov;
    ca1_dim2 = *nov;
    faa_dim1 = *nov;
    a2_dim1 = *nov;
    a2_dim2 = *nov;
    a1_dim1 = *nov;
    a1_dim2 = *nov;
    cc_dim1 = *nov;
    cc_dim2 = *nrc;

    nbcp1 = *nbc + 1;
    nap1 = *na + 1;
    nam1 = *na - 1;
    xkwt = (double)(*kwt);
    {
        double tmp2 = r_lg10(xkwt) / r_lg10(2.0);
        nlev = i_nint(&tmp2);
    }
    notsend = true;

    // At each recursive level determine the master node (holding the pivot
    /* row after swapping), which will send the pivot row to the worker node
     */
    // at distance 2**(K-1) from the master. Here K is the recursion level.

    if (*par) {
        for (int64 i = 0; i < nlev; ++i) {
            master[i] = false;
            worker[i] = false;
            k1 = pow_ii(2, i);
            niam = *iam / k1;
            if (notsend) {
                if (niam % 2 == 0) {
                    master[i] = true;
                    notsend = false;
                    ism[i] = (i + 1) + *iam + 10000;
                    irm[i] = ism[i] + k1;
                    myright[i] = *iam + k1;
                } else {
                    worker[i] = true;
                    ism[i] = (i + 1) + *iam + 10000;
                    irm[i] = ism[i] - k1;
                    myleft[i] = *iam - k1;
                }
            }
        }
    }

    // Reduce concurrently in each node
    for (i1 = 0; i1 < nam1; ++i1) {
        i2 = i1 + 1;
        for (ic = 0; ic < *nov; ++ic) {
            icp1 = ic + 1;
            ipiv1 = ARRAY2D(ipr, ic, i1);
            if (ipiv1 <= *nov) {
                tmp = ARRAY2D(faa, ic, i1);
                ARRAY2D(faa, ic, i1) = ARRAY2D(faa, (ipiv1 - 1), i1);
                ARRAY2D(faa, (ipiv1 - 1), i1) = tmp;
            } else {
                l1 = (ipiv1 - *nov) - 1;
                tmp = ARRAY2D(faa, ic, i1);
                ARRAY2D(faa, ic, i1) = ARRAY2D(faa, l1, i2);
                ARRAY2D(faa, l1, i2) = tmp;
            }
            for (ir = icp1; ir < *nov; ++ir) {
                l1 = ARRAY2D(icf2, ic, i1) - 1;
                rm = ARRAY3D(a2, ir, l1, i1);
                ARRAY2D(faa, ir, i1) -= rm*ARRAY2D(faa, ic, i1);
            }
            for (ir = 0; ir < *nov; ++ir) {
                l1 = ARRAY2D(icf1, ic, i2) - 1;
                rm = ARRAY3D(a1, ir, l1, i2);
                ARRAY2D(faa, ir, i2) -= rm*ARRAY2D(faa, ic, i1);
            }
            for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
                l1 = ARRAY2D(icf2, ic, i1) - 1;
                rm = ARRAY3D(cc, l1, ir, i2);
                fc[ir] -= rm*ARRAY2D(faa, ic, i1);
            }
        }
    }

    // Inter-node reduction needs communication between nodes
    if (*par) {
        for (int64 i = 0; i < nlev; ++i) {
            for (ic = 0; ic < *nov; ++ic) {
                icp1 = ic + 1;
                if (master[i]) {
                    ipiv1 = ARRAY2D(ipr, ic, (*na - 1));
                    if (ipiv1 <= *nov) {
                        buf[0] = ARRAY2D(faa, (ipiv1 - 1), (*na - 1));
                        ARRAY2D(faa, (ipiv1 - 1), *na) = ARRAY2D(faa, ic, (*na - 1));
                        ARRAY2D(faa, ic, (*na - 1)) = buf[0];
                        buf[1] = -1.;
                        csend();
                    } else {
                        buf[0] = ARRAY2D(faa, ic, (*na - 1));
                        buf[1] = (double)(ARRAY2D(ipr, ic, (*na - 1)) - *nov);
                        csend();
                        crecv();
                    }

                    for (ir = icp1; ir < *nov; ++ir) {
                        l1 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
                        rm = ARRAY3D(a2, ir, l1, (*na - 1));
                        ARRAY2D(faa, ir, (*na - 1)) -= rm*ARRAY2D(faa, ic, (*na - 1));
                    }
                    for (ir = nbcp1 - 1; ir < *nrc; ++ir) {
                        l1 = ARRAY2D(icf2, ic, (*na - 1)) - 1;
                        rm = ARRAY3D(cc, l1, ir, (nap1 - 1));
                        fc[ir] -= rm*ARRAY2D(faa, ic, (*na - 1));
                    }
                }

                if (worker[i]) {
                    crecv();
                    ipiv2 = i_dnnt(&buf[1]);
                    if (ipiv2 < 0) {
                        tmp = buf[0];
                    } else {
                        tmp = ARRAY2D(faa, (ipiv2 - 1), (*na - 1));
                        ARRAY2D(faa, (ipiv2 - 1), (*na - 1)) = buf[0];
                        csend();
                    }

                    for (ir = 0; ir < *nov; ++ir) {
                        l1 = ARRAY2D(icf11, ic, i) - 1;
                        rm = ARRAY3D(ca1, ir, l1, i);
                        ARRAY2D(faa, ir, (*na - 1)) -= rm*tmp;
                    }
                }
            }
            /*           **Synchronization at each recursion level among all n
                         odes */
        }

        l1 = *nrc - *nbc;
        gdsum();
    }

    return 0;
}

int32
dimrge(int64 *iam, int64 *kwt, int64 *par, double *e, double *cc, double *d, double *fc,
       int64 *ifst, int64 *na, int64 *nrc, int64 *nov, int64 *ncb, int64 *idb, int64 *nllv,
       double *fcc, double *p0, double *p1, double *det, double *s, double *a2, double *faa,
       double *bb) {
    int64 e_dim1, cc_dim1, cc_dim2, d_dim1, p0_dim1, p1_dim1, s_dim1, s_dim2, faa_dim1, a2_dim1,
        a2_dim2, bb_dim1, bb_dim2;

    int64 k;

    int64 novpi;
    int64 novpj;
    int64 k1;
    int64 k2;

    int64 novpj2;
    int64 kc;
    int64 kr;
    int64 ncrloc;
    int64 msglen1;
    int64 msglen2;
    int64 nap1;

    double *xe;

    (void)ifst;

    xe = xmalloc(sizeof(*xe)*(usize)(*nov + *nrc));

    faa_dim1 = *nov;
    a2_dim1 = *nov;
    a2_dim2 = *nov;
    s_dim1 = *nov;
    s_dim2 = *nov;
    p1_dim1 = *nov;
    p0_dim1 = *nov;
    cc_dim1 = *nov;
    cc_dim2 = *nrc;
    e_dim1 = *nov + *nrc;
    bb_dim1 = *nov;
    bb_dim2 = *ncb;
    d_dim1 = *ncb;

    nap1 = *na + 1;
    msglen1 = (*nrc*8)**nov;
    // Computing 2nd power
    msglen2 = (*nov + *nrc + ((*nov**nov)*2) + 1)*8;
    ncrloc = *nrc + *nov;

    // Send CC(1:NOV,1:NRC,1) in node 0 to node KWT-1

    if (*par) {
        if (*iam == 0) {
            csend();
        }
        if (*iam == *kwt - 1) {
            crecv();
        }
    }

    // Copy
    if (*iam == *kwt - 1) {
        for (int64 i = 0; i < *nov; ++i) {
            for (int64 j = 0; j < *nov; ++j) {
                novpj = *nov + j;
                ARRAY2D(e, i, j) = ARRAY3D(s, i, j, (*na - 1));
                ARRAY2D(p0, i, j) = ARRAY3D(s, i, j, (*na - 1));
                ARRAY2D(e, i, novpj) = ARRAY3D(a2, i, j, (*na - 1));
                ARRAY2D(p1, i, j) = ARRAY3D(a2, i, j, (*na - 1));
            }
            for (int64 j = 0; j < *ncb; ++j) {
                novpj2 = (*nov*2) + j;
                ARRAY2D(e, i, novpj2) = ARRAY3D(bb, i, j, (*na - 1));
            }
        }

        for (int64 i = 0; i < *nrc; ++i) {
            novpi = *nov + i;
            for (int64 j = 0; j < *nov; ++j) {
                novpj = *nov + j;
                ARRAY2D(e, novpi, j) = ARRAY3D(cc, j, i, 0);
                ARRAY2D(e, novpi, novpj) = ARRAY3D(cc, j, i, (nap1 - 1));
            }
            for (int64 j = 0; j < *ncb; ++j) {
                novpj2 = (*nov*2) + j;
                ARRAY2D(e, novpi, novpj2) = ARRAY2D(d, j, i);
            }
        }

        for (int64 i = 0; i < *nov; ++i) {
            xe[i] = ARRAY2D(faa, i, (*na - 1));
        }

        for (int64 i = 0; i < *nrc; ++i) {
            novpi = *nov + i;
            xe[novpi] = fc[i];
        }

        if (*idb >= 3) {
            fprintf(fp9, " Residuals of reduced system:\n");

            fprintf(fp9, " ");
            for (int64 i = 0; i < ncrloc; ++i) {
                fprintf(fp9, "%11.3E", xe[i]);
                if ((i + 1) % 10 == 0) {
                    fprintf(fp9, "\n ");
                }
            }
            fprintf(fp9, "\n");
        }

        if (*idb >= 4) {
            fprintf(fp9, " Reduced Jacobian matrix:\n");

            for (int64 i = 0; i < ncrloc; ++i) {
                int32 total_printed = 0;
                for (int64 j = 0; j < ncrloc; ++j) {
                    if ((total_printed != 0) && (total_printed % 10 == 0)) {
                        fprintf(fp9, "\n");
                    }
                    fprintf(fp9, " %11.3E", ARRAY2D(e, i, j));
                    total_printed++;
                }
                fprintf(fp9, "\n");
            }
        }

        // Solve for FCC
        if (*nllv == 0) {
            ge(ncrloc, ncrloc, e, 1, ncrloc, fcc, ncrloc, xe, det);
        } else if (*nllv > 0) {
            nlvc(ncrloc, ncrloc, *nllv, e, fcc);
        } else {
            for (int64 i = 0; i < ncrloc - 1; ++i) {
                xe[i] = 0.;
            }
            xe[-1 + ncrloc] = 1.;
            ge(ncrloc, ncrloc, e, 1, ncrloc, fcc, ncrloc, xe, det);
        }
        if (*idb >= 4) {
            fprintf(fp9, " Solution vector:\n");

            for (int64 i = 0; i < ncrloc; ++i) {
                if ((i != 0) && (i % 7 == 0)) {
                    fprintf(fp9, "\n");
                }
                fprintf(fp9, " %11.3E", fcc[i]);
            }
            fprintf(fp9, "\n");
        }

        k1 = ncrloc;
        // Computing 2nd power
        k2 = k1 + (*nov)*(*nov);
        for (kr = 0; kr < *nov; ++kr) {
            for (kc = 0; kc < *nov; ++kc) {
                k = kr**nov + kc;
                fcc[k1 + k] = ARRAY2D(p0, kr, kc);
                fcc[k2 + k] = ARRAY2D(p1, kr, kc);
            }
        }
        // Computing 2nd power
        fcc[ncrloc + ((*nov)*(*nov)*2)] = *det;
    }

    // Broadcast FCC from node KWT-1. The matrices P0 and P1 are
    // buffered in the tail of FCC so all nodes receive them.
    if (*par) {
        if (*iam == *kwt - 1) {
            csend();
        } else {
            crecv();
        }
    }

    for (int64 i = 0; i < *nrc; ++i) {
        fc[i] = fcc[*nov + i];
    }

    if (*iam < *kwt - 1) {
        k1 = ncrloc;
        // Computing 2nd power
        k2 = k1 + (*nov)*(*nov);
        for (kr = 1; kr <= *nov; ++kr) {
            for (kc = 1; kc <= *nov; ++kc) {
                k = kr**nov + kc;
                ARRAY2D(p0, kr, kc) = fcc[k1 + k];
                ARRAY2D(p1, kr, kc) = fcc[k2 + k];
            }
        }
        // Computing 2nd power
        *det = fcc[ncrloc + ((*nov)*(*nov)*2)];
    }
    // free the memory
    /* Not the we have modified these parameter before, so
       we undo the modifications here and then free them. */
    free(xe);

    return 0;
}

int32
bcksub(int64 *iam, int64 *kwt, int64 *par, double *s1, double *s2, double *a2, double *bb,
       double *faa, double *fc, double *fcc, double *sol1, double *sol2, double *sol3, int64 *na,
       int64 *nov, int64 *ncb, int64 *icf2) {
    int64 icf2_dim1, s1_dim1, s1_dim2, s2_dim1, s2_dim2, a2_dim1, a2_dim2, bb_dim1, bb_dim2,
        sol1_dim1, sol2_dim1, sol3_dim1, faa_dim1;

    int64 niam;
    int64 ibuf;
    int64 even = false;
    int64 nlev;
    int64 hasright;
    double xkwt;
    int64 rmsgtype;
    int64 smsgtype;
    int64 k;
    int64 l;

    int64 nlist[2];
    int64 itest;
    int64 l1;
    int64 l2;
    double sm;
    int64 msglen;

    int64 master[KREDO];
    int64 myleft;
    int64 kp1;
    int64 odd = false;
    int64 ism;
    int64 irm;
    int64 hasleft;
    int64 notsend;
    int64 nam1;
    int64 myright;
    int64 nov2;
    int64 nov3;
    double *buf = NULL;

    icf2_dim1 = *nov;
    sol3_dim1 = *nov;
    sol2_dim1 = *nov;
    sol1_dim1 = *nov;
    faa_dim1 = *nov;
    a2_dim1 = *nov;
    a2_dim2 = *nov;
    s2_dim1 = *nov;
    s2_dim2 = *nov;
    s1_dim1 = *nov;
    s1_dim2 = *nov;
    bb_dim1 = *nov;
    bb_dim2 = *ncb;

    xkwt = (double)(*kwt);
    {
        double tmp = d_lg10(&xkwt) / r_lg10(2.0);
        nlev = i_dnnt(&tmp);
    }
    nov2 = *nov*2;
    nov3 = *nov*3;
    ibuf = (nov3 + 1)*8;

    // The backsubstitution in the reduction process is recursive.
    notsend = true;

    /*At each recursion level determine the sender nodes (called MASTER here).
     */
    if (*par) {
        for (int64 i = 0; i < nlev; ++i) {
            master[i] = false;
            niam = *iam / pow_ii(2, i);
            if (notsend) {
                if (niam % 2 == 0) {
                    master[i] = true;
                    notsend = false;
                }
            }
        }
    }

    if (*par) {
        /*Initialization for the master or sender node at the last recursion l
          evel.*/
        if (master[nlev - 1]) {
            for (l = 0; l < *nov; ++l) {
                ARRAY2D(sol1, l, (*na - 1)) = fcc[l];
                ARRAY2D(sol3, l, (*na - 1)) = fc[l];
            }
        }

        for (int64 i = nlev - 1; i >= 0; --i) {
            if (master[i]) {
                ism = i + nlev + (*kwt*4);
                irm = ism + 1;
                k = pow_ii(2, i - 1);
                //              **Compute the ID of the receiving node
                nlist[0] = *iam - k;
                nlist[1] = *iam + k;
                //              **Receive solutions from previous level
                if ((i + 1) < nlev) {
                    crecv();
                    niam = i_dnnt(&buf[nov3 + 1]);
                    if (*iam < niam) {
                        for (l = 0; l < *nov; ++l) {
                            ARRAY2D(sol1, l, (*na - 1)) = buf[l + 1];
                            ARRAY2D(sol3, l, (*na - 1)) = buf[*nov + l + 1];
                        }
                    } else {
                        for (l = 0; l < *nov; ++l) {
                            ARRAY2D(sol1, l, (*na - 1)) = buf[*nov + l + 1];
                            ARRAY2D(sol3, l, (*na - 1)) = buf[nov2 + l + 1];
                        }
                    }
                }
                //              **Backsubstitute
                for (k = *nov - 1; k >= 0; --k) {
                    kp1 = k + 1;
                    sm = 0.;
                    for (l = 0; l < *nov; ++l) {
                        sm += ARRAY3D(s1, k, l, (*na - 1))*ARRAY2D(sol1, l, (*na - 1));
                        sm += ARRAY3D(s2, k, l, (*na - 1))*ARRAY2D(sol3, l, (*na - 1));
                    }
                    for (l = 0; l < *ncb; ++l) {
                        sm += ARRAY3D(bb, k, l, (*na - 1))*fc[*nov + l];
                    }
                    for (l = kp1; l < *nov; ++l) {
                        l1 = ARRAY2D(icf2, l, (*na - 1)) - 1;
                        sm += ARRAY2D(sol2, l1, (*na - 1))*ARRAY3D(a2, k, l1, (*na - 1));
                    }
                    l2 = ARRAY2D(icf2, k, (*na - 1)) - 1;
                    ARRAY2D(sol2, l2, (*na - 1)) =
                        (ARRAY2D(faa, k, (*na - 1)) - sm) / ARRAY3D(a2, k, l2, (*na - 1));
                }
                //              **Send solutions to the next level
                if (i + 1 > 1) {
                    for (l = 0; l < *nov; ++l) {
                        buf[l + 1] = ARRAY2D(sol1, l, (*na - 1));
                        buf[*nov + l + 1] = ARRAY2D(sol2, l, (*na - 1));
                        buf[nov2 + l + 1] = ARRAY2D(sol3, l, (*na - 1));
                    }
                    buf[nov3 + 1] = (double)(*iam);
                    gsendx();
                }
            }
            //           **Synchronization at each recursion level
        }

        // Define odd and even nodes
        if (*iam % 2 == 0) {
            even = true;
        } else {
            odd = true;
        }

        // Determine whether I have a right neighbor
        if (*iam == *kwt - 1) {
            hasright = false;
        } else {
            hasright = true;
        }

        // Determine whether I have a left neighbor
        if (*iam == 0) {
            hasleft = false;
        } else {
            hasleft = true;
        }

        // Define send message type
        smsgtype = *iam + 1000;

        // Define receive message type
        rmsgtype = smsgtype - 1;

        // Define my right neighbor
        myleft = *iam - 1;
        myright = *iam + 1;
        msglen = *nov << 3;

        // May only need odd sends to even
        itest = 0;
        if (itest == 1) {
            if (odd && hasright) {
                csend();
            }
            if (even && hasleft) {
                crecv();
            }
        }

        // Even nodes send and odd nodes receive
        if (even && hasright) {
            csend();
        }
        if (odd && hasleft) {
            crecv();
        }

    } else {

        for (l = 0; l < *nov; ++l) {
            ARRAY2D(sol1, l, (*na - 1)) = fcc[l];
            ARRAY2D(sol2, l, (*na - 1)) = fc[l];
        }
    }

    if (*iam == *kwt - 1) {
        for (l = 0; l < *nov; ++l) {
            ARRAY2D(sol2, l, (*na - 1)) = fc[l];
        }
    }

    if (*na > 1) {
        for (l = 0; l < *nov; ++l) {
            ARRAY2D(sol1, l, (*na - 2)) = ARRAY2D(sol1, l, (*na - 1));
            ARRAY2D(sol3, l, (*na - 2)) = ARRAY2D(sol2, l, (*na - 1));
        }
    }

    // Backsubstitution process; concurrently in each node.
    nam1 = *na - 1;
    for (int64 i = nam1 - 1; i >= 0; --i) {
        for (k = *nov - 1; k >= 0; --k) {
            sm = 0.;
            for (l = 0; l < *nov; ++l) {
                sm += ARRAY2D(sol1, l, i)*ARRAY3D(s1, k, l, i);
                sm += ARRAY2D(sol3, l, i)*ARRAY3D(s2, k, l, i);
            }
            for (l = 0; l < *ncb; ++l) {
                sm += fc[*nov + l]*ARRAY3D(bb, k, l, i);
            }
            for (l = k + 1; l < *nov; ++l) {
                l1 = ARRAY2D(icf2, l, i) - 1;
                sm += ARRAY2D(sol2, l1, i)*ARRAY3D(a2, k, l1, i);
            }
            l2 = ARRAY2D(icf2, k, i) - 1;
            ARRAY2D(sol2, l2, i) = (ARRAY2D(faa, k, i) - sm) / ARRAY3D(a2, k, l2, i);
        }
        for (l = 0; l < *nov; ++l) {
            ARRAY2D(sol1, l, (i + 1)) = ARRAY2D(sol2, l, i);
            if (i + 1 > 1) {
                ARRAY2D(sol3, l, (i - 1)) = ARRAY2D(sol2, l, i);
                ARRAY2D(sol1, l, (i - 1)) = ARRAY2D(sol1, l, i);
            }
        }
    }

    return 0;
}

int32
infpar(int64 *iam, int64 *par, double *a, double *b, double *fa, double *sol1, double *sol2,
       double *fc, int64 *na, int64 *nov, int64 *nra, int64 *nca, int64 *ncb, int64 *irf,
       int64 *icf) {
    int64 irf_dim1, icf_dim1, a_dim1, a_dim2, b_dim1, b_dim2, fa_dim1, sol1_dim1, sol2_dim1;

    int64 nram;
    int64 icfj1;
    double *x;
    int64 nrapj;
    int64 irfir;
    int64 j1;
    int64 novpj;
    int64 icfnovpir;
    int64 ir;
    double sm;
    int64 novpir;
    int64 irp1;

    (void)iam;
    (void)par;

    x = xmalloc(sizeof(*x)*(usize)(*nra));

    // Determine the local varables by backsubstitition.

    sol2_dim1 = *nov;
    sol1_dim1 = *nov;
    irf_dim1 = *nra;
    fa_dim1 = *nra;
    icf_dim1 = *nca;
    a_dim1 = *nca;
    a_dim2 = *nra;
    b_dim1 = *ncb;
    b_dim2 = *nra;

    nram = *nra - *nov;

    /* Backsubstitution in the condensation of parameters; no communication.
     */
    for (int64 i = 0; i < *na; ++i) {
        for (ir = nram - 1; ir >= 0; --ir) {
            irp1 = ir + 1;
            sm = 0.;
            irfir = ARRAY2D(irf, ir, i) - 1;
            for (int64 j = 0; j < *nov; ++j) {
                nrapj = *nra + j;
                sm += ARRAY3D(a, j, irfir, i)*ARRAY2D(sol1, j, i);
                sm += ARRAY3D(a, nrapj, irfir, i)*ARRAY2D(sol2, j, i);
            }
            for (int64 j = 0; j < *ncb; ++j) {
                novpj = *nov + j;
                sm += ARRAY3D(b, j, irfir, i)*fc[novpj];
            }
            for (int64 j = irp1; j < nram; ++j) {
                j1 = j + *nov;
                icfj1 = ARRAY2D(icf, j1, i) - 1;
                sm += ARRAY3D(a, icfj1, irfir, i)*x[icfj1];
            }
            novpir = *nov + ir;
            icfnovpir = ARRAY2D(icf, novpir, i) - 1;
            x[icfnovpir] = (ARRAY2D(fa, irfir, i) - sm) / ARRAY3D(a, icfnovpir, irfir, i);
        }
        //        **Copy SOL1 and X into FA
        for (int64 j = 0; j < *nov; ++j) {
            ARRAY2D(fa, j, i) = ARRAY2D(sol1, j, i);
        }
        for (int64 j = *nov; j < *nra; ++j) {
            ARRAY2D(fa, j, i) = x[j];
        }
    }
    free(x);

    return 0;
}

int32
rd0(int64 *iam, int64 *kwt, double *d, int64 *nrc) {
    int64 niam;
    int64 even[KREDO];
    double xkwt;
    int64 n;

    int64 nredo;
    int64 msglen;
    int64 rmtype[KREDO];
    int64 smtype[KREDO];
    int64 odd[KREDO];

    double *buf;

    int64 notsend;
    int64 myright[KREDO];

    buf = xmalloc(sizeof(*buf)*(usize)(*nrc));

    //     RECURSIVE DOUBLING PROCEDURE TO GET
    //     THE GLOBAL SUM OF VECTORS FROM
    //     EACH NODE. THE GLOBAL SUM IS ONLY AVAILABLE
    //     IN THE LAST NODE

    // Copying
    xkwt = (double)(*kwt);

    // Determine the recursion level
    {
        double tmp = log(xkwt) / log((double)2.);
        nredo = i_dnnt(&tmp);
    }

    // At each recursion level determine the odd and even nodes
    notsend = true;
    for (n = 0; n < nredo; ++n) {
        smtype[n] = n + 1000 + *iam + 1;
        rmtype[n] = smtype[n] - pow_ii(2, n);
        myright[n] = *iam + pow_ii(2, n);
        even[n] = false;
        odd[n] = false;
        niam = *iam / pow_ii(2, n);
        if (notsend) {
            if (niam % 2 == 0) {
                even[n] = true;
                notsend = false;
            } else {
                odd[n] = true;
            }
        }
    }

    niam = *nrc;
    msglen = niam*8;
    for (n = 0; n < nredo; ++n) {
        /*        **Even nodes send and odd nodes receive from left to right
         */
        if (even[n]) {
            csend();
        }
        if (odd[n]) {
            crecv();
            /*          ** Accumulate the partial sum in the current receiving
                        node */
            for (int64 i = 0; i < niam; ++i) {
                d[i] += buf[i];
            }
        }
    }
    free(buf);
    return 0;
}

int32
print1(int64 *nov, int64 *na, int64 *nra, int64 *nca, int64 *ncb, int64 *nrc, double *a, double *b,
       double *c, double *d, double *fa, double *fc) {
    int64 a_dim1;
    int64 a_dim2;
    int64 b_dim1;
    int64 b_dim2;
    int64 c_dim1;
    int64 c_dim2;
    int64 d_dim1;
    int64 fa_dim1;

    int64 ic;
    int64 ir;

    (void)nov;

    fa_dim1 = *nra;
    a_dim1 = *nca;
    a_dim2 = *nra;
    d_dim1 = *ncb;
    b_dim1 = *ncb;
    b_dim2 = *nra;
    c_dim1 = *nca;
    c_dim2 = *nrc;

    fprintf(fp9, "AA , BB , FA (Full dimension) :\n");
    // should be 10.3f
    for (int64 i = 0; i < *na; ++i) {
        fprintf(fp9, "I=%3ld\n", i + 1);
        for (ir = 0; ir < *nra; ++ir) {
            int32 total_written = 0;
            for (ic = 0; ic < *nca; ++ic) {
                if ((total_written != 0) && (total_written % 12 == 0)) {
                    fprintf(fp9, "\n");
                }
                fprintf(fp9, " %10.3E", ARRAY3D(a, ic, ir, i));
                total_written++;
            }
            for (ic = 0; ic < *ncb; ++ic) {
                if ((total_written != 0) && (total_written % 12 == 0)) {
                    fprintf(fp9, "\n");
                }
                fprintf(fp9, " %10.3E", ARRAY3D(b, ic, ir, i));
                total_written++;
            }
            if ((total_written != 0) && (total_written % 12 == 0)) {
                fprintf(fp9, "\n");
            }
            fprintf(fp9, " %10.3E", ARRAY2D(fa, ir, i));
            fprintf(fp9, "\n");
        }
    }

    fprintf(fp9, "CC (Full dimension) :\n");

    for (int64 i = 0; i < *na; ++i) {
        fprintf(fp9, "I=%3ld\n", i + 1);
        for (ir = 0; ir < *nrc; ++ir) {
            int32 total_written = 0;
            for (ic = 0; ic < *nca; ++ic) {
                if ((total_written != 0) && (total_written % 12 == 0)) {
                    fprintf(fp9, "\n");
                }
                fprintf(fp9, " %10.3E", ARRAY3D(c, ic, ir, i));
                total_written++;
            }
            fprintf(fp9, "\n");
        }
    }

    fprintf(fp9, "DD , FC\n");

    for (ir = 0; ir < *nrc; ++ir) {
        int32 total_written = 0;
        for (ic = 0; ic < *ncb; ++ic) {
            if ((total_written != 0) && (total_written % 12 == 0)) {
                fprintf(fp9, "\n");
            }
            fprintf(fp9, " %10.3E", ARRAY2D(d, ic, ir));
            total_written++;
        }
        fprintf(fp9, " %10.3E\n", fc[ir]);
    }

    return 0;
}

/* ----------------------------------------------------------------------- */
/*         Dummy Routines for the Sequential Version */
/* ----------------------------------------------------------------------- */
int64
mynode(void) {
    int64 ret_val;
    ret_val = 0;
    return ret_val;
}

int64
numnodes(void) {
    int64 ret_val;
    ret_val = 1;
    return ret_val;
}

int32
csend(void) {
    return 0;
}

int32
crecv(void) {
    return 0;
}

int32
gdsum(void) {
    return 0;
}

int32
gsendx(void) {
    return 0;
}

int32
gcol(void) {
    return 0;
}

int32
led(void) {
    return 0;
}
