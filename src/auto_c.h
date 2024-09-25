#include <stdio.h>
#include <unistd.h>
#ifdef PTHREADS
#include <pthread.h>
#endif
#ifdef MPI
#include <mpi.h>
#include "auto_mpi.h"
#endif

#include "auto_f2c.h"
#include "somemath.h"

#include "integers.h"

#ifndef AUTO_C_H
#define AUTO_C_H

#define NPARX (36) /*get rid of*/
#define NBIFX (20)
#define KREDO (1)          /*get rid of*/
#define NPARX2 (NPARX*2) /*get rid of*/
#define HMACH (1.0e-7)
#define RSMALL (1.0e-30)
#define RLARGE (1.0e+30)
#define HMACH1 (HMACH + 1.0e0)
#define M1SB (NBIFX)
#define LEFT (1)
#define RIGHT (2)
#define QZMATZ (.FALSE.)
#define QZEPS1 (0.0E0)
#define HMACHHO (1.0e-13)

extern FILE *fp2;
extern FILE *fp3;
extern FILE *fp7;
extern FILE *fp8;
extern FILE *fp9;
extern FILE *fp10;
extern FILE *fp12;

#define CONPAR_DEFAULT 0
#define CONPAR_PTHREADS 1
#define CONPAR_MPI 2
#define SETUBV_DEFAULT 0
#define SETUBV_PTHREADS 1
#define SETUBV_MPI 2

extern int32 global_conpar_type;
extern int32 global_setubv_type;
extern int32 global_num_procs;
extern int32 global_verbose_flag;

typedef struct {
    /* 1 */ int64 ndim;
    /* 2 */ int64 ips;
    /* 3 */ int64 irs;
    /* 4 */ int64 ilp;
    /* 5 */ int64 ntst;
    /* 6 */ int64 ncol;
    /* 7 */ int64 iad;
    /* 8 */ int64 iads;
    /* 9 */ int64 isp;
    /* 10 */ int64 isw;
    /* 11 */ int64 iplt;
    /* 12 */ int64 nbc;
    /* 13 */ int64 nint;
    /* 14 */ int64 nmx;
    /* 15 */ int64 nuzr;
    /* 16 */ int64 npr;
    /* 17 */ int64 mxbf;
    /* 18 */ int64 iid;
    /* 19 */ int64 itmx;
    /* 20 */ int64 itnw;
    /* 21 */ int64 nwtn;
    /* 22 */ int64 jac;
    /* 23 */ int64 ndm;
    /* 24 */ int64 nbc0;
    /* 25 */ int64 nnt0;
    /* 26 */ int64 iuzr;
    /* 27 */ int64 itp;
    /* 28 */ int64 itpst;
    /* 29 */ int64 nfpr;
    /* 30 */ int64 ibr;
    /* 31 */ int64 nit;
    /* 32 */ int64 ntot;
    /* 33 */ int64 nins;
    /* 34 */ int64 istop;
    /* 35 */ int64 nbif;
    /* 36 */ int64 ipos;
    /* 37 */ int64 lab;
    /* 41 */ int64 nicp;
    /* The following are not set in init_.
       They have to do with the old parallel version. */
    /* 38 */ int64 mynode;
    /* 39 */ int64 numnodes;
    /* 40 */ int64 parallel_flag;
} iap_type;

typedef struct {
    /* 1 */ double ds;
    /* 2 */ double dsmin;
    /* 3 */ double dsmax;
    /* There is no 4 */
    /* 5 */ double dsold;
    /* 6 */ double rl0;
    /* 7 */ double rl1;
    /* 8 */ double a0;
    /* 9 */ double a1;
    /* 10 */ double amp;
    /* 11 */ double epsl;
    /* 12 */ double epsu;
    /* 13 */ double epss;
    /* 14 */ double det;
    /* 15 */ double tivp;
    /* 16 */ double fldf;
    /* 17 */ double hbff;
    /* 18 */ double biff;
    /* 19 */ double spbf;
} rap_type;

/*This is the type for all functions which can be used as "funi" the function
  which evaluates the right hand side of the equations and generates
  the Jacobian*/
#define FUNI_TYPE(X)                                                           \
    int32 X(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold, \
            int64 *icp, double *par, int64 ijac, double *f, double *dfdu,      \
            double *dfdp)

/*This is the type for all functions which can be used as "bcni" the function
  which evaluates the boundary conditions */
#define BCNI_TYPE(X)                                                           \
    int32 X(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, \
            int64 nbc, double *u0, double *u1, double *f, int64 ijac,          \
            double *dbc)

/*This is the type for all functions which can be used as "icni" the function
  which evaluates kernel of the integral constraints */
#define ICNI_TYPE(X)                                                           \
    int32 X(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp, \
            int64 nint, double *u, double *uold, double *udot, double *upold,  \
            double *f, int64 ijac, double *dint)

/*This is the type for all functions which can be used as additional
  output functions for algebraic problems */
#define PVLI_TYPE_AE(X)                                                        \
    int32 X(iap_type *iap, rap_type *rap, double *u, double *par)

/*This is the type for all functions which can be used as additional
  output functions for BVPs */
#define PVLI_TYPE_BVP(X)                                                       \
    int32 X(iap_type *iap, rap_type *rap, int64 *icp, double *dtm,             \
            int64 *ndxloc, double *ups, int64 *ndim, double *p0, double *p1,   \
            double *par)

/* This is the type for all functions that can be used at starting points
   for algebraic problems */
#define STPNT_TYPE_AE(X)                                                       \
    int32 X(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *u)

/* This is the type for all functions that can be used at starting points
   for BVPs */
#define STPNT_TYPE_BVP(X)                                                      \
    int32 X(iap_type *iap, rap_type *rap, double *par, int64 *icp,             \
            int64 *ntsrs, int64 *ncolrs, double *rlcur, double *rldot,         \
            int64 *ndxloc, double *ups, double *udotps, double *upoldp,        \
            double *tm, double *dtm, int64 *nodir, double *thl, double *thu)

/*This is the type for all functions which can be used to detect
  special points for algebraic problems */
#define FNCS_TYPE_AE(X)                                                        \
    double X(iap_type *iap, rap_type *rap, double *par, int64 *icp,            \
             int64 *chng, FUNI_TYPE((*funi)), int64 *m1aaloc, double *aa,      \
             double *rlcur, double *rlold, double *rldot, double *u,           \
             double *uold, double *udot, double *rhs, double *dfdu,            \
             double *dfdp, int64 *iuz, double *vuz)

/*This is the type for all functions which can be used to detect
  special points for BVPS */
#define FNCS_TYPE_BVP(X)                                                       \
    double X(iap_type *iap, rap_type *rap, double *par, int64 *icp,            \
             int64 *chng, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)),              \
             ICNI_TYPE((*icni)), double *p0, double *p1, doublecomplex *ev,    \
             double *rlcur, double *rlold, double *rldot, int64 *ndxloc,       \
             double *ups, double *uoldps, double *udotps, double *upoldp,      \
             double *fa, double *fc, double *dups, double *tm, double *dtm,    \
             double *thl, double *thu, int64 *iuz, double *vuz)

#define AUTOAE 0
#define AUTOBV 1

typedef struct {
    FUNI_TYPE((*funi));
    BCNI_TYPE((*bcni));
    ICNI_TYPE((*icni));
    STPNT_TYPE_BVP((*stpnt));
    PVLI_TYPE_BVP((*pvli));
} autobv_function_list;

typedef struct {
    FUNI_TYPE((*funi));
    STPNT_TYPE_AE((*stpnt));
    PVLI_TYPE_AE((*pvli));
} autoae_function_list;

typedef struct {
    int32 type;
    autobv_function_list bvlist;
    autoae_function_list aelist;
} function_list;

double conpar2_time_start(void);
double conpar2_time_end(double);
void allocate_global_memory(iap_type);
int32 init(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *thl,
           double **thu_pointer, int64 *iuz, double *vuz);
int32 autlib_check_dimensions(iap_type *iap);
int32 autoae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), STPNT_TYPE_AE((*stpnt)), PVLI_TYPE_AE((*pvli)),
             double *thl, double *thu, int64 *iuz, double *vuz);

int32 autobv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
             STPNT_TYPE_BVP((*stpnt)), PVLI_TYPE_BVP((*pvli)), double *thl,
             double *thu, int64 *iuz, double *vuz);
int32 autlib1_init(iap_type *iap, rap_type *rap, int64 *icp, double *par);
int32 cnrlae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), STPNT_TYPE_AE((*stpnt)), PVLI_TYPE_AE((*pvli)),
             double *thl, double *thu, int64 *iuz, double *vuz);
STPNT_TYPE_AE(stpnus);
STPNT_TYPE_AE(stpnae);
int32 stprae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), double *rds, int64 *m1aaloc, double *aa,
             double *rhs, double *rlcur, double *rlold, double *rldot,
             double *u, double *du, double *uold, double *udot, double *f,
             double *dfdu, double *dfdp, double *thl, double *thu);
int32 contae(iap_type *iap, rap_type *rap, double *rds, double *rlcur,
             double *rlold, double *rldot, double *u, double *uold,
             double *udot);
int32 solvae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), double *rds, int64 *m1aaloc, double *aa,
             double *rhs, double *rlcur, double *rlold, double *rldot,
             double *u, double *du, double *uold, double *udot, double *f,
             double *dfdu, double *dfdp, double *thl, double *thu);
int32 lcspae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FNCS_TYPE_AE((*fncs)), FUNI_TYPE((*funi)), int64 *m1aaloc,
             double *aa, double *rhs, double *rlcur, double *rlold,
             double *rldot, double *u, double *du, double *uold, double *udot,
             double *f, double *dfdu, double *dfdp, double *q, double *thl,
             double *thu, int64 *iuz, double *vuz);
int32 mueller(double *q0, double *q1, double *q, double *s0, double *s1,
              double *s, double *rds);
FNCS_TYPE_AE(fnbpae);
FNCS_TYPE_AE(fnlpae);
FNCS_TYPE_AE(fnhbae);
FNCS_TYPE_AE(fnuzae);
int32 stbif(iap_type *iap, rap_type *rap, double *par, int64 *icp,
            int64 *m1aaloc, double *aa, int64 m1sbloc, double *stud,
            double *stu, double *stla, double *stld, double *rlcur,
            double *rlold, double *rldot, double *u, double *du, double *udot,
            double *dfdu, double *dfdp, double *thl, double *thu);
int32 swpnt(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *rds,
            int64 m1sbloc, double *stud, double *stu, double *stla,
            double *stld, double *rlcur, double *rlold, double *rldot,
            double *u, double *udot);
int32 sthd(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *thl,
           double *thu);
int32 headng(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 iunit,
             int64 *n1, int64 *n2);
int32 stplae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             double *rlcur, double *u);
int32 wrline(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *icu,
             int64 *ibr, int64 *ntot, int64 *lab, double *vaxis, double *u);
int32 wrtsp8(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *lab,
             double *rlcur, double *u);
int32 wrjac(iap_type *iap, int64 *n, int64 *m1aaloc, double *aa, double *rhs);
int32 msh(iap_type *iap, double *tm);
int32 genwts(int64 ncol, int64 n1, double *wt, double *wp);
int32 cpnts(int64 ncol, double *zm);
int32 cntdif(int64 *n, double *d);
int32 wint(int64 n, double *wi);
int32 adptds(iap_type *iap, rap_type *rap, double *rds);
int32 adapt(iap_type *iap, rap_type *rap, int64 *nold, int64 *ncold,
            int64 *nnew, int64 *ncnew, double *tm, double *dtm, int64 *ndxloc,
            double *ups, double *vps);
int32 interp(iap_type *iap, rap_type *rap, int64 *ndim, int64 *n, int64 *nc,
             double *tm, int64 *ndxloc, double *ups, int64 *n1, int64 *nc1,
             double *tm1, double *ups1, double *tm2, int64 *itm1);
int32 newmsh(iap_type *iap, rap_type *rap, int64 *ndxloc, double *ups,
             int64 *nold, int64 *ncold, double *tmold, double *dtmold,
             int64 *nnew, double *tmnew, int64 *iper);
int32 ordr(iap_type *iap, rap_type *rap, int64 *n, double *tm, int64 *n1,
           double *tm1, int64 *itm1);
int32 intwts(iap_type *iap, rap_type *rap, int64 *n, double *z__, double *x,
             double *wts);
int32 eqdf(iap_type *iap, rap_type *rap, int64 *ntst, int64 *ndim, int64 *ncol,
           double *dtm, int64 *ndxloc, double *ups, double *eqf, int64 *iper);
int32 eig(iap_type *iap, int64 *ndim, int64 *m1a, double *a, doublecomplex *ev,
          int64 *ier);
int32 nlvc(int64 n, int64 m, int64 k, double *a, double *u);
int32 nrmlz(int64 *ndim, double *v);
double pi(double r__);
int32 ge(int64 n, int64 m1a, double *a, int64 nrhs, int64 ndxloc, double *u,
         int64 m1f, double *f, double *det);
int32 newlab(iap_type *iap);
int32 findlb(iap_type *iap, rap_type *rap, int64 irs, int64 *nfpr,
             int64 *found);
int32 readlb(double *u, double *par);
int32 skip3(int64 *nskip, int64 *eof3);
double rinpr(iap_type *iap, int64 *ndim1, int64 *ndxloc, double *ups,
             double *vps, double *dtm, double *thu);
double rnrmsq(iap_type *iap, int64 *ndim1, int64 *ndxloc, double *ups,
              double *dtm, double *thu);
double rintg(iap_type *iap, int64 *ndxloc, int64 ic, double *ups, double *dtm);
double rnrm2(iap_type *iap, int64 *ndxloc, int64 *ic, double *ups, double *dtm);
double rmxups(iap_type *iap, int64 *ndxloc, int64 *i__, double *ups);
double rmnups(iap_type *iap, int64 *ndxloc, int64 *i__, double *ups);
int32 scaleb(iap_type *iap, int64 *icp, int64 *ndxloc, double *dvps,
             double *rld, double *dtm, double *thl, double *thu);
int32 cnrlbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
             STPNT_TYPE_BVP((*stpnt)), PVLI_TYPE_BVP((*pvli)), double *thl,
             double *thu, int64 *iuz, double *vuz);
int32 contbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), double *rds, double *rlcur, double *rlold,
             double *rldot, int64 *ndxloc, double *ups, double *uoldps,
             double *udotps, double *upoldp, double *dtm, double *thl,
             double *thu);
int32 extrbv(iap_type *iap, rap_type *rap, FUNI_TYPE((*funi)), double *rds,
             double *rlcur, double *rlold, double *rldot, int64 *ndxloc,
             double *ups, double *uoldps, double *udotps);
int32 stupbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), double *rlcur, double *rlold, double *rldot,
             int64 *ndxloc, double *ups, double *uoldps, double *upoldp);
int32 stepbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
             PVLI_TYPE_BVP((*pvli)), double *rds, double *rlcur, double *rlold,
             double *rldot, int64 *ndxloc, double *ups, double *dups,
             double *uoldps, double *udotps, double *upoldp, double *fa,
             double *fc, double *tm, double *dtm, double *p0, double *p1,
             double *thl, double *thu);
int32 rsptbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), STPNT_TYPE_BVP((*stpnt)), double *rds,
             double *rlcur, double *rlold, double *rldot, int64 *ndxloc,
             double *ups, double *uoldps, double *udotps, double *upoldp,
             double *dups, double *tm, double *dtm, doublecomplex *ev,
             int64 *nodir, double *thl, double *thu);
STPNT_TYPE_BVP(stpnbv);
STPNT_TYPE_BVP(stpnub);
int32 setrtn(iap_type *iap, int64 *ntst, int64 *ndxloc, double *ups,
             double *par);
int32 stdrbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
             double *rlcur, double *rlold, double *rldot, int64 ndxloc,
             double *ups, double *dups, double *uoldps, double *udotps,
             double *upoldp, double *fa, double *fc, double *dtm, int64 iperp,
             double *p0, double *p1, double *thl, double *thu);
int32 lcspbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FNCS_TYPE_BVP((*fncs)), FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)),
             ICNI_TYPE((*icni)), PVLI_TYPE_BVP((*pvli)), double *q,
             double *rlcur, double *rlold, double *rldot, int64 *ndxloc,
             double *ups, double *dups, double *uoldps, double *udotps,
             double *upoldp, double *fa, double *fc, double *tm, double *dtm,
             double *p0, double *p1, doublecomplex *ev, double *thl,
             double *thu, int64 *iuz, double *vuz);
FNCS_TYPE_BVP(fnlpbv);
FNCS_TYPE_BVP(fnbpbv);
FNCS_TYPE_BVP(fnspbv);
FNCS_TYPE_BVP(fnuzbv);
int32 tpspbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             doublecomplex *ev);
int32 stplbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             double *rldot, int64 *ndxloc, double *ups, double *udotps,
             double *tm, double *dtm, double *thl, double *thu);
int32 wrtbv8(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             double *rldot, int64 *ndxloc, double *ups, double *udotps,
             double *tm, double *dtm);
int32 wrtbv9(iap_type *iap, rap_type *rap, double *par, int64 *icp,
             double *rlcur, int64 *ndxloc, double *ups, double *tm, double *dtm,
             double *thl, double *thu);
PVLI_TYPE_AE(pvlsae);
PVLI_TYPE_BVP(pvlsbv);
int32 setpae(iap_type *iap, rap_type *rap);
int32 setpbv(iap_type *iap, rap_type *rap, double *dtm);
int32 autim0(double *t);
int32 autim1(double *t);
double getp(char *code, int64 *ic, double *ups, int64 code_len);
int32 solvbv(int64 *ifst, iap_type *iap, rap_type *rap, double *par, int64 *icp,
             FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
             double *rds, int64 *nllv, double *rlcur, double *rlold,
             double *rldot, int64 *ndxloc, double *ups, double *dups,
             double *uoldps, double *udotps, double *upoldp, double *dtm,
             double *fa, double *fc, double *p0, double *p1, double *thl,
             double *thu);
int32 setfcdd(int64 *ifst, double *dd, double *fc, int64 *ncb, int64 *nrc);
int32 faft(double *ff, double *fa, int64 *ntst, int64 *nrow, int64 *ndxloc);
int32 partition(int64 *n, int64 *kwt, int64 *m);
int64 mypart(int64 *iam, int64 *np);
int32 setrhs(int64 *ndim, int64 *ips, int64 *na, int64 *ntst, int64 *np,
             int64 *ncol, int64 *nbc, int64 *nint, int64 *ncb, int64 *nrc,
             int64 *nra, int64 *nca, int64 *iam, int64 *kwt, int64 *ipar,
             FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
             int64 *ndxloc, iap_type *iap, rap_type *rap, double *par,
             int64 *icp, double *rds, double *fa, double *fc, double *rlcur,
             double *rlold, double *rldot, double *ups, double *uoldps,
             double *udotps, double *upoldp, double *dups, double *dtm,
             double *thl, double *thu, double *p0, double *p1);
int32 brbd(double *a, double *b, double *c, double *d, double *fa, double *fc,
           double *p0, double *p1, int64 *ifst, int64 *idb, int64 *nllv,
           double *det, int64 *nov, int64 *na, int64 *nbc, int64 *nra,
           int64 *nca, int64 *ncb, int64 *nrc, int64 *iam, int64 *kwt,
           int64 *par, double *a1, double *a2, double *bb, double *cc,
           double *faa, double *ca1, double *s1, double *s2, int64 *icf11,
           int64 *ipr, int64 *icf1, int64 *icf2, int64 *irf, int64 *icf);
int32 setzero(double *fa, double *fc, int64 *na, int64 *nra, int64 *nrc);
int32 conrhs(int64 *nov, int64 *na, int64 *nra, int64 *nca, double *a,
             int64 *nbc, int64 *nrc, double *c, double *fa, double *fc,
             int64 *irf, int64 *icf, int64 *iam);
int32 copycp(int64 *iam, int64 *kwt, int64 *na, int64 *nov, int64 *nra,
             int64 *nca, double *a, int64 *ncb, double *b, int64 *nrc,
             double *c, double *a1, double *a2, double *bb, double *cc,
             int64 *irf);
int32 cpyrhs(int64 *na, int64 *nov, int64 *nra, double *faa, double *fa,
             int64 *irf);
int32 reduce(int64 *iam, int64 *kwt, int64 *par, double *a1, double *a2,
             double *bb, double *cc, double *dd, int64 *na, int64 *nov,
             int64 *ncb, int64 *nrc, double *s1, double *s2, double *ca1,
             int64 *icf1, int64 *icf2, int64 *icf11, int64 *ipr, int64 *nbc);
int32 redrhs(int64 *iam, int64 *kwt, int64 *par, double *a1, double *a2,
             double *cc, double *faa, double *fc, int64 *na, int64 *nov,
             int64 *ncb, int64 *nrc, double *ca1, int64 *icf1, int64 *icf2,
             int64 *icf11, int64 *ipr, int64 *nbc);
int32 dimrge(int64 *iam, int64 *kwt, int64 *par, double *e, double *cc,
             double *d, double *fc, int64 *ifst, int64 *na, int64 *nrc,
             int64 *nov, int64 *ncb, int64 *idb, int64 *nllv, double *fcc,
             double *p0, double *p1, double *det, double *s, double *a2,
             double *faa, double *bb);
int32 bcksub(int64 *iam, int64 *kwt, int64 *par, double *s1, double *s2,
             double *a2, double *bb, double *faa, double *fc, double *fcc,
             double *sol1, double *sol2, double *sol3, int64 *na, int64 *nov,
             int64 *ncb, int64 *icf2);
int32 infpar(int64 *iam, int64 *par, double *a, double *b, double *fa,
             double *sol1, double *sol2, double *fc, int64 *na, int64 *nov,
             int64 *nra, int64 *nca, int64 *ncb, int64 *irf, int64 *icf);
int32 rd0(int64 *iam, int64 *kwt, double *d, int64 *nrc);
int32 print1(int64 *nov, int64 *na, int64 *nra, int64 *nca, int64 *ncb,
             int64 *nrc, double *a, double *b, double *c, double *d, double *fa,
             double *fc);
int64 mynode(void);
int64 numnodes(void);
int32 csend(void);
int32 crecv(void);
int32 gcol(void);
int32 led(void);
FUNI_TYPE(fnlp);
int32 fflp(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, double *f, int64 ndm, double *dfdu,
           double *dfdp);
STPNT_TYPE_AE(stpnlp);
FUNI_TYPE(fnc1);
STPNT_TYPE_AE(stpnc1);
FUNI_TYPE(fnc2);
int32 ffc2(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, double *f, int64 ndm, double *dfdu,
           double *dfdp);
STPNT_TYPE_AE(stpnc2);
FUNI_TYPE(fnds);
FUNI_TYPE(fnti);
FUNI_TYPE(fnhd);
int32 ffhd(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, double *f, int64 ndm, double *dfdu,
           double *dfdp);
STPNT_TYPE_AE(stpnhd);
FUNI_TYPE(fnhb);
int32 ffhb(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, double *f, int64 ndm, double *dfdu,
           double *dfdp);
STPNT_TYPE_AE(stpnhb);
FUNI_TYPE(fnhw);
int32 ffhw(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, double *f, int64 ndm, double *dfdu,
           double *dfdp);
STPNT_TYPE_AE(stpnhw);
FUNI_TYPE(fnps);
BCNI_TYPE(bcps);
ICNI_TYPE(icps);
int32 pdble(iap_type *iap, rap_type *rap, int64 *ndim, int64 *ntst, int64 *ncol,
            int64 *ndxloc, double *ups, double *udotps, double *tm,
            double *par);
STPNT_TYPE_BVP(stpnps);
FUNI_TYPE(fnws);
int32 ffws(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, int64 ijac, double *f, double *dfdu,
           double *dfdp, int64 ndm, double *dfu, double *dfp);
FUNI_TYPE(fnwp);
int32 stpnwp(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr,
             int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,
             double *ups, double *udotps, double *upoldp, double *tm,
             double *dtm, int64 *nodir, double *thl, double *thu);
FUNI_TYPE(fnsp);
int32 ffsp(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, int64 ijac, double *f, double *dfdu,
           double *dfdp, int64 ndm, double *dfu, double *dfp);
FUNI_TYPE(fnpe);
int32 ffpe(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, int64 ijac, double *f, double *dfdu,
           double *dfdp, int64 ndm, double *dfu, double *dfp);
ICNI_TYPE(icpe);
FUNI_TYPE(fnpl);
int32 ffpl(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, double *f, int64 ndm, double *dfdu,
           double *dfdp);
BCNI_TYPE(bcpl);
ICNI_TYPE(icpl);
int32 stpnpl(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr,
             int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,
             double *ups, double *udotps, double *upoldp, double *tm,
             double *dtm, int64 *nodir, double *thl, double *thu);
FUNI_TYPE(fnpd);
int32 ffpd(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, double *f, int64 ndm, double *dfdu,
           double *dfdp);
BCNI_TYPE(bcpd);
ICNI_TYPE(icpd);
int32 stpnpd(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr,
             int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,
             double *ups, double *udotps, double *upoldp, double *tm,
             double *dtm, int64 *nodir, double *thl, double *thu);
FUNI_TYPE(fntr);
int32 fftr(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, double *f, int64 ndm, double *dfdu,
           double *dfdp);
int32 bctr(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp,
           int64 nbc, double *u0, double *u1, double *f, int64 ijac,
           double *dbc);
ICNI_TYPE(ictr);
int32 stpntr(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr,
             int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,
             double *ups, double *udotps, double *upoldp, double *tm,
             double *dtm, int64 *nodir, double *thl, double *thu);
FUNI_TYPE(fnpo);
int32 ffpo(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           double *upold, int64 *icp, double *par, double *f, int64 ndm,
           double *dfdu, double *dfdp);
BCNI_TYPE(bcpo);
ICNI_TYPE(icpo);
int32 fipo(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp,
           int64 nint, int64 nnt0, double *u, double *uold, double *udot,
           double *upold, double *fi, double *dint, int64 ndmt, double *dfdu,
           double *dfdp);
STPNT_TYPE_BVP(stpnpo);
FUNI_TYPE(fnbl);
int32 ffbl(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, double *f, int64 ndm, double *dfdu,
           double *dfdp);
BCNI_TYPE(bcbl);
int32 fbbl(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp,
           int64 nbc, int64 nbc0, double *u0, double *u1, double *f,
           double *dbc);
ICNI_TYPE(icbl);
int32 fibl(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp,
           int64 nint, int64 nnt0, double *u, double *uold, double *udot,
           double *upold, double *f, double *dint);
STPNT_TYPE_BVP(stpnbl);
FUNI_TYPE(funi);
BCNI_TYPE(bcni);
ICNI_TYPE(icni);
int32 fopi(iap_type *iap, rap_type *rap, int64 ndim, double *u, int64 *icp,
           double *par, int64 ijac, double *f, double *dfdu, double *dfdp);
int32 flowkm(int64 *ndim, double *c0, double *c1, int64 *iid, double *rwork,
             doublecomplex *ev);
int32 dhhpr(int64 *k, int64 *j, int64 *n, double *x, int64 *incx, double *beta,
            double *v);
int32 dhhap(int64 *k, int64 *j, int64 *n, int64 *q, double *beta, double *v,
            int64 *job, double *a, int64 *lda);
FUNI_TYPE(fnho);
int32 ffho(iap_type *iap, rap_type *rap, int64 ndim, double *u, double *uold,
           int64 *icp, double *par, double *f, int64 ndm, double *dfdu,
           double *dfdp);
BCNI_TYPE(bcho);
int32 fbho(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp,
           int64 nbc, int64 nbc0, double *u0, double *u1, double *fb,
           double *dbc);
ICNI_TYPE(icho);
int32 fiho(iap_type *iap, rap_type *rap, int64 ndim, double *par, int64 *icp,
           int64 nint, int64 nnt0, double *u, double *uold, double *udot,
           double *upold, double *fi, double *dint);
int32 inho(iap_type *iap, int64 *icp, double *par);
int32 preho(int64 *ndx, int64 *ntsr, int64 *nar, int64 *ndim, int64 *ncolrs,
            double *ups, double *udotps, double *tm, double *par);
int32 stpnho(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr,
             int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,
             double *ups, double *udotps, double *upoldp, double *tm,
             double *dtm, int64 *nodir, double *thl, double *thu);
int32 stpho(iap_type *iap, int64 *icp, double *u, double *par, double *t);
PVLI_TYPE_BVP(pvlsho);
double psiho(iap_type *iap, int64 is, double *rr, double *ri, double *v,
             double *vt, int64 *icp, double *par);
int32 eighi(int64 isign, int64 itrans, double *rr, double *ri, double *vret,
            double *xequib, int64 *icp, double *par, int64 *ndm);
int32 eigho(int64 *isign, int64 *itrans, double *rr, double *ri, double *vret,
            double *xequib, int64 *icp, double *par, int64 *ndm, double *dfdu,
            double *dfdp, double *zz);
int32 prjcti(double *bound, double *xequib, int64 *icp, double *par, int64 imfd,
             int64 is, int64 itrans, int64 *ndm);
int32 prjctn(double *bound, double *xequib, int64 *icp, double *par,
             int64 *imfd, int64 *is, int64 *itrans, int64 *ndm, double *dfdu,
             double *dfdp);
int32 rg(int64 nm, int64 n, double *a, double *wr, double *wi, int64 matz,
         double *z__, int64 *iv1, double *fv1, int64 *ierr);
int32 hqr(int64 *nm, int64 *n, int64 *low, int64 *igh, double *h__, double *wr,
          double *wi, int64 *ierr);
int32 hqr2(int64 *nm, int64 *n, int64 *low, int64 *igh, double *h__, double *wr,
           double *wi, double *z__, int64 *ierr);
int32 cdiv(double *ar, double *ai, double *br, double *bi, double *cr,
           double *ci);
int32 balanc(int64 *nm, int64 *n, double *a, int64 *low, int64 *igh,
             double *scale);
int32 balbak(int64 *nm, int64 *n, int64 *low, int64 *igh, double *scale,
             int64 *m, double *z__);
int32 elmhes(int64 *nm, int64 *n, int64 *low, int64 *igh, double *a,
             int64 *int__);
int32 eltran(int64 *nm, int64 *n, int64 *low, int64 *igh, double *a,
             int64 *int__, double *z__);
int32 qzhes(int64 nm, int64 n, double *a, double *b, int64 matz, double *z__);
int32 qzit(int64 nm, int64 n, double *a, double *b, double eps1, int64 matz,
           double *z__, int64 *ierr);
int32 qzval(int64 nm, int64 n, double *a, double *b, double *alfr, double *alfi,
            double *beta, int64 matz, double *z__);
double epslon(double x);
double dnrm2(int64 *n, double *dx, int64 *incx);
double ddot(int64 *n, double *dx, int64 *incx, double *dy, int64 *incy);
int32 dscal(int64 *n, double *da, double *dx, int64 *incx);
int64 idamax(int64 *n, double *dx, int64 *incx);
int32 daxpy(int64 *n, double *da, double *dx, int64 *incx, double *dy,
            int64 *incy);
int32 drot(int64 *n, double *dx, int64 *incx, double *dy, int64 *incy,
           double *c, double *s);
int32 dswap(int64 *n, double *dx, int64 *incx, double *dy, int64 *incy);
int32 dgemc(int64 *m, int64 *n, double *a, int64 *lda, double *b, int64 *ldb,
            int64 *trans);
int32 xerbla(char *srname, int64 *info, int64 srname_len);
int64 lsame(char *ca, char *cb, int64 ca_len, int64 cb_len);
int32 dgemm(char *transa, char *transb, int64 *m, int64 *n, int64 *k,
            double *alpha, double *a, int64 *lda, double *b, int64 *ldb,
            double *beta, double *c, int64 *ldc, int64 transa_len,
            int64 transb_len);
int32 ezsvd(double *x, int64 *ldx, int64 *n, int64 *p, double *s, double *e,
            double *u, int64 *ldu, double *v, int64 *ldv, double *work,
            int64 *job, int64 *info, double *tol);
int32 ndrotg(double *f, double *g, double *cs, double *sn);
int32 ndsvd(double *x, int64 *ldx, int64 *n, int64 *p, double *s, double *e,
            double *u, int64 *ldu, double *v, int64 *ldv, double *work,
            int64 *job, int64 *info, int64 *maxitr, double *tol, int64 *idbg,
            int64 *ifull, int64 *kount, int64 *kount1, int64 *kount2,
            int64 *skip, int64 *limshf, double *maxsin, int64 *iidir);
int32 prse(int64 *ll, int64 *m, int64 *nrow, int64 *ncol, double *s, double *e);
int32 sig22(double *a, double *b, double *c, double *sigmin, double *sigmax,
            double *snr, double *csr, double *snl, double *csl);
double sigmin(double *a, double *b, double *c);
int32 sndrtg(double *f, double *g, double *cs, double *sn);
int32 hqr3lc(double *a, double *v, int64 *n, int64 *nlow, int64 *nup,
             double *eps, double *er, double *ei, int64 *type__, int64 *na,
             int64 *nv, int64 *imfd);
int32 split(double *a, double *v, int64 *n, int64 *l, double *e1, double *e2,
            int64 *na, int64 *nv);
int32 exchng(double *a, double *v, int64 *n, int64 *l, int64 *b1, int64 *b2,
             double *eps, int64 *fail, int64 *na, int64 *nv);
int32 qrstep(double *a, double *v, double *p, double *q, double *r__, int64 *nl,
             int64 *nu, int64 *n, int64 *na, int64 *nv);
int32 orthes(int64 *nm, int64 *n, int64 *low, int64 *igh, double *a,
             double *ort);
int32 ortran(int64 *nm, int64 *n, int64 *low, int64 *igh, double *a,
             double *ort, double *z__);

/* problem defined functions*/
int32 func(int64 ndim, double *u, int64 *icp, double *par, int64 ijac,
           double *f, double *dfdu, double *dfdp);
int32 stpnt(int64 ndim, double t, double *u, double *par);
int32 bcnd(int64 ndim, double *par, int64 *icp, int64 nbc, double *u0,
           double *u1, int64 ijac, double *f, double *dbc);
int32 icnd(int64 ndim, double *par, int64 *icp, int64 nint, double *u,
           double *uold, double *udot, double *upold, int64 ijac, double *fi,
           double *dint);
int32 fopt(int64 ndim, double *u, int64 *icp, double *par, int64 ijac,
           double *fs, double *dfdu, double *dfdp);
int32 pvls(int64 ndim, double *u, double *par);
void *conpar2_process(void *);
int32 conpar2(int64 *nov, int64 *na, int64 *nra, int64 *nca, double *a,
              int64 *ncb, double *b, int64 *nbc, int64 *nrc, double *c,
              double *d, int64 *irf, int64 *icf);
/*setubv.c */
#include "auto_types.h"
void *setubv_make_aa_bb_cc(void *);
int32 setubv(int64 ndim, int64 ips, int64 na, int64 ncol, int64 nbc, int64 nint,
             int64 ncb, int64 nrc, int64 nra, int64 nca, FUNI_TYPE((*funi)),
             BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), int64 ndxloc,
             iap_type *iap, rap_type *rap, double *par, int64 *icp, double rds,
             double *aa, double *bb, double *cc, double *dd, double *fa,
             double *fc, double *rlcur, double *rlold, double *rldot,
             double *ups, double *uoldps, double *udotps, double *upoldp,
             double *dups, double *dtm, double *thl, double *thu, double *p0,
             double *p1);
void setubv_parallel_arglist_copy(setubv_parallel_arglist *output,
                                  setubv_parallel_arglist input);
void setubv_parallel_arglist_constructor(
    int64 ndim, int64 ips, int64 na, int64 ncol, int64 nbc, int64 nint,
    int64 ncb, int64 nrc, int64 nra, int64 nca, FUNI_TYPE((*funi)),
    ICNI_TYPE((*icni)), int64 ndxloc, iap_type *iap, rap_type *rap, double *par,
    int64 *icp, double *aa, double *bb, double *cc, double *dd, double *fa,
    double *fc, double *ups, double *uoldps, double *udotps, double *upoldp,
    double *dtm, double *wp, double *wt, double *wi, double *thu, double *thl,
    double *rldot, BCNI_TYPE((*bcni)), setubv_parallel_arglist *data);
void setubv_make_fa(setubv_parallel_arglist larg);
void setubv_make_fc_dd(setubv_parallel_arglist larg, double *dups,
                       double *rlcur, double *rlold, double rds);

#include "auto_types.h"
int32 set_funi_and_icni(iap_type *, setubv_parallel_arglist *);
int32 set_function_pointers(iap_type, function_list *);

#ifdef AUTO_CONSTRUCT_DESCTRUCT
int32 user_construct(int32 argc, char **argv);
int32 user_destruct();
#endif

#endif
