#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#ifdef PTHREADS
#include <pthread.h>
#endif
#ifdef MPI
#include <mpi.h>
#include "auto_mpi.h"
#endif
#include <string.h>

#ifndef __AUTO_C_H__
#define __AUTO_C_H__

#define NPARX (36) /*get rid of*/
#define NBIFX (20)
#define KREDO (1)          /*get rid of*/
#define NPARX2 (NPARX * 2) /*get rid of*/
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

extern int global_conpar_type;
extern int global_setubv_type;
extern int global_num_procs;
extern int global_verbose_flag;

typedef struct {
    /* 1 */ integer ndim;
    /* 2 */ integer ips;
    /* 3 */ integer irs;
    /* 4 */ integer ilp;
    /* 5 */ integer ntst;
    /* 6 */ integer ncol;
    /* 7 */ integer iad;
    /* 8 */ integer iads;
    /* 9 */ integer isp;
    /* 10 */ integer isw;
    /* 11 */ integer iplt;
    /* 12 */ integer nbc;
    /* 13 */ integer nint;
    /* 14 */ integer nmx;
    /* 15 */ integer nuzr;
    /* 16 */ integer npr;
    /* 17 */ integer mxbf;
    /* 18 */ integer iid;
    /* 19 */ integer itmx;
    /* 20 */ integer itnw;
    /* 21 */ integer nwtn;
    /* 22 */ integer jac;
    /* 23 */ integer ndm;
    /* 24 */ integer nbc0;
    /* 25 */ integer nnt0;
    /* 26 */ integer iuzr;
    /* 27 */ integer itp;
    /* 28 */ integer itpst;
    /* 29 */ integer nfpr;
    /* 30 */ integer ibr;
    /* 31 */ integer nit;
    /* 32 */ integer ntot;
    /* 33 */ integer nins;
    /* 34 */ integer istop;
    /* 35 */ integer nbif;
    /* 36 */ integer ipos;
    /* 37 */ integer lab;
    /* 41 */ integer nicp;
    /* The following are not set in init_.
       They have to do with the old parallel version. */
    /* 38 */ integer mynode;
    /* 39 */ integer numnodes;
    /* 40 */ integer parallel_flag;
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
    int X(const iap_type *iap, const rap_type *rap, integer ndim,              \
          const double *u, const double *uold, const integer *icp,     \
          double *par, integer ijac, double *f, double *dfdu,      \
          double *dfdp)

/*This is the type for all functions which can be used as "bcni" the function
  which evaluates the boundary conditions */
#define BCNI_TYPE(X)                                                           \
    int X(const iap_type *iap, const rap_type *rap, integer ndim,              \
          double *par, const integer *icp, integer nbc,                    \
          const double *u0, const double *u1, double *f,           \
          integer ijac, double *dbc)

/*This is the type for all functions which can be used as "icni" the function
  which evaluates kernel of the integral constraints */
#define ICNI_TYPE(X)                                                           \
    int X(const iap_type *iap, const rap_type *rap, integer ndim,              \
          double *par, const integer *icp, integer nint,                   \
          const double *u, const double *uold, const double *udot, \
          const double *upold, double *f, integer ijac,                \
          double *dint)

/*This is the type for all functions which can be used as additional
  output functions for algebraic problems */
#define PVLI_TYPE_AE(X)                                                        \
    int X(iap_type *iap, rap_type *rap, double *u, double *par)

/*This is the type for all functions which can be used as additional
  output functions for BVPs */
#define PVLI_TYPE_BVP(X)                                                       \
    int X(iap_type *iap, rap_type *rap, integer *icp, double *dtm,         \
          integer *ndxloc, double *ups, integer *ndim, double *p0,     \
          double *p1, double *par)

/* This is the type for all functions that can be used at starting points
   for algebraic problems */
#define STPNT_TYPE_AE(X)                                                       \
    int X(iap_type *iap, rap_type *rap, double *par, integer *icp,         \
          double *u)

/* This is the type for all functions that can be used at starting points
   for BVPs */
#define STPNT_TYPE_BVP(X)                                                      \
    int X(iap_type *iap, rap_type *rap, double *par, integer *icp,         \
          integer *ntsrs, integer *ncolrs, double *rlcur,                  \
          double *rldot, integer *ndxloc, double *ups,                 \
          double *udotps, double *upoldp, double *tm,              \
          double *dtm, integer *nodir, double *thl, double *thu)

/*This is the type for all functions which can be used to detect
  special points for algebraic problems */
#define FNCS_TYPE_AE(X)                                                        \
    double X(iap_type *iap, rap_type *rap, double *par, integer *icp,  \
                 logical *chng, FUNI_TYPE((*funi)), integer *m1aaloc,          \
                 double *aa, double *rlcur, double *rlold,         \
                 double *rldot, double *u, double *uold,           \
                 double *udot, double *rhs, double *dfdu,          \
                 double *dfdp, integer *iuz, double *vuz)

/*This is the type for all functions which can be used to detect
  special points for BVPS */
#define FNCS_TYPE_BVP(X)                                                       \
    double X(iap_type *iap, rap_type *rap, double *par, integer *icp,  \
                 logical *chng, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)),        \
                 ICNI_TYPE((*icni)), double *p0, double *p1,           \
                 doublecomplex *ev, double *rlcur, double *rlold,      \
                 double *rldot, integer *ndxloc, double *ups,          \
                 double *uoldps, double *udotps, double *upoldp,   \
                 double *fa, double *fc, double *dups,             \
                 double *tm, double *dtm, double *thl,             \
                 double *thu, integer *iuz, double *vuz)

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
    int type;
    autobv_function_list bvlist;
    autoae_function_list aelist;
} function_list;

/* main.c */
int main();

/* autlib1.c */
double time_start(void);
double time_end(double);
void allocate_global_memory(const iap_type);
int init(iap_type *iap, rap_type *rap, double *par, integer *icp,
         double *thl, double **thu_pointer, integer *iuz,
         double *vuz);
int chdim(iap_type *iap);
int autoae(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), STPNT_TYPE_AE((*stpnt)), PVLI_TYPE_AE((*pvli)),
           double *thl, double *thu, integer *iuz, double *vuz);

int autobv(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           STPNT_TYPE_BVP((*stpnt)), PVLI_TYPE_BVP((*pvli)), double *thl,
           double *thu, integer *iuz, double *vuz);
int init1(iap_type *iap, rap_type *rap, integer *icp, double *par);
int cnrlae(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), STPNT_TYPE_AE((*stpnt)), PVLI_TYPE_AE((*pvli)),
           double *thl, double *thu, integer *iuz, double *vuz);
STPNT_TYPE_AE(stpnus);
STPNT_TYPE_AE(stpnae);
int stprae(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), double *rds, integer *m1aaloc,
           double *aa, double *rhs, double *rlcur,
           double *rlold, double *rldot, double *u, double *du,
           double *uold, double *udot, double *f, double *dfdu,
           double *dfdp, double *thl, double *thu);
int contae(iap_type *iap, rap_type *rap, double *rds, double *rlcur,
           double *rlold, double *rldot, double *u,
           double *uold, double *udot);
int solvae(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), double *rds, integer *m1aaloc,
           double *aa, double *rhs, double *rlcur,
           double *rlold, double *rldot, double *u, double *du,
           double *uold, double *udot, double *f, double *dfdu,
           double *dfdp, double *thl, double *thu);
int lcspae(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FNCS_TYPE_AE((*fncs)), FUNI_TYPE((*funi)), integer *m1aaloc,
           double *aa, double *rhs, double *rlcur,
           double *rlold, double *rldot, double *u, double *du,
           double *uold, double *udot, double *f, double *dfdu,
           double *dfdp, double *q, double *thl, double *thu,
           integer *iuz, double *vuz);
int mueller(double *q0, double *q1, double *q, double *s0,
            double *s1, double *s, double *rds);
FNCS_TYPE_AE(fnbpae);
FNCS_TYPE_AE(fnlpae);
FNCS_TYPE_AE(fnhbae);
FNCS_TYPE_AE(fnuzae);
int stbif(iap_type *iap, rap_type *rap, double *par, integer *icp,
          integer *m1aaloc, double *aa, integer m1sbloc, double *stud,
          double *stu, double *stla, double *stld,
          double *rlcur, double *rlold, double *rldot,
          double *u, double *du, double *udot, double *dfdu,
          double *dfdp, double *thl, double *thu);
int swpnt(iap_type *iap, rap_type *rap, double *par, integer *icp,
          double *rds, integer m1sbloc, double *stud, double *stu,
          double *stla, double *stld, double *rlcur,
          double *rlold, double *rldot, double *u,
          double *udot);
int swprc(iap_type *iap, rap_type *rap, double *par, integer *icp,
          FUNI_TYPE((*funi)), integer *m1aaloc, double *aa, double *rhs,
          double *rlcur, double *rlold, double *rldot,
          double *u, double *du, double *uold, double *udot,
          double *f, double *dfdu, double *dfdp, double *rds,
          double *thl, double *thu);
int sthd(iap_type *iap, rap_type *rap, double *par, integer *icp,
         double *thl, double *thu);
int headng(iap_type *iap, rap_type *rap, double *par, integer *icp,
           integer iunit, integer *n1, integer *n2);
int stplae(iap_type *iap, rap_type *rap, double *par, integer *icp,
           double *rlcur, double *u);
int wrline(iap_type *iap, rap_type *rap, double *par, integer *icp,
           integer *icu, integer *ibr, integer *ntot, integer *lab,
           double *vaxis, double *u);
int wrtsp8(iap_type *iap, rap_type *rap, double *par, integer *icp,
           integer *lab, double *rlcur, double *u);
int wrjac(iap_type *iap, integer *n, integer *m1aaloc, double *aa,
          double *rhs);
int msh(const iap_type *iap, const rap_type *rap, double *tm);
int genwts(const integer ncol, const integer n1, double *wt,
           double *wp);
int cpnts(const integer ncol, double *zm);
int cntdif(integer *n, double *d);
int wint(const integer n, double *wi);
int adptds(iap_type *iap, rap_type *rap, double *rds);
int adapt(iap_type *iap, rap_type *rap, integer *nold, integer *ncold,
          integer *nnew, integer *ncnew, double *tm, double *dtm,
          integer *ndxloc, double *ups, double *vps);
int interp(iap_type *iap, rap_type *rap, integer *ndim, integer *n, integer *nc,
           double *tm, integer *ndxloc, double *ups, integer *n1,
           integer *nc1, double *tm1, double *ups1, double *tm2,
           integer *itm1);
int newmsh(iap_type *iap, rap_type *rap, integer *ndxloc, double *ups,
           integer *nold, integer *ncold, double *tmold, double *dtmold,
           integer *nnew, double *tmnew, integer *iper);
int ordr(iap_type *iap, rap_type *rap, integer *n, double *tm, integer *n1,
         double *tm1, integer *itm1);
int intwts(iap_type *iap, rap_type *rap, integer *n, double *z__,
           double *x, double *wts);
int eqdf(iap_type *iap, rap_type *rap, integer *ntst, integer *ndim,
         integer *ncol, double *dtm, integer *ndxloc, double *ups,
         double *eqf, integer *iper);
int eig(iap_type *iap, integer *ndim, integer *m1a, double *a,
        doublecomplex *ev, integer *ier);
int nlvc(integer n, integer m, integer k, double *a, double *u);
int nrmlz(integer *ndim, double *v);
double pi(double r__);
int ge(integer n, integer m1a, double *a, integer nrhs, integer ndxloc,
       double *u, integer m1f, double *f, double *det);
int newlab(iap_type *iap, rap_type *rap);
int findlb(iap_type *iap, const rap_type *rap, integer irs, integer *nfpr,
           logical *found);
int readlb(const iap_type *iap, const rap_type *rap, double *u,
           double *par);
int skip3(integer *nskip, logical *eof3);
double rinpr(iap_type *iap, integer *ndim1, integer *ndxloc,
                 double *ups, double *vps, double *dtm,
                 double *thu);
double rnrmsq(iap_type *iap, integer *ndim1, integer *ndxloc,
                  double *ups, double *dtm, double *thu);
double rintg(iap_type *iap, integer *ndxloc, integer ic, double *ups,
                 double *dtm);
double rnrm2(iap_type *iap, integer *ndxloc, integer *ic, double *ups,
                 double *dtm);
double rmxups(iap_type *iap, integer *ndxloc, integer *i__,
                  double *ups);
double rmnups(iap_type *iap, integer *ndxloc, integer *i__,
                  double *ups);
int scaleb(iap_type *iap, integer *icp, integer *ndxloc, double *dvps,
           double *rld, double *dtm, double *thl, double *thu);
int cnrlbv(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           STPNT_TYPE_BVP((*stpnt)), PVLI_TYPE_BVP((*pvli)), double *thl,
           double *thu, integer *iuz, double *vuz);
int contbv(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), double *rds, double *rlcur,
           double *rlold, double *rldot, integer *ndxloc,
           double *ups, double *uoldps, double *udotps,
           double *upoldp, double *dtm, double *thl,
           double *thu);
int extrbv(iap_type *iap, rap_type *rap, FUNI_TYPE((*funi)), double *rds,
           double *rlcur, double *rlold, double *rldot,
           integer *ndxloc, double *ups, double *uoldps,
           double *udotps);
int stupbv(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), double *rlcur, double *rlold,
           double *rldot, integer *ndxloc, double *ups,
           double *uoldps, double *upoldp);
int stepbv(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           PVLI_TYPE_BVP((*pvli)), double *rds, double *rlcur,
           double *rlold, double *rldot, integer *ndxloc,
           double *ups, double *dups, double *uoldps,
           double *udotps, double *upoldp, double *fa,
           double *fc, double *tm, double *dtm, double *p0,
           double *p1, double *thl, double *thu);
int rsptbv(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), STPNT_TYPE_BVP((*stpnt)), double *rds,
           double *rlcur, double *rlold, double *rldot,
           integer *ndxloc, double *ups, double *uoldps,
           double *udotps, double *upoldp, double *dups,
           double *tm, double *dtm, doublecomplex *ev, integer *nodir,
           double *thl, double *thu);
STPNT_TYPE_BVP(stpnbv);
STPNT_TYPE_BVP(stpnub);
int setrtn(iap_type *iap, integer *ntst, integer *ndxloc, double *ups,
           double *par);
int stdrbv(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           double *rlcur, double *rlold, double *rldot,
           integer ndxloc, double *ups, double *dups,
           double *uoldps, double *udotps, double *upoldp,
           double *fa, double *fc, double *dtm, integer iperp,
           double *p0, double *p1, double *thl, double *thu);
int lcspbv(iap_type *iap, rap_type *rap, double *par, integer *icp,
           FNCS_TYPE_BVP((*fncs)), FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)),
           ICNI_TYPE((*icni)), PVLI_TYPE_BVP((*pvli)), double *q,
           double *rlcur, double *rlold, double *rldot,
           integer *ndxloc, double *ups, double *dups,
           double *uoldps, double *udotps, double *upoldp,
           double *fa, double *fc, double *tm, double *dtm,
           double *p0, double *p1, doublecomplex *ev, double *thl,
           double *thu, integer *iuz, double *vuz);
FNCS_TYPE_BVP(fnlpbv);
FNCS_TYPE_BVP(fnbpbv);
FNCS_TYPE_BVP(fnspbv);
FNCS_TYPE_BVP(fnuzbv);
int tpspbv(iap_type *iap, rap_type *rap, double *par, integer *icp,
           doublecomplex *ev);
int stplbv(iap_type *iap, rap_type *rap, double *par, integer *icp,
           double *rldot, integer *ndxloc, double *ups,
           double *udotps, double *tm, double *dtm, double *thl,
           double *thu);
int wrtbv8(iap_type *iap, rap_type *rap, double *par, integer *icp,
           double *rldot, integer *ndxloc, double *ups,
           double *udotps, double *tm, double *dtm);
int wrtbv9(iap_type *iap, rap_type *rap, double *par, integer *icp,
           double *rlcur, integer *ndxloc, double *ups, double *tm,
           double *dtm, double *thl, double *thu);
PVLI_TYPE_AE(pvlsae);
PVLI_TYPE_BVP(pvlsbv);
int setpae(iap_type *iap, rap_type *rap);
int setpbv(iap_type *iap, rap_type *rap, double *dtm);
int autim0(double *t);
int autim1(double *t);
double getp(char *code, integer *ic, double *ups, integer code_len);
/* autlib2.c */
int solvbv(integer *ifst, iap_type *iap, rap_type *rap, double *par,
           integer *icp, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)),
           ICNI_TYPE((*icni)), double *rds, integer *nllv,
           double *rlcur, double *rlold, double *rldot,
           integer *ndxloc, double *ups, double *dups,
           double *uoldps, double *udotps, double *upoldp,
           double *dtm, double *fa, double *fc, double *p0,
           double *p1, double *thl, double *thu);
int setfcdd(integer *ifst, double *dd, double *fc, integer *ncb,
            integer *nrc);
int faft(double *ff, double *fa, integer *ntst, integer *nrow,
         integer *ndxloc);
int partition(integer *n, integer *kwt, integer *m);
integer mypart(integer *iam, integer *np);
int setrhs(integer *ndim, integer *ips, integer *na, integer *ntst, integer *np,
           integer *ncol, integer *nbc, integer *nint, integer *ncb,
           integer *nrc, integer *nra, integer *nca, integer *iam, integer *kwt,
           logical *ipar, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)),
           ICNI_TYPE((*icni)), integer *ndxloc, iap_type *iap, rap_type *rap,
           double *par, integer *icp, double *rds, double *fa,
           double *fc, double *rlcur, double *rlold,
           double *rldot, double *ups, double *uoldps,
           double *udotps, double *upoldp, double *dups,
           double *dtm, double *thl, double *thu, double *p0,
           double *p1);
int brbd(double *a, double *b, double *c, double *d,
         double *fa, double *fc, double *p0, double *p1,
         integer *ifst, integer *idb, integer *nllv, double *det,
         integer *nov, integer *na, integer *nbc, integer *nra, integer *nca,
         integer *ncb, integer *nrc, integer *iam, integer *kwt, logical *par,
         double *a1, double *a2, double *bb, double *cc,
         double *faa, double *ca1, double *s1, double *s2,
         integer *icf11, integer *ipr, integer *icf1, integer *icf2,
         integer *irf, integer *icf);
int setzero(double *fa, double *fc, integer *na, integer *nra,
            integer *nrc);
int conrhs(integer *nov, integer *na, integer *nra, integer *nca, double *a,
           integer *nbc, integer *nrc, double *c, double *fa,
           double *fc, integer *irf, integer *icf, integer *iam);
int copycp(integer *iam, integer *kwt, integer *na, integer *nov, integer *nra,
           integer *nca, double *a, integer *ncb, double *b,
           integer *nrc, double *c, double *a1, double *a2,
           double *bb, double *cc, integer *irf);
int cpyrhs(integer *na, integer *nov, integer *nra, double *faa,
           double *fa, integer *irf);
int reduce(integer *iam, integer *kwt, logical *par, double *a1,
           double *a2, double *bb, double *cc, double *dd,
           integer *na, integer *nov, integer *ncb, integer *nrc,
           double *s1, double *s2, double *ca1, integer *icf1,
           integer *icf2, integer *icf11, integer *ipr, integer *nbc);
int redrhs(integer *iam, integer *kwt, logical *par, double *a1,
           double *a2, double *cc, double *faa, double *fc,
           integer *na, integer *nov, integer *ncb, integer *nrc,
           double *ca1, integer *icf1, integer *icf2, integer *icf11,
           integer *ipr, integer *nbc);
int dimrge(integer *iam, integer *kwt, logical *par, double *e,
           double *cc, double *d, double *fc, integer *ifst,
           integer *na, integer *nrc, integer *nov, integer *ncb, integer *idb,
           integer *nllv, double *fcc, double *p0, double *p1,
           double *det, double *s, double *a2, double *faa,
           double *bb);
int bcksub(integer *iam, integer *kwt, logical *par, double *s1,
           double *s2, double *a2, double *bb, double *faa,
           double *fc, double *fcc, double *sol1, double *sol2,
           double *sol3, integer *na, integer *nov, integer *ncb,
           integer *icf2);
int infpar(integer *iam, logical *par, double *a, double *b,
           double *fa, double *sol1, double *sol2, double *fc,
           integer *na, integer *nov, integer *nra, integer *nca, integer *ncb,
           integer *irf, integer *icf);
int rd0(integer *iam, integer *kwt, double *d, integer *nrc);
int print1(integer *nov, integer *na, integer *nra, integer *nca, integer *ncb,
           integer *nrc, double *a, double *b, double *c,
           double *d, double *fa, double *fc);
integer mynode(void);
integer numnodes(void);
int gsync(void);
double dclock(void);
int csend(void);
int crecv(void);
int gdsum(void);
int gsendx(void);
int gcol(void);
int led(void);
int setiomode(void);
/* autlib3.c */
FUNI_TYPE(fnlp);
int fflp(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, double *f, integer ndm, double *dfdu,
         double *dfdp);
STPNT_TYPE_AE(stpnlp);
FUNI_TYPE(fnc1);
STPNT_TYPE_AE(stpnc1);
FUNI_TYPE(fnc2);
int ffc2(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, double *f, integer ndm, double *dfdu,
         double *dfdp);
STPNT_TYPE_AE(stpnc2);
FUNI_TYPE(fnds);
FUNI_TYPE(fnti);
FUNI_TYPE(fnhd);
int ffhd(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, double *f, integer ndm, double *dfdu,
         double *dfdp);
STPNT_TYPE_AE(stpnhd);
FUNI_TYPE(fnhb);
int ffhb(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, double *f, integer ndm, double *dfdu,
         double *dfdp);
STPNT_TYPE_AE(stpnhb);
FUNI_TYPE(fnhw);
int ffhw(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, double *f, integer ndm, double *dfdu,
         double *dfdp);
STPNT_TYPE_AE(stpnhw);
FUNI_TYPE(fnps);
BCNI_TYPE(bcps);
ICNI_TYPE(icps);
int pdble(const iap_type *iap, const rap_type *rap, integer *ndim,
          integer *ntst, integer *ncol, integer *ndxloc, double *ups,
          double *udotps, double *tm, double *par);
STPNT_TYPE_BVP(stpnps);
FUNI_TYPE(fnws);
int ffws(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, integer ijac, double *f, double *dfdu,
         double *dfdp, integer ndm, double *dfu, double *dfp);
FUNI_TYPE(fnwp);
int stpnwp(iap_type *iap, rap_type *rap, double *par, integer *icp,
           integer *ntsr, integer *ncolrs, double *rlcur, double *rldot,
           integer *ndxloc, double *ups, double *udotps,
           double *upoldp, double *tm, double *dtm, integer *nodir,
           double *thl, double *thu);
FUNI_TYPE(fnsp);
int ffsp(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, integer ijac, double *f, double *dfdu,
         double *dfdp, integer ndm, double *dfu, double *dfp);
FUNI_TYPE(fnpe);
int ffpe(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, integer ijac, double *f, double *dfdu,
         double *dfdp, integer ndm, double *dfu, double *dfp);
ICNI_TYPE(icpe);
FUNI_TYPE(fnpl);
int ffpl(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, double *f, integer ndm, double *dfdu,
         double *dfdp);
BCNI_TYPE(bcpl);
ICNI_TYPE(icpl);
int stpnpl(iap_type *iap, rap_type *rap, double *par, integer *icp,
           integer *ntsr, integer *ncolrs, double *rlcur, double *rldot,
           integer *ndxloc, double *ups, double *udotps,
           double *upoldp, double *tm, double *dtm, integer *nodir,
           double *thl, double *thu);
FUNI_TYPE(fnpd);
int ffpd(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, double *f, integer ndm, double *dfdu,
         double *dfdp);
BCNI_TYPE(bcpd);
ICNI_TYPE(icpd);
int stpnpd(iap_type *iap, rap_type *rap, double *par, integer *icp,
           integer *ntsr, integer *ncolrs, double *rlcur, double *rldot,
           integer *ndxloc, double *ups, double *udotps,
           double *upoldp, double *tm, double *dtm, integer *nodir,
           double *thl, double *thu);
FUNI_TYPE(fntr);
int fftr(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, double *f, integer ndm, double *dfdu,
         double *dfdp);
int bctr(const iap_type *iap, const rap_type *rap, integer ndim,
         double *par, const integer *icp, integer nbc, const double *u0,
         const double *u1, double *f, integer ijac, double *dbc);
ICNI_TYPE(ictr);
int stpntr(iap_type *iap, rap_type *rap, double *par, integer *icp,
           integer *ntsr, integer *ncolrs, double *rlcur, double *rldot,
           integer *ndxloc, double *ups, double *udotps,
           double *upoldp, double *tm, double *dtm, integer *nodir,
           double *thl, double *thu);
FUNI_TYPE(fnpo);
int ffpo(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const double *upold,
         const integer *icp, double *par, double *f, integer ndm,
         double *dfdu, double *dfdp);
BCNI_TYPE(bcpo);
ICNI_TYPE(icpo);
int fipo(const iap_type *iap, const rap_type *rap, integer ndim,
         double *par, const integer *icp, integer nint, integer nnt0,
         const double *u, const double *uold, const double *udot,
         const double *upold, double *fi, double *dint,
         integer ndmt, double *dfdu, double *dfdp);
STPNT_TYPE_BVP(stpnpo);
FUNI_TYPE(fnbl);
int ffbl(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, double *f, integer ndm, double *dfdu,
         double *dfdp);
BCNI_TYPE(bcbl);
int fbbl(const iap_type *iap, const rap_type *rap, integer ndim,
         double *par, const integer *icp, integer nbc, integer nbc0,
         const double *u0, const double *u1, double *f,
         double *dbc);
ICNI_TYPE(icbl);
int fibl(const iap_type *iap, const rap_type *rap, integer ndim,
         double *par, const integer *icp, integer nint, integer nnt0,
         const double *u, const double *uold, const double *udot,
         const double *upold, double *f, double *dint);
STPNT_TYPE_BVP(stpnbl);
FUNI_TYPE(funi);
BCNI_TYPE(bcni);
ICNI_TYPE(icni);
int fopi(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const integer *icp, double *par, integer ijac,
         double *f, double *dfdu, double *dfdp);
/* autlib4.c */
int flowkm(integer *ndim, double *c0, double *c1, integer *iid,
           double *rwork, doublecomplex *ev);
int dhhpr(integer *k, integer *j, integer *n, double *x, integer *incx,
          double *beta, double *v);
int dhhap(integer *k, integer *j, integer *n, integer *q, double *beta,
          double *v, integer *job, double *a, integer *lda);
/* autlib5.c */
FUNI_TYPE(fnho);
int ffho(const iap_type *iap, const rap_type *rap, integer ndim,
         const double *u, const double *uold, const integer *icp,
         double *par, double *f, integer ndm, double *dfdu,
         double *dfdp);
BCNI_TYPE(bcho);
int fbho(const iap_type *iap, const rap_type *rap, integer ndim,
         double *par, const integer *icp, integer nbc, integer nbc0,
         const double *u0, const double *u1, double *fb,
         double *dbc);
ICNI_TYPE(icho);
int fiho(const iap_type *iap, const rap_type *rap, integer ndim,
         double *par, const integer *icp, integer nint, integer nnt0,
         const double *u, const double *uold, const double *udot,
         const double *upold, double *fi, double *dint);
int inho(iap_type *iap, integer *icp, double *par);
int preho(integer *ndx, integer *ntsr, integer *nar, integer *ndim,
          integer *ncolrs, double *ups, double *udotps, double *tm,
          double *par);
int stpnho(iap_type *iap, rap_type *rap, double *par, integer *icp,
           integer *ntsr, integer *ncolrs, double *rlcur, double *rldot,
           integer *ndxloc, double *ups, double *udotps,
           double *upoldp, double *tm, double *dtm, integer *nodir,
           double *thl, double *thu);
int stpho(iap_type *iap, integer *icp, double *u, double *par,
          double *t);
PVLI_TYPE_BVP(pvlsho);
double psiho(const iap_type *iap, integer is, double *rr,
                 double *ri, double *v, double *vt,
                 const integer *icp, double *par);
int eighi(integer isign, integer itrans, double *rr, double *ri,
          double *vret, double *xequib, const integer *icp,
          double *par, integer *ndm);
int eigho(integer *isign, integer *itrans, double *rr, double *ri,
          double *vret, double *xequib, const integer *icp,
          double *par, integer *ndm, double *dfdu, double *dfdp,
          double *zz);
int prjcti(double *bound, double *xequib, const integer *icp,
           double *par, integer imfd, integer is, integer itrans,
           integer *ndm);
int prjctn(double *bound, double *xequib, const integer *icp,
           double *par, integer *imfd, integer *is, integer *itrans,
           integer *ndm, double *dfdu, double *dfdp);
/* eispack.c */
int rg(integer nm, integer n, double *a, double *wr, double *wi,
       integer matz, double *z__, integer *iv1, double *fv1,
       integer *ierr);
int hqr(integer *nm, integer *n, integer *low, integer *igh, double *h__,
        double *wr, double *wi, integer *ierr);
int hqr2(integer *nm, integer *n, integer *low, integer *igh, double *h__,
         double *wr, double *wi, double *z__, integer *ierr);
int cdiv(double *ar, double *ai, double *br, double *bi,
         double *cr, double *ci);
int balanc(integer *nm, integer *n, double *a, integer *low, integer *igh,
           double *scale);
int balbak(integer *nm, integer *n, integer *low, integer *igh,
           double *scale, integer *m, double *z__);
int elmhes(integer *nm, integer *n, integer *low, integer *igh, double *a,
           integer *int__);
int eltran(integer *nm, integer *n, integer *low, integer *igh, double *a,
           integer *int__, double *z__);
int qzhes(integer nm, integer n, double *a, double *b, logical matz,
          double *z__);
int qzit(integer nm, integer n, double *a, double *b, double eps1,
         logical matz, double *z__, integer *ierr);
int qzval(integer nm, integer n, double *a, double *b, double *alfr,
          double *alfi, double *beta, logical matz, double *z__);
double epslon(double x);
double dnrm2(integer *n, double *dx, integer *incx);
double ddot(integer *n, double *dx, integer *incx, double *dy,
                integer *incy);
int dscal(integer *n, double *da, double *dx, integer *incx);
integer idamax(integer *n, double *dx, integer *incx);
int daxpy(integer *n, double *da, double *dx, integer *incx,
          double *dy, integer *incy);
int drot(integer *n, double *dx, integer *incx, double *dy,
         integer *incy, double *c, double *s);
int dswap(integer *n, double *dx, integer *incx, double *dy,
          integer *incy);
int dgemc(integer *m, integer *n, double *a, integer *lda, double *b,
          integer *ldb, logical *trans);
int xerbla(char *srname, integer *info, integer srname_len);
logical lsame(char *ca, char *cb, integer ca_len, integer cb_len);
int dgemm(char *transa, char *transb, integer *m, integer *n, integer *k,
          double *alpha, double *a, integer *lda, double *b,
          integer *ldb, double *beta, double *c, integer *ldc,
          integer transa_len, integer transb_len);
int ezsvd(double *x, integer *ldx, integer *n, integer *p, double *s,
          double *e, double *u, integer *ldu, double *v,
          integer *ldv, double *work, integer *job, integer *info,
          double *tol);
int ndrotg(double *f, double *g, double *cs, double *sn);
int ndsvd(double *x, integer *ldx, integer *n, integer *p, double *s,
          double *e, double *u, integer *ldu, double *v,
          integer *ldv, double *work, integer *job, integer *info,
          integer *maxitr, double *tol, integer *idbg, integer *ifull,
          integer *kount, integer *kount1, integer *kount2, integer *skip,
          integer *limshf, double *maxsin, integer *iidir);
int prse(integer *ll, integer *m, integer *nrow, integer *ncol, double *s,
         double *e);
int sig22(double *a, double *b, double *c, double *sigmin,
          double *sigmax, double *snr, double *csr, double *snl,
          double *csl);
double sigmin(double *a, double *b, double *c);
int sndrtg(double *f, double *g, double *cs, double *sn);
int hqr3lc(double *a, double *v, integer *n, integer *nlow,
           integer *nup, double *eps, double *er, double *ei,
           integer *type__, integer *na, integer *nv, integer *imfd);
int split(double *a, double *v, integer *n, integer *l, double *e1,
          double *e2, integer *na, integer *nv);
int exchng(double *a, double *v, integer *n, integer *l, integer *b1,
           integer *b2, double *eps, logical *fail, integer *na,
           integer *nv);
int qrstep(double *a, double *v, double *p, double *q,
           double *r__, integer *nl, integer *nu, integer *n, integer *na,
           integer *nv);
int orthes(integer *nm, integer *n, integer *low, integer *igh, double *a,
           double *ort);
int ortran(integer *nm, integer *n, integer *low, integer *igh, double *a,
           double *ort, double *z__);

/* problem defined functions*/
int func(integer ndim, const double *u, const integer *icp,
         const double *par, integer ijac, double *f, double *dfdu,
         double *dfdp);
int stpnt(integer ndim, double t, double *u, double *par);
int bcnd(integer ndim, const double *par, const integer *icp, integer nbc,
         const double *u0, const double *u1, integer ijac,
         double *f, double *dbc);
int icnd(integer ndim, const double *par, const integer *icp, integer nint,
         const double *u, const double *uold, const double *udot,
         const double *upold, integer ijac, double *fi,
         double *dint);
int fopt(integer ndim, const double *u, const integer *icp,
         const double *par, integer ijac, double *fs, double *dfdu,
         double *dfdp);
int pvls(integer ndim, const double *u, double *par);
/* conpar.c */
void *conpar_process(void *);
int conpar(integer *nov, integer *na, integer *nra, integer *nca, double *a,
           integer *ncb, double *b, integer *nbc, integer *nrc,
           double *c, double *d, integer *irf, integer *icf);
/*setubv.c */
#include "auto_types.h"
void *setubv_make_aa_bb_cc(void *);
int setubv(integer ndim, integer ips, integer na, integer ncol, integer nbc,
           integer nint, integer ncb, integer nrc, integer nra, integer nca,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           integer ndxloc, iap_type *iap, rap_type *rap, double *par,
           integer *icp, double rds, double *aa, double *bb,
           double *cc, double *dd, double *fa, double *fc,
           double *rlcur, double *rlold, double *rldot,
           double *ups, double *uoldps, double *udotps,
           double *upoldp, double *dups, double *dtm,
           double *thl, double *thu, double *p0, double *p1);
void setubv_parallel_arglist_copy(setubv_parallel_arglist *output,
                                  const setubv_parallel_arglist input);
void setubv_parallel_arglist_constructor(
    integer ndim, integer ips, integer na, integer ncol, integer nbc,
    integer nint, integer ncb, integer nrc, integer nra, integer nca,
    FUNI_TYPE((*funi)), ICNI_TYPE((*icni)), integer ndxloc, iap_type *iap,
    rap_type *rap, double *par, integer *icp, double *aa,
    double *bb, double *cc, double *dd, double *fa,
    double *fc, double *ups, double *uoldps, double *udotps,
    double *upoldp, double *dtm, double *wp, double *wt,
    double *wi, double *thu, double *thl, double *rldot,
    BCNI_TYPE((*bcni)), setubv_parallel_arglist *data);
void setubv_make_fa(setubv_parallel_arglist larg);
void setubv_make_fc_dd(setubv_parallel_arglist larg, double *dups,
                       double *rlcur, double *rlold, double rds);

/*worker.c*/
int mpi_worker();
int mpi_setubv_worker();
int mpi_conpar_worker();
#include "auto_types.h"
int set_funi_and_icni(iap_type *, setubv_parallel_arglist *);
int set_function_pointers(const iap_type, function_list *);

#ifdef AUTO_CONSTRUCT_DESCTRUCT
int user_construct(int argc, char **argv);
int user_destruct();
#endif

#endif
