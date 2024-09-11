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
    int X(const iap_type *iap, const rap_type *rap, int64 ndim,                \
          const double *u, const double *uold, const int64 *icp, double *par,  \
          int64 ijac, double *f, double *dfdu, double *dfdp)

/*This is the type for all functions which can be used as "bcni" the function
  which evaluates the boundary conditions */
#define BCNI_TYPE(X)                                                           \
    int X(const iap_type *iap, const rap_type *rap, int64 ndim, double *par,   \
          const int64 *icp, int64 nbc, const double *u0, const double *u1,     \
          double *f, int64 ijac, double *dbc)

/*This is the type for all functions which can be used as "icni" the function
  which evaluates kernel of the integral constraints */
#define ICNI_TYPE(X)                                                           \
    int X(const iap_type *iap, const rap_type *rap, int64 ndim, double *par,   \
          const int64 *icp, int64 nint, const double *u, const double *uold,   \
          const double *udot, const double *upold, double *f, int64 ijac,      \
          double *dint)

/*This is the type for all functions which can be used as additional
  output functions for algebraic problems */
#define PVLI_TYPE_AE(X)                                                        \
    int X(iap_type *iap, rap_type *rap, double *u, double *par)

/*This is the type for all functions which can be used as additional
  output functions for BVPs */
#define PVLI_TYPE_BVP(X)                                                       \
    int X(iap_type *iap, rap_type *rap, int64 *icp, double *dtm,               \
          int64 *ndxloc, double *ups, int64 *ndim, double *p0, double *p1,     \
          double *par)

/* This is the type for all functions that can be used at starting points
   for algebraic problems */
#define STPNT_TYPE_AE(X)                                                       \
    int X(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *u)

/* This is the type for all functions that can be used at starting points
   for BVPs */
#define STPNT_TYPE_BVP(X)                                                      \
    int X(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsrs, \
          int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,          \
          double *ups, double *udotps, double *upoldp, double *tm,             \
          double *dtm, int64 *nodir, double *thl, double *thu)

/*This is the type for all functions which can be used to detect
  special points for algebraic problems */
#define FNCS_TYPE_AE(X)                                                        \
    double X(iap_type *iap, rap_type *rap, double *par, int64 *icp,            \
             logical *chng, FUNI_TYPE((*funi)), int64 *m1aaloc, double *aa,    \
             double *rlcur, double *rlold, double *rldot, double *u,           \
             double *uold, double *udot, double *rhs, double *dfdu,            \
             double *dfdp, int64 *iuz, double *vuz)

/*This is the type for all functions which can be used to detect
  special points for BVPS */
#define FNCS_TYPE_BVP(X)                                                       \
    double X(iap_type *iap, rap_type *rap, double *par, int64 *icp,            \
             logical *chng, FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)),            \
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
int init(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *thl,
         double **thu_pointer, int64 *iuz, double *vuz);
int chdim(iap_type *iap);
int autoae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), STPNT_TYPE_AE((*stpnt)), PVLI_TYPE_AE((*pvli)),
           double *thl, double *thu, int64 *iuz, double *vuz);

int autobv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           STPNT_TYPE_BVP((*stpnt)), PVLI_TYPE_BVP((*pvli)), double *thl,
           double *thu, int64 *iuz, double *vuz);
int init1(iap_type *iap, rap_type *rap, int64 *icp, double *par);
int cnrlae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), STPNT_TYPE_AE((*stpnt)), PVLI_TYPE_AE((*pvli)),
           double *thl, double *thu, int64 *iuz, double *vuz);
STPNT_TYPE_AE(stpnus);
STPNT_TYPE_AE(stpnae);
int stprae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), double *rds, int64 *m1aaloc, double *aa,
           double *rhs, double *rlcur, double *rlold, double *rldot, double *u,
           double *du, double *uold, double *udot, double *f, double *dfdu,
           double *dfdp, double *thl, double *thu);
int contae(iap_type *iap, rap_type *rap, double *rds, double *rlcur,
           double *rlold, double *rldot, double *u, double *uold, double *udot);
int solvae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), double *rds, int64 *m1aaloc, double *aa,
           double *rhs, double *rlcur, double *rlold, double *rldot, double *u,
           double *du, double *uold, double *udot, double *f, double *dfdu,
           double *dfdp, double *thl, double *thu);
int lcspae(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FNCS_TYPE_AE((*fncs)), FUNI_TYPE((*funi)), int64 *m1aaloc,
           double *aa, double *rhs, double *rlcur, double *rlold, double *rldot,
           double *u, double *du, double *uold, double *udot, double *f,
           double *dfdu, double *dfdp, double *q, double *thl, double *thu,
           int64 *iuz, double *vuz);
int mueller(double *q0, double *q1, double *q, double *s0, double *s1,
            double *s, double *rds);
FNCS_TYPE_AE(fnbpae);
FNCS_TYPE_AE(fnlpae);
FNCS_TYPE_AE(fnhbae);
FNCS_TYPE_AE(fnuzae);
int stbif(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *m1aaloc,
          double *aa, int64 m1sbloc, double *stud, double *stu, double *stla,
          double *stld, double *rlcur, double *rlold, double *rldot, double *u,
          double *du, double *udot, double *dfdu, double *dfdp, double *thl,
          double *thu);
int swpnt(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *rds,
          int64 m1sbloc, double *stud, double *stu, double *stla, double *stld,
          double *rlcur, double *rlold, double *rldot, double *u, double *udot);
int swprc(iap_type *iap, rap_type *rap, double *par, int64 *icp,
          FUNI_TYPE((*funi)), int64 *m1aaloc, double *aa, double *rhs,
          double *rlcur, double *rlold, double *rldot, double *u, double *du,
          double *uold, double *udot, double *f, double *dfdu, double *dfdp,
          double *rds, double *thl, double *thu);
int sthd(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *thl,
         double *thu);
int headng(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 iunit,
           int64 *n1, int64 *n2);
int stplae(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *rlcur,
           double *u);
int wrline(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *icu,
           int64 *ibr, int64 *ntot, int64 *lab, double *vaxis, double *u);
int wrtsp8(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *lab,
           double *rlcur, double *u);
int wrjac(iap_type *iap, int64 *n, int64 *m1aaloc, double *aa, double *rhs);
int msh(const iap_type *iap, const rap_type *rap, double *tm);
int genwts(const int64 ncol, const int64 n1, double *wt, double *wp);
int cpnts(const int64 ncol, double *zm);
int cntdif(int64 *n, double *d);
int wint(const int64 n, double *wi);
int adptds(iap_type *iap, rap_type *rap, double *rds);
int adapt(iap_type *iap, rap_type *rap, int64 *nold, int64 *ncold, int64 *nnew,
          int64 *ncnew, double *tm, double *dtm, int64 *ndxloc, double *ups,
          double *vps);
int interp(iap_type *iap, rap_type *rap, int64 *ndim, int64 *n, int64 *nc,
           double *tm, int64 *ndxloc, double *ups, int64 *n1, int64 *nc1,
           double *tm1, double *ups1, double *tm2, int64 *itm1);
int newmsh(iap_type *iap, rap_type *rap, int64 *ndxloc, double *ups,
           int64 *nold, int64 *ncold, double *tmold, double *dtmold,
           int64 *nnew, double *tmnew, int64 *iper);
int ordr(iap_type *iap, rap_type *rap, int64 *n, double *tm, int64 *n1,
         double *tm1, int64 *itm1);
int intwts(iap_type *iap, rap_type *rap, int64 *n, double *z__, double *x,
           double *wts);
int eqdf(iap_type *iap, rap_type *rap, int64 *ntst, int64 *ndim, int64 *ncol,
         double *dtm, int64 *ndxloc, double *ups, double *eqf, int64 *iper);
int eig(iap_type *iap, int64 *ndim, int64 *m1a, double *a, doublecomplex *ev,
        int64 *ier);
int nlvc(int64 n, int64 m, int64 k, double *a, double *u);
int nrmlz(int64 *ndim, double *v);
double pi(double r__);
int ge(int64 n, int64 m1a, double *a, int64 nrhs, int64 ndxloc, double *u,
       int64 m1f, double *f, double *det);
int newlab(iap_type *iap, rap_type *rap);
int findlb(iap_type *iap, const rap_type *rap, int64 irs, int64 *nfpr,
           logical *found);
int readlb(const iap_type *iap, const rap_type *rap, double *u, double *par);
int skip3(int64 *nskip, logical *eof3);
double rinpr(iap_type *iap, int64 *ndim1, int64 *ndxloc, double *ups,
             double *vps, double *dtm, double *thu);
double rnrmsq(iap_type *iap, int64 *ndim1, int64 *ndxloc, double *ups,
              double *dtm, double *thu);
double rintg(iap_type *iap, int64 *ndxloc, int64 ic, double *ups, double *dtm);
double rnrm2(iap_type *iap, int64 *ndxloc, int64 *ic, double *ups, double *dtm);
double rmxups(iap_type *iap, int64 *ndxloc, int64 *i__, double *ups);
double rmnups(iap_type *iap, int64 *ndxloc, int64 *i__, double *ups);
int scaleb(iap_type *iap, int64 *icp, int64 *ndxloc, double *dvps, double *rld,
           double *dtm, double *thl, double *thu);
int cnrlbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           STPNT_TYPE_BVP((*stpnt)), PVLI_TYPE_BVP((*pvli)), double *thl,
           double *thu, int64 *iuz, double *vuz);
int contbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), double *rds, double *rlcur, double *rlold,
           double *rldot, int64 *ndxloc, double *ups, double *uoldps,
           double *udotps, double *upoldp, double *dtm, double *thl,
           double *thu);
int extrbv(iap_type *iap, rap_type *rap, FUNI_TYPE((*funi)), double *rds,
           double *rlcur, double *rlold, double *rldot, int64 *ndxloc,
           double *ups, double *uoldps, double *udotps);
int stupbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), double *rlcur, double *rlold, double *rldot,
           int64 *ndxloc, double *ups, double *uoldps, double *upoldp);
int stepbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           PVLI_TYPE_BVP((*pvli)), double *rds, double *rlcur, double *rlold,
           double *rldot, int64 *ndxloc, double *ups, double *dups,
           double *uoldps, double *udotps, double *upoldp, double *fa,
           double *fc, double *tm, double *dtm, double *p0, double *p1,
           double *thl, double *thu);
int rsptbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), STPNT_TYPE_BVP((*stpnt)), double *rds,
           double *rlcur, double *rlold, double *rldot, int64 *ndxloc,
           double *ups, double *uoldps, double *udotps, double *upoldp,
           double *dups, double *tm, double *dtm, doublecomplex *ev,
           int64 *nodir, double *thl, double *thu);
STPNT_TYPE_BVP(stpnbv);
STPNT_TYPE_BVP(stpnub);
int setrtn(iap_type *iap, int64 *ntst, int64 *ndxloc, double *ups, double *par);
int stdrbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           double *rlcur, double *rlold, double *rldot, int64 ndxloc,
           double *ups, double *dups, double *uoldps, double *udotps,
           double *upoldp, double *fa, double *fc, double *dtm, int64 iperp,
           double *p0, double *p1, double *thl, double *thu);
int lcspbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FNCS_TYPE_BVP((*fncs)), FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)),
           ICNI_TYPE((*icni)), PVLI_TYPE_BVP((*pvli)), double *q, double *rlcur,
           double *rlold, double *rldot, int64 *ndxloc, double *ups,
           double *dups, double *uoldps, double *udotps, double *upoldp,
           double *fa, double *fc, double *tm, double *dtm, double *p0,
           double *p1, doublecomplex *ev, double *thl, double *thu, int64 *iuz,
           double *vuz);
FNCS_TYPE_BVP(fnlpbv);
FNCS_TYPE_BVP(fnbpbv);
FNCS_TYPE_BVP(fnspbv);
FNCS_TYPE_BVP(fnuzbv);
int tpspbv(iap_type *iap, rap_type *rap, double *par, int64 *icp,
           doublecomplex *ev);
int stplbv(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *rldot,
           int64 *ndxloc, double *ups, double *udotps, double *tm, double *dtm,
           double *thl, double *thu);
int wrtbv8(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *rldot,
           int64 *ndxloc, double *ups, double *udotps, double *tm, double *dtm);
int wrtbv9(iap_type *iap, rap_type *rap, double *par, int64 *icp, double *rlcur,
           int64 *ndxloc, double *ups, double *tm, double *dtm, double *thl,
           double *thu);
PVLI_TYPE_AE(pvlsae);
PVLI_TYPE_BVP(pvlsbv);
int setpae(iap_type *iap, rap_type *rap);
int setpbv(iap_type *iap, rap_type *rap, double *dtm);
int autim0(double *t);
int autim1(double *t);
double getp(char *code, int64 *ic, double *ups, int64 code_len);
/* autlib2.c */
int solvbv(int64 *ifst, iap_type *iap, rap_type *rap, double *par, int64 *icp,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           double *rds, int64 *nllv, double *rlcur, double *rlold,
           double *rldot, int64 *ndxloc, double *ups, double *dups,
           double *uoldps, double *udotps, double *upoldp, double *dtm,
           double *fa, double *fc, double *p0, double *p1, double *thl,
           double *thu);
int setfcdd(int64 *ifst, double *dd, double *fc, int64 *ncb, int64 *nrc);
int faft(double *ff, double *fa, int64 *ntst, int64 *nrow, int64 *ndxloc);
int partition(int64 *n, int64 *kwt, int64 *m);
int64 mypart(int64 *iam, int64 *np);
int setrhs(int64 *ndim, int64 *ips, int64 *na, int64 *ntst, int64 *np,
           int64 *ncol, int64 *nbc, int64 *nint, int64 *ncb, int64 *nrc,
           int64 *nra, int64 *nca, int64 *iam, int64 *kwt, logical *ipar,
           FUNI_TYPE((*funi)), BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)),
           int64 *ndxloc, iap_type *iap, rap_type *rap, double *par, int64 *icp,
           double *rds, double *fa, double *fc, double *rlcur, double *rlold,
           double *rldot, double *ups, double *uoldps, double *udotps,
           double *upoldp, double *dups, double *dtm, double *thl, double *thu,
           double *p0, double *p1);
int brbd(double *a, double *b, double *c, double *d, double *fa, double *fc,
         double *p0, double *p1, int64 *ifst, int64 *idb, int64 *nllv,
         double *det, int64 *nov, int64 *na, int64 *nbc, int64 *nra, int64 *nca,
         int64 *ncb, int64 *nrc, int64 *iam, int64 *kwt, logical *par,
         double *a1, double *a2, double *bb, double *cc, double *faa,
         double *ca1, double *s1, double *s2, int64 *icf11, int64 *ipr,
         int64 *icf1, int64 *icf2, int64 *irf, int64 *icf);
int setzero(double *fa, double *fc, int64 *na, int64 *nra, int64 *nrc);
int conrhs(int64 *nov, int64 *na, int64 *nra, int64 *nca, double *a, int64 *nbc,
           int64 *nrc, double *c, double *fa, double *fc, int64 *irf,
           int64 *icf, int64 *iam);
int copycp(int64 *iam, int64 *kwt, int64 *na, int64 *nov, int64 *nra,
           int64 *nca, double *a, int64 *ncb, double *b, int64 *nrc, double *c,
           double *a1, double *a2, double *bb, double *cc, int64 *irf);
int cpyrhs(int64 *na, int64 *nov, int64 *nra, double *faa, double *fa,
           int64 *irf);
int reduce(int64 *iam, int64 *kwt, logical *par, double *a1, double *a2,
           double *bb, double *cc, double *dd, int64 *na, int64 *nov,
           int64 *ncb, int64 *nrc, double *s1, double *s2, double *ca1,
           int64 *icf1, int64 *icf2, int64 *icf11, int64 *ipr, int64 *nbc);
int redrhs(int64 *iam, int64 *kwt, logical *par, double *a1, double *a2,
           double *cc, double *faa, double *fc, int64 *na, int64 *nov,
           int64 *ncb, int64 *nrc, double *ca1, int64 *icf1, int64 *icf2,
           int64 *icf11, int64 *ipr, int64 *nbc);
int dimrge(int64 *iam, int64 *kwt, logical *par, double *e, double *cc,
           double *d, double *fc, int64 *ifst, int64 *na, int64 *nrc,
           int64 *nov, int64 *ncb, int64 *idb, int64 *nllv, double *fcc,
           double *p0, double *p1, double *det, double *s, double *a2,
           double *faa, double *bb);
int bcksub(int64 *iam, int64 *kwt, logical *par, double *s1, double *s2,
           double *a2, double *bb, double *faa, double *fc, double *fcc,
           double *sol1, double *sol2, double *sol3, int64 *na, int64 *nov,
           int64 *ncb, int64 *icf2);
int infpar(int64 *iam, logical *par, double *a, double *b, double *fa,
           double *sol1, double *sol2, double *fc, int64 *na, int64 *nov,
           int64 *nra, int64 *nca, int64 *ncb, int64 *irf, int64 *icf);
int rd0(int64 *iam, int64 *kwt, double *d, int64 *nrc);
int print1(int64 *nov, int64 *na, int64 *nra, int64 *nca, int64 *ncb,
           int64 *nrc, double *a, double *b, double *c, double *d, double *fa,
           double *fc);
int64 mynode(void);
int64 numnodes(void);
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
int fflp(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, double *f,
         int64 ndm, double *dfdu, double *dfdp);
STPNT_TYPE_AE(stpnlp);
FUNI_TYPE(fnc1);
STPNT_TYPE_AE(stpnc1);
FUNI_TYPE(fnc2);
int ffc2(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, double *f,
         int64 ndm, double *dfdu, double *dfdp);
STPNT_TYPE_AE(stpnc2);
FUNI_TYPE(fnds);
FUNI_TYPE(fnti);
FUNI_TYPE(fnhd);
int ffhd(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, double *f,
         int64 ndm, double *dfdu, double *dfdp);
STPNT_TYPE_AE(stpnhd);
FUNI_TYPE(fnhb);
int ffhb(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, double *f,
         int64 ndm, double *dfdu, double *dfdp);
STPNT_TYPE_AE(stpnhb);
FUNI_TYPE(fnhw);
int ffhw(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, double *f,
         int64 ndm, double *dfdu, double *dfdp);
STPNT_TYPE_AE(stpnhw);
FUNI_TYPE(fnps);
BCNI_TYPE(bcps);
ICNI_TYPE(icps);
int pdble(const iap_type *iap, const rap_type *rap, int64 *ndim, int64 *ntst,
          int64 *ncol, int64 *ndxloc, double *ups, double *udotps, double *tm,
          double *par);
STPNT_TYPE_BVP(stpnps);
FUNI_TYPE(fnws);
int ffws(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, int64 ijac,
         double *f, double *dfdu, double *dfdp, int64 ndm, double *dfu,
         double *dfp);
FUNI_TYPE(fnwp);
int stpnwp(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr,
           int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,
           double *ups, double *udotps, double *upoldp, double *tm, double *dtm,
           int64 *nodir, double *thl, double *thu);
FUNI_TYPE(fnsp);
int ffsp(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, int64 ijac,
         double *f, double *dfdu, double *dfdp, int64 ndm, double *dfu,
         double *dfp);
FUNI_TYPE(fnpe);
int ffpe(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, int64 ijac,
         double *f, double *dfdu, double *dfdp, int64 ndm, double *dfu,
         double *dfp);
ICNI_TYPE(icpe);
FUNI_TYPE(fnpl);
int ffpl(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, double *f,
         int64 ndm, double *dfdu, double *dfdp);
BCNI_TYPE(bcpl);
ICNI_TYPE(icpl);
int stpnpl(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr,
           int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,
           double *ups, double *udotps, double *upoldp, double *tm, double *dtm,
           int64 *nodir, double *thl, double *thu);
FUNI_TYPE(fnpd);
int ffpd(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, double *f,
         int64 ndm, double *dfdu, double *dfdp);
BCNI_TYPE(bcpd);
ICNI_TYPE(icpd);
int stpnpd(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr,
           int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,
           double *ups, double *udotps, double *upoldp, double *tm, double *dtm,
           int64 *nodir, double *thl, double *thu);
FUNI_TYPE(fntr);
int fftr(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, double *f,
         int64 ndm, double *dfdu, double *dfdp);
int bctr(const iap_type *iap, const rap_type *rap, int64 ndim, double *par,
         const int64 *icp, int64 nbc, const double *u0, const double *u1,
         double *f, int64 ijac, double *dbc);
ICNI_TYPE(ictr);
int stpntr(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr,
           int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,
           double *ups, double *udotps, double *upoldp, double *tm, double *dtm,
           int64 *nodir, double *thl, double *thu);
FUNI_TYPE(fnpo);
int ffpo(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const double *upold, const int64 *icp, double *par,
         double *f, int64 ndm, double *dfdu, double *dfdp);
BCNI_TYPE(bcpo);
ICNI_TYPE(icpo);
int fipo(const iap_type *iap, const rap_type *rap, int64 ndim, double *par,
         const int64 *icp, int64 nint, int64 nnt0, const double *u,
         const double *uold, const double *udot, const double *upold,
         double *fi, double *dint, int64 ndmt, double *dfdu, double *dfdp);
STPNT_TYPE_BVP(stpnpo);
FUNI_TYPE(fnbl);
int ffbl(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, double *f,
         int64 ndm, double *dfdu, double *dfdp);
BCNI_TYPE(bcbl);
int fbbl(const iap_type *iap, const rap_type *rap, int64 ndim, double *par,
         const int64 *icp, int64 nbc, int64 nbc0, const double *u0,
         const double *u1, double *f, double *dbc);
ICNI_TYPE(icbl);
int fibl(const iap_type *iap, const rap_type *rap, int64 ndim, double *par,
         const int64 *icp, int64 nint, int64 nnt0, const double *u,
         const double *uold, const double *udot, const double *upold, double *f,
         double *dint);
STPNT_TYPE_BVP(stpnbl);
FUNI_TYPE(funi);
BCNI_TYPE(bcni);
ICNI_TYPE(icni);
int fopi(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const int64 *icp, double *par, int64 ijac, double *f, double *dfdu,
         double *dfdp);
/* autlib4.c */
int flowkm(int64 *ndim, double *c0, double *c1, int64 *iid, double *rwork,
           doublecomplex *ev);
int dhhpr(int64 *k, int64 *j, int64 *n, double *x, int64 *incx, double *beta,
          double *v);
int dhhap(int64 *k, int64 *j, int64 *n, int64 *q, double *beta, double *v,
          int64 *job, double *a, int64 *lda);
/* autlib5.c */
FUNI_TYPE(fnho);
int ffho(const iap_type *iap, const rap_type *rap, int64 ndim, const double *u,
         const double *uold, const int64 *icp, double *par, double *f,
         int64 ndm, double *dfdu, double *dfdp);
BCNI_TYPE(bcho);
int fbho(const iap_type *iap, const rap_type *rap, int64 ndim, double *par,
         const int64 *icp, int64 nbc, int64 nbc0, const double *u0,
         const double *u1, double *fb, double *dbc);
ICNI_TYPE(icho);
int fiho(const iap_type *iap, const rap_type *rap, int64 ndim, double *par,
         const int64 *icp, int64 nint, int64 nnt0, const double *u,
         const double *uold, const double *udot, const double *upold,
         double *fi, double *dint);
int inho(iap_type *iap, int64 *icp, double *par);
int preho(int64 *ndx, int64 *ntsr, int64 *nar, int64 *ndim, int64 *ncolrs,
          double *ups, double *udotps, double *tm, double *par);
int stpnho(iap_type *iap, rap_type *rap, double *par, int64 *icp, int64 *ntsr,
           int64 *ncolrs, double *rlcur, double *rldot, int64 *ndxloc,
           double *ups, double *udotps, double *upoldp, double *tm, double *dtm,
           int64 *nodir, double *thl, double *thu);
int stpho(iap_type *iap, int64 *icp, double *u, double *par, double *t);
PVLI_TYPE_BVP(pvlsho);
double psiho(const iap_type *iap, int64 is, double *rr, double *ri, double *v,
             double *vt, const int64 *icp, double *par);
int eighi(int64 isign, int64 itrans, double *rr, double *ri, double *vret,
          double *xequib, const int64 *icp, double *par, int64 *ndm);
int eigho(int64 *isign, int64 *itrans, double *rr, double *ri, double *vret,
          double *xequib, const int64 *icp, double *par, int64 *ndm,
          double *dfdu, double *dfdp, double *zz);
int prjcti(double *bound, double *xequib, const int64 *icp, double *par,
           int64 imfd, int64 is, int64 itrans, int64 *ndm);
int prjctn(double *bound, double *xequib, const int64 *icp, double *par,
           int64 *imfd, int64 *is, int64 *itrans, int64 *ndm, double *dfdu,
           double *dfdp);
/* eispack.c */
int rg(int64 nm, int64 n, double *a, double *wr, double *wi, int64 matz,
       double *z__, int64 *iv1, double *fv1, int64 *ierr);
int hqr(int64 *nm, int64 *n, int64 *low, int64 *igh, double *h__, double *wr,
        double *wi, int64 *ierr);
int hqr2(int64 *nm, int64 *n, int64 *low, int64 *igh, double *h__, double *wr,
         double *wi, double *z__, int64 *ierr);
int cdiv(double *ar, double *ai, double *br, double *bi, double *cr,
         double *ci);
int balanc(int64 *nm, int64 *n, double *a, int64 *low, int64 *igh,
           double *scale);
int balbak(int64 *nm, int64 *n, int64 *low, int64 *igh, double *scale, int64 *m,
           double *z__);
int elmhes(int64 *nm, int64 *n, int64 *low, int64 *igh, double *a,
           int64 *int__);
int eltran(int64 *nm, int64 *n, int64 *low, int64 *igh, double *a, int64 *int__,
           double *z__);
int qzhes(int64 nm, int64 n, double *a, double *b, logical matz, double *z__);
int qzit(int64 nm, int64 n, double *a, double *b, double eps1, logical matz,
         double *z__, int64 *ierr);
int qzval(int64 nm, int64 n, double *a, double *b, double *alfr, double *alfi,
          double *beta, logical matz, double *z__);
double epslon(double x);
double dnrm2(int64 *n, double *dx, int64 *incx);
double ddot(int64 *n, double *dx, int64 *incx, double *dy, int64 *incy);
int dscal(int64 *n, double *da, double *dx, int64 *incx);
int64 idamax(int64 *n, double *dx, int64 *incx);
int daxpy(int64 *n, double *da, double *dx, int64 *incx, double *dy,
          int64 *incy);
int drot(int64 *n, double *dx, int64 *incx, double *dy, int64 *incy, double *c,
         double *s);
int dswap(int64 *n, double *dx, int64 *incx, double *dy, int64 *incy);
int dgemc(int64 *m, int64 *n, double *a, int64 *lda, double *b, int64 *ldb,
          logical *trans);
int xerbla(char *srname, int64 *info, int64 srname_len);
logical lsame(char *ca, char *cb, int64 ca_len, int64 cb_len);
int dgemm(char *transa, char *transb, int64 *m, int64 *n, int64 *k,
          double *alpha, double *a, int64 *lda, double *b, int64 *ldb,
          double *beta, double *c, int64 *ldc, int64 transa_len,
          int64 transb_len);
int ezsvd(double *x, int64 *ldx, int64 *n, int64 *p, double *s, double *e,
          double *u, int64 *ldu, double *v, int64 *ldv, double *work,
          int64 *job, int64 *info, double *tol);
int ndrotg(double *f, double *g, double *cs, double *sn);
int ndsvd(double *x, int64 *ldx, int64 *n, int64 *p, double *s, double *e,
          double *u, int64 *ldu, double *v, int64 *ldv, double *work,
          int64 *job, int64 *info, int64 *maxitr, double *tol, int64 *idbg,
          int64 *ifull, int64 *kount, int64 *kount1, int64 *kount2, int64 *skip,
          int64 *limshf, double *maxsin, int64 *iidir);
int prse(int64 *ll, int64 *m, int64 *nrow, int64 *ncol, double *s, double *e);
int sig22(double *a, double *b, double *c, double *sigmin, double *sigmax,
          double *snr, double *csr, double *snl, double *csl);
double sigmin(double *a, double *b, double *c);
int sndrtg(double *f, double *g, double *cs, double *sn);
int hqr3lc(double *a, double *v, int64 *n, int64 *nlow, int64 *nup, double *eps,
           double *er, double *ei, int64 *type__, int64 *na, int64 *nv,
           int64 *imfd);
int split(double *a, double *v, int64 *n, int64 *l, double *e1, double *e2,
          int64 *na, int64 *nv);
int exchng(double *a, double *v, int64 *n, int64 *l, int64 *b1, int64 *b2,
           double *eps, logical *fail, int64 *na, int64 *nv);
int qrstep(double *a, double *v, double *p, double *q, double *r__, int64 *nl,
           int64 *nu, int64 *n, int64 *na, int64 *nv);
int orthes(int64 *nm, int64 *n, int64 *low, int64 *igh, double *a, double *ort);
int ortran(int64 *nm, int64 *n, int64 *low, int64 *igh, double *a, double *ort,
           double *z__);

/* problem defined functions*/
int func(int64 ndim, const double *u, const int64 *icp, const double *par,
         int64 ijac, double *f, double *dfdu, double *dfdp);
int stpnt(int64 ndim, double t, double *u, double *par);
int bcnd(int64 ndim, const double *par, const int64 *icp, int64 nbc,
         const double *u0, const double *u1, int64 ijac, double *f,
         double *dbc);
int icnd(int64 ndim, const double *par, const int64 *icp, int64 nint,
         const double *u, const double *uold, const double *udot,
         const double *upold, int64 ijac, double *fi, double *dint);
int fopt(int64 ndim, const double *u, const int64 *icp, const double *par,
         int64 ijac, double *fs, double *dfdu, double *dfdp);
int pvls(int64 ndim, const double *u, double *par);
/* conpar.c */
void *conpar_process(void *);
int conpar(int64 *nov, int64 *na, int64 *nra, int64 *nca, double *a, int64 *ncb,
           double *b, int64 *nbc, int64 *nrc, double *c, double *d, int64 *irf,
           int64 *icf);
/*setubv.c */
#include "auto_types.h"
void *setubv_make_aa_bb_cc(void *);
int setubv(int64 ndim, int64 ips, int64 na, int64 ncol, int64 nbc, int64 nint,
           int64 ncb, int64 nrc, int64 nra, int64 nca, FUNI_TYPE((*funi)),
           BCNI_TYPE((*bcni)), ICNI_TYPE((*icni)), int64 ndxloc, iap_type *iap,
           rap_type *rap, double *par, int64 *icp, double rds, double *aa,
           double *bb, double *cc, double *dd, double *fa, double *fc,
           double *rlcur, double *rlold, double *rldot, double *ups,
           double *uoldps, double *udotps, double *upoldp, double *dups,
           double *dtm, double *thl, double *thu, double *p0, double *p1);
void setubv_parallel_arglist_copy(setubv_parallel_arglist *output,
                                  const setubv_parallel_arglist input);
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
