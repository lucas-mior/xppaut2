/* This structure contains all of the input data for the setubv routine
   Those values which are arrays and those
   which are input and output are markered as such*/
#ifndef AUTO_TYPES_H
#define AUTO_TYPES_H
#include "integers.h"
typedef struct {
    int64 ndim, ips, ncol, nbc, nint, ncb, nrc, nra, nca, na; /*scalar input */
    FUNI_TYPE((*funi));                                       /*scalar input*/
    ICNI_TYPE((*icni));                                       /*scalar input*/
    int64 ndxloc;                                             /*scalar input*/
    iap_type *iap; /*array input size: NIAP*/
    rap_type *rap; /*array input size: NRAP*/
    double *par;   /*array input size: NPARX2*/
    int64 *icp;    /*array input size:  NPARX2*/
    double *aa;    /*array output (but must be initialized to 0) size: *nca X
                    *nra X *na */
    double *bb;    /*array output (but must be initialized to 0) size: *ncb X
                    *nra X *na */
    double *cc;    /*array output (but must be initialized to 0) size: *nca X
                    *nrc X *na */
    double
        *dd; /*array output (but must be initialized to 0) size: *ncb X *nrc */
    double
        *fa; /*array output (but must be initialized to 0) size: *nra X *na */
    double *fc;     /*array output (but must be initialized to 0) size: *nrc */
    double *ups;    /*array input size: *ndxloc X (*ndim X *ncol) */
    double *uoldps; /*array input size: *ndxloc X (*ndim X *ncol) */
    double *udotps; /*array input size: *ndxloc X (*ndim X *ncol) */
    double *upoldp; /*array input size: *ndxloc X (*ndim X *ncol) */
    double *dtm;    /*array input size: *na */
    int64 loop_start;  /*scalar input*/
    int64 loop_end;    /*scalar input*/
    int64 loop_offset; /*scalar input*/
    double *wp;        /*array input size: MCL2*MCL1 */
    double *wt;        /*array input size: MCL2*MCL1 */
    double *wi;    /*array input size: MCL2*MCL1??? Not sure of this one yet */
    double *thu;   /*array input size: ndim*8 */
    double *thl;   /*array input size: NPARX */
    double *rldot; /*array input size: NPARX */
    BCNI_TYPE((*bcni));
} setubv_parallel_arglist;

/* This structure contains all of the input data for the conpar2 routine
   Those values which are arrays and those
   which are input and output are markered as such*/
typedef struct {
    int64 *nov, *nra, *nca; /*scalars input*/
    double *a;              /*array input and output size: nca X nra X na */
    int64 *ncb;             /*scalar input */
    double *b;              /*array input and output size: ncb X nra X na*/
    int64 *nbc, *nrc;       /*scalar input */
    double *c;              /*array input and output size: nca X nrc X *na*/
    double *d;              /*array input and output size: ncb X nrc*/
    int64 *irf;             /*array input size: na X nra*/
    int64 *icf;             /*array input: na X nca*/
    int64 loop_start;       /*scalar input*/
    int64 loop_end;         /*scalar output*/
} conpar_parallel_arglist;
#endif
