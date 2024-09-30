/******************************************************************
 *                                                                *
 * File          : band.h                                         *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the header file for a generic BAND linear solver       * package.
 * There are two sets of band solver routines listed in  * this file: one set
 * uses type BandMat defined below and the     * other set uses the type double
 * ** for band matrix arguments.     * The two sets of band solver routines make
 * it easy to work      * with two types of band matrices:
 * *
 *                                                                *
 * (1) The BandMat type is intended for use with large            *
 *     band matrices whose elements/columns may be stored in      *
 *     non-contiguous memory locations or even distributed into   *
 *     different processor memories. This type may be modified to *
 *     include such distribution information. If this is done,    *
 *     then all the routines that use BandMat must be modified to *
 *     reflect the new data structure.                            *
 *                                                                *
 * (2) The set of routines that use double ** (and NOT the BandMat  *
 *     type) is intended for use with small matrices which can    *
 *     easily be allocated within a contiguous block of memory    *
 *     on a single processor.                                     *
 *                                                                *
 * Routines that work with the type BandMat begin with "Band".    *
 * The band_alloc_mat function allocates a band matrix for use in   *
 * the other matrix routines listed in this file. Matrix storage  *
 * details are given in the documentation for the type BandMat.   *
 * The BandAllocPiv function allocates memory for pivot           *
 * information. The storage allocated by band_alloc_mat and         *
 * BandAllocPiv is deallocated by the routines band_free_mat and    *
 * band_free_piv, respectively. The band_factor and band_back_solve    *
 * routines perform the actual solution of a band linear system.  *
 * Note that the band_back_solve routine has a parameter b of type  *
 * Vector. The current implementation makes use of a machine    *
 * environment specific macro (N_VDATA) which may not exist for   *
 * other implementations of the type Vector. Thus, the          *
 * implementation of band_back_solve may need to change if the      *
 * type Vector is changed.                                      *
 *                                                                *
 * Routines that work with double ** begin with "band" (except for  *
 * the factor and solve routines which are called gbfa and gbsl,  *
 * respectively). The underlying matrix storage is described in   *
 * the documentation for bandalloc.                               *
 *                                                                *
 ******************************************************************/

#ifndef BAND_H
#define BAND_H

#include "vector.h"

/******************************************************************
 * Type: BandMat                                                  *
 *----------------------------------------------------------------*
 * The type BandMat is the type of a large (possibly distributed) *
 * band matrix. It is defined to be a pointer to a structure      *
 * with the following fields:                                     *
 *                                                                *
 * size is the number of columns (== number of rows)              *
 * mu   is the upper bandwidth, 0 <= mu <= size-1                 *
 * ml   is the lower bandwidth, 0 <= ml <= size-1                 *
 * smu  is the storage upper bandwidth, mu <= smu <= size-1.      *
 *         The band_factor routine writes the LU factors           *
 *         into the storage for A. The upper triangular factor U, *
 *         however, may have an upper bandwidth as big as         *
 *         MIN(size-1,mu+ml) because of partial pivoting. The smu *
 *         field holds the upper bandwidth allocated for A.       *
 * data is a two dimensional array used for component storage.    *
 *         The elements of a band matrix of type BandMat are      *
 *         stored columnwise (i.e. columns are stored one on top  *
 *         of the other in memory). Only elements within the      *
 *         specified bandwidths are stored.                       *
 *                                                                *
 * If we number rows and columns in the band matrix starting      *
 * from 0, then                                                   *
 *                                                                *
 * data[0] is a pointer to (smu+ml+1)*size contiguous locations   *
 *            which hold the elements within the band of A        *
 * data[j] is a pointer to the uppermost element within the band  *
 *            in the jth column. This pointer may be treated as   *
 *            an array indexed from smu-mu (to access the         *
 *            uppermost element within the band in the jth        *
 *            column) to smu+ml (to access the lowest element     *
 *            within the band in the jth column). (Indices from 0 *
 *            to smu-mu-1 give access to extra storage elements   *
 *            required by band_factor.)                            *
 * data[j][i-j+smu] is the (i,j)th element, j-mu <= i <= j+ml.    *
 *                                                                *
 * The macros below allow a user to access individual matrix      *
 * elements without writing out explicit data structure           *
 * references and without knowing too much about the underlying   *
 * element storage. The only storage assumption needed is that    *
 * elements are stored columnwise and that a pointer into the jth *
 * column of elements can be obtained via the BAND_COL macro. The *
 * BAND_COL_ELEM macro selects an element from a column which has *
 * already been isolated via BAND_COL. BAND_COL_ELEM allows the   *
 * user to avoid the translation from the matrix location (i,j)   *
 * to the index in the array returned by BAND_COL at which the    *
 * (i,j)th element is stored. See the documentation for BAND_COL  *
 * and BAND_COL_ELEM for usage details. Users should use these    *
 * macros whenever possible.                                      *
 *                                                                *
 ******************************************************************/
typedef struct {
    int64 size;
    int64 mu;
    int64 ml;
    int64 smu;
    double **data;
} *BandMat;

/* BandMat accessor macros */

/******************************************************************
 * Macro : BAND_ELEM                                              *
 * Usage : BAND_ELEM(A,i,j) = a_ij;  OR                           *
 *         a_ij = BAND_ELEM(A,i,j);                               *
 *----------------------------------------------------------------*
 * BAND_ELEM(A,i,j) references the (i,j)th element of the         *
 * N by N band matrix A, where 0 <= i,j <= N-1. The location      *
 * (i,j) should further satisfy j-(A->mu) <= i <= j+(A->ml).      *
 ******************************************************************/
#define BAND_ELEM(A, i, j) ((A->data)[j][i - j + (A->smu)])

/******************************************************************
 * Macro : BAND_COL                                               *
 * Usage : col_j = BAND_COL(A,j);                                 *
 *----------------------------------------------------------------*
 * BAND_COL(A,j) references the diagonal element of the jth       *
 * column of the N by N band matrix A, 0 <= j <= N-1. The type of *
 * the expression BAND_COL(A,j) is double *. The pointer returned   *
 * by the call BAND_COL(A,j) can be treated as an array which is  *
 * indexed from -(A->mu) to (A->ml).                              *
 ******************************************************************/
#define BAND_COL(A, j) (((A->data)[j]) + (A->smu))

/******************************************************************
 * Macro : BAND_COL_ELEM                                          *
 * Usage : col_j = BAND_COL(A,j);                                 *
 *         BAND_COL_ELEM(col_j,i,j) = a_ij;  OR                   *
 *         a_ij = BAND_COL_ELEM(col_j,i,j);                       *
 *----------------------------------------------------------------*
 * This macro references the (i,j)th entry of the band matrix A   *
 * when used in conjunction with BAND_COL as shown above. The     *
 * index (i,j) should satisfy j-(A->mu) <= i <= j+(A->ml).        *
 ******************************************************************/
#define BAND_COL_ELEM(col_j, i, j) (col_j[i - j])

/* Functions that use the BandMat representation for a band matrix */

/******************************************************************
 * Function : band_alloc_mat                                        *
 * Usage    : A = band_alloc_mat(N, mu, ml, smu);                   *
 *            if (A == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * band_alloc_mat allocates memory for an N by N band matrix with   *
 * upper bandwidth mu, lower bandwidth ml, and storage upper      *
 * bandwidth smu. Pass smu as follows depending on whether A will *
 * be factored by band_factor:                                     *
 *                                                                *
 * (1) Pass smu = mu if A will not be factored.                   *
 * (2) Pass smu = MIN(N-1,mu+ml) if A will be factored.           *
 *                                                                *
 * band_alloc_mat returns the storage allocated (type BandMat) or   *
 * NULL if the request for matrix storage cannot be satisfied.    *
 * See the documentation for the type BandMat for matrix storage  *
 * details.                                                       *
 ******************************************************************/
BandMat band_alloc_mat(int64 N, int64 mu, int64 ml, int64 smu);

/******************************************************************
 * Function : BandAllocPiv                                        *
 * Usage    : p = band_alloc_piv(N);                                *
 *            if (p == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * BandAllocPiv allocates memory for pivot information to be      *
 * filled in by the band_factor routine during the factorization   *
 * of an N by N band matrix. The underlying type for pivot        *
 * information is an array of N integers and this routine returns *
 * the pointer to the memory it allocates. If the request for     *
 * pivot storage cannot be satisfied, BandAllocPiv returns NULL.  *
 ******************************************************************/
int64 *band_alloc_piv(int64 N);

/******************************************************************
 * Function : band_factor                                          *
 * Usage    : ier = band_factor(A, p);                             *
 *            if (ier != 0) ... A is singular                     *
 *----------------------------------------------------------------*
 * band_factor performs the LU factorization of the N by N band    *
 * matrix A. This is done using standard Gaussian elimination     *
 * with partial pivoting.                                         *
 *                                                                *
 * A successful LU factorization leaves the "matrix" A and the    *
 * pivot array p with the following information:                  *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., N-1.  *
 * (2) If the unique LU factorization of A is given by PA = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of A     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of A contains the multipliers, I-L.        *
 *                                                                *
 * band_factor returns 0 if successful. Otherwise it encountered   *
 * a zero diagonal element during the factorization. In this case *
 * it returns the column index (numbered from one) at which       *
 * it encountered the zero.                                       *
 *                                                                *
 * Important Note: A must allocated to accomodate the increase    *
 * in upper bandwidth that occurs during factorization. If,       *
 * mathematically, A is a band matrix with upper bandwidth mu and *
 * lower bandwidth ml, then the upper triangular factor U can     *
 * have upper bandwidth as big as smu=MIN(n-1,mu+ml). The lower   *
 * triangular factor L has lower bandwidth ml. Allocate A with    *
 * call A = band_alloc_mat(N,mu,ml,smu), where mu, ml, and smu are  *
 * as defined above. The user does not have to zero the "extra"   *
 * storage allocated for the purpose of factorization. This will  *
 * handled by the band_factor routine.                             *
 ******************************************************************/
int64 band_factor(BandMat A, int64 *p);

/******************************************************************
 * Function : band_back_solve                                       *
 * Usage    : band_back_solve(A, p, b);                             *
 *----------------------------------------------------------------*
 * band_back_solve solves the N-dimensional system A x = b using    *
 * the LU factorization in A and the pivot information in p       *
 * computed in band_factor. The solution x is returned in b. This  *
 * routine cannot fail if the corresponding call to band_factor    *
 * did not fail.                                                  *
 ******************************************************************/
void band_back_solve(BandMat A, int64 *p, Vector b);

/******************************************************************
 * Function : band_zero                                            *
 * Usage    : band_zero(A);                                        *
 *----------------------------------------------------------------*
 * A(i,j) <- 0.0,    j-(A->mu) <= i <= j+(A->ml).                 *
 ******************************************************************/
void band_zero(BandMat A);

/******************************************************************
 * Function : band_copy                                            *
 * Usage    : band_copy(A, B, copymu, copyml);                     *
 *----------------------------------------------------------------*
 * band_copy copies the submatrix with upper and lower bandwidths  *
 * copymu, copyml of the N by N band matrix A into the N by N     *
 * band matrix B.                                                 *
 ******************************************************************/
void band_copy(BandMat A, BandMat B, int64 copymu, int64 copyml);

/******************************************************************
 * Function: band_scale                                            *
 * Usage   : band_scale(c, A);                                     *
 *----------------------------------------------------------------*
 * A(i,j) <- c*A(i,j),   j-(A->mu) <= i <= j+(A->ml).             *
 ******************************************************************/
void band_scale(double c, BandMat A);

/******************************************************************
 * Function : band_add_i                                            *
 * Usage    : band_add_i(A);                                        *
 *----------------------------------------------------------------*
 * A(j,j) <- A(j,j)+1.0,   0 <= j <= (A->size)-1.                 *
 ******************************************************************/
void band_add_i(BandMat A);

/******************************************************************
 * Function : band_free_mat                                         *
 * Usage    : band_free_mat(A);                                     *
 *----------------------------------------------------------------*
 * band_free_mat frees the memory allocated by band_alloc_mat for     *
 * the band matrix A.                                             *
 ******************************************************************/
void band_free_mat(BandMat A);

/******************************************************************
 * Function : band_free_piv                                         *
 * Usage    : band_free_piv(p);                                     *
 *----------------------------------------------------------------*
 * band_free_piv frees the memory allocated by BandAllocPiv for     *
 * the pivot information array p.                                 *
 ******************************************************************/
void band_free_piv(int64 *p);

/******************************************************************
 * Function : band_print                                           *
 * Usage    : band_print(A);                                       *
 *----------------------------------------------------------------*
 * This routine prints the N by N band matrix A (upper and lower  *
 * bandwidths A->mu and A->ml, respectively) to standard output   *
 * as it would normally appear on paper. It is intended as a      *
 * debugging tool with small values of N. The elements are        *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 ******************************************************************/
void band_print(BandMat A);

/* Functions that use the double ** representation for a band matrix */

/******************************************************************
 * Function : bandalloc                                           *
 * Usage    : double **a;                                           *
 *            a = band_alloc2(n, smu, ml);                          *
 *            if (a == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * band_alloc2(n, smu, ml) allocates storage for an n by n band     *
 * matrix A with storage upper bandwidth smu and lower bandwidth  *
 * ml. It returns a pointer to the newly allocated storage if     *
 * successful. If the memory request cannot be satisfied, then    *
 * bandalloc returns NULL. If, mathematically, A has upper and    *
 * lower bandwidths mu and ml, respectively, then the value       *
 * passed to bandalloc for smu may need to be greater than mu.    *
 * The gbfa routine writes the LU factors into the storage (named *
 * "a" in the above usage documentation) for A (thus destroying   *
 * the original elements of A). The upper triangular factor U,    *
 * however, may have a larger upper bandwidth than the upper      *
 * bandwidth mu of A. Thus some "extra" storage for A must be     *
 * allocated if A is to be factored by gbfa. Pass smu as follows: *
 *                                                                *
 * (1) Pass smu = mu if A will not be factored.                   *
 * (2) Pass smu = MIN(n-1,mu+ml) if A will be factored.           *
 *                                                                *
 * The underlying type of the band matrix returned is double **. If *
 * we allocate a band matrix A in double **a by                     *
 * a = band_alloc2(n,smu,ml), then a[0] is a pointer to             *
 * n*(smu + ml + 1) contiguous storage locations and a[j] is a  *
 * pointer to the uppermost element in the storage for the jth    *
 * column. The expression a[j][i-j+smu] references the (i,j)th    *
 * element of A, where 0 <= i,j <= n-1 and j-mu <= i <= j+ml.     *
 * (The elements a[j][0], a[j][1], ..., a[j][smu-mu-1] are used   *
 * by gbfa and gbsl.)                                             *
 ******************************************************************/
double **band_alloc2(int64 n, int64 smu, int64 ml);

/******************************************************************
 * Function : bandallocpiv                                        *
 * Usage    : int64 *pivot;                                     *
 *            pivot = band_alloc_piv2(n);                            *
 *            if (pivot == NULL) ... memory request failed        *
 *----------------------------------------------------------------*
 * band_alloc_piv2(n) allocates an array of n integers. It returns a *
 * pointer to the first element in the array if successful. It    *
 * returns NULL if the memory request could not be satisfied.     *
 ******************************************************************/
int64 *band_alloc_piv2(int64 n);

/******************************************************************
 * Function : gbfa                                                *
 * Usage    : int64 ier;                                        *
 *            ier = band_gbfa(a,n,mu,ml,smu,p);                        *
 *            if (ier > 0) ... zero element encountered during    *
 *                             the factorization                  *
 *----------------------------------------------------------------*
 * band_gbfa(a,n,mu,ml,smu,p) factors the n by n band matrix A (upper  *
 * and lower bandwidths mu and ml, storage upper bandwidth smu)   *
 * stored in "a". It overwrites the elements of A with the LU     *
 * factors and it keeps track of the pivot rows chosen in the     *
 * pivot array p.                                                 *
 *                                                                *
 * A successful LU factorization leaves a and pivot array p with  *
 * the following information:                                     *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., n-1.  *
 * (2) If the unique LU factorization of A is given by PA = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of A     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of A contains the multipliers, I-L.        *
 *                                                                *
 * gbfa returns 0 if successful. Otherwise it encountered a zero  *
 * diagonal element during the factorization. In this case it     *
 * returns the column index (numbered from one) at which it       *
 * encountered the zero.                                          *
 *                                                                *
 * IMPORTANT NOTE: Suppose A is a band matrix with upper          *
 * bandwidth mu and lower bandwidth ml, then the upper triangular *
 * factor U can have upper bandwidth as big as MIN(n-1,mu+ml)     *
 * because of partial pivoting. The lower triangular factor L has *
 * lower bandwidth ml. Thus, if A is to be factored and           *
 * backsolved using gbfa and gbsl, then it should be allocated    *
 * as a = band_alloc2(n,smu,ml), where smu = MIN(n-1,mu+ml). The    *
 * call to gbfa is ier = band_gbfa(a,n,mu,ml,smu,p). The corresponding *
 * call to gbsl is band_gbsl(a,n,smu,ml,p,b). The user does not need   *
 * to zero the "extra" storage allocated for the purpose of       *
 * factorization. This is handled by the gbfa routine. If A is    *
 * not going to be factored and backsolved, then it can be        *
 * allocated as a = band_alloc2(n,smu,ml). In either case, all      *
 * routines in this section use the parameter name smu for a      *
 * parameter which must be the "storage upper bandwidth" which    *
 * was passed to bandalloc.                                       *
 ******************************************************************/
int64 band_gbfa(double **a, int64 n, int64 mu, int64 ml, int64 smu, int64 *p);

/******************************************************************
 * Function : gbsl                                                *
 * Usage    : double *b;                                            *
 *            ier = band_gbfa(a,n,mu,ml,smu,p);                        *
 *            if (ier == 0) band_gbsl(a,n,smu,ml,p,b);                 *
 *----------------------------------------------------------------*
 * band_gbsl(a,n,smu,ml,p,b) solves the n by n linear system           *
 * Ax = b, where A is band matrix stored in "a" with storage      *
 * upper bandwidth smu and lower bandwidth ml. It assumes that A  *
 * has been LU factored and the pivot array p has been set by a   *
 * successful call band_gbfa(a,n,mu,ml,smu,p). The solution x is       *
 * written into the b array.                                      *
 ******************************************************************/
void band_gbsl(double **a, int64 n, int64 smu, int64 ml, int64 *p, double *b);

/******************************************************************
 * Function : band_zero2                                            *
 * Usage    : band_zero2(a,n,mu,ml,smu);                            *
 *----------------------------------------------------------------*
 * a(i,j) <- 0.0,   0 <= i,j <= n-1, j-mu <= i <= j+ml.           *
 ******************************************************************/
void band_zero2(double **a, int64 n, int64 mu, int64 ml, int64 smu);

/******************************************************************
 * Function : band_copy2                                            *
 * Usage    : band_copy2(a,b,n,a_smu,b_smu,copymu,copyml);          *
 *----------------------------------------------------------------*
 * b(i,j) <- a(i,j), 0 <= i,j <= n-1, j-copymu <= i <= j+copyml.  *
 ******************************************************************/
void band_copy2(double **a, double **b, int64 n, int64 a_smu, int64 b_smu, int64 copymu,
                int64 copyml);

/******************************************************************
 * Function : band_scale2                                           *
 * Usage    : band_scale2(c,a,n,mu,ml);                             *
 *----------------------------------------------------------------*
 * a(i,j) <- c*a(i,j),   0 <= i,j <= n-1, j-mu <= i <= j+ml.      *
 ******************************************************************/
void band_scale2(double c, double **a, int64 n, int64 mu, int64 ml, int64 smu);

/******************************************************************
 * Function : band_add_i2                                            *
 * Usage    : band_add_i2(a,n,smu);                                  *
 *----------------------------------------------------------------*
 * a(j,j) <- a(j,j)+1.0,   0 <= j <= n-1.                         *
 ******************************************************************/
void band_add_i2(double **a, int64 n, int64 smu);

/******************************************************************
 * Function : band_free_pvi2                                         *
 * Usage    : band_free_pvi2(p);                                     *
 *----------------------------------------------------------------*
 * band_free_pvi2(p) frees the pivot array p allocated by            *
 * bandallocpiv.                                                  *
 ******************************************************************/
void band_free_pvi2(int64 *p);

/******************************************************************
 * Function : band_free2                                            *
 * Usage    : band_free2(a);                                        *
 *----------------------------------------------------------------*
 * band_free2(a) frees the band matrix a allocated by bandalloc.    *
 ******************************************************************/
void band_free2(double **a);

/******************************************************************
 * Function : band_print2                                           *
 * Usage    : band_print2(a,n,mu,ml,smu);                           *
 *----------------------------------------------------------------*
 * band_print2(a,n,mu,ml,smu) prints the n by n band matrix stored  *
 * in a (with upper bandwidth mu and lower bandwidth ml) to       *
 * standard output as it would normally appear on paper. It is    *
 * intended as a debugging tool with small values of n. The       *
 * elements are printed using the %g option. A blank line is      *
 * printed before and after the matrix.                           *
 ******************************************************************/
void band_print2(double **a, int64 n, int64 mu, int64 ml, int64 smu);

#endif
