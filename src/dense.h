/******************************************************************
 *                                                                *
 * File          : dense.h                                        *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the header file for a generic DENSE linear solver      *
 * package. There are two sets of dense solver routines listed in *
 * this file: one set uses type DenseMat defined below and the    *
 * other set uses the type double ** for dense matrix arguments.    *
 * The two sets of dense solver routines make it easy to work     *
 * with two types of dense matrices:                              *
 *                                                                *
 * (1) The DenseMat type is intended for use with large dense     *
 *     matrices whose elements/columns may be stored in           *
 *     non-contiguous memory locations or even distributed into   *
 *     different processor memories. This type may be modified to *
 *     include such distribution information. If this is done,    *
 *     then all the routines that use DenseMat must be modified   *
 *     to reflect the new data structure.                         *
 *                                                                *
 * (2) The set of routines that use double ** (and NOT the DenseMat *
 *     type) is intended for use with small matrices which can    *
 *     easily be allocated within a contiguous block of memory    *
 *     on a single processor.                                     *
 *                                                                *
 * Routines that work with the type DenseMat begin with "Dense".  *
 * The dense_alloc_mat function allocates a dense matrix for use in *
 * the other DenseMat routines listed in this file. Matrix        *
 * storage details are given in the documentation for the type    *
 * DenseMat. The DenseAllocPiv function allocates memory for      *
 * pivot information. The storage allocated by dense_alloc_mat and  *
 * DenseAllocPiv is deallocated by the routines dense_free_mat and  *
 * dense_free_piv, respectively. The dense_factor and dense_back_solve *
 * routines perform the actual solution of a dense linear system. *
 * Note that the dense_back_solve routine has a parameter b of type *
 * Vector. The current implementation makes use of a machine    *
 * environment-specific macro (N_VDATA) which may not exist for   *
 * other implementations of the type Vector. Thus, the          *
 * implementation of dense_back_solve may need to change if the     *
 * type Vector is changed.                                      *
 *                                                                *
 * Routines that work with double ** begin with "den" (except for   *
 * the factor and solve routines which are called gefa and gesl,  *
 * respectively). The underlying matrix storage is described in   *
 * the documentation for denalloc.                                *
 *                                                                *
 ******************************************************************/

#ifndef DENSE_H
#define DENSE_H

#include "vector.h"

/******************************************************************
 *                                                                *
 * Type: DenseMat                                                 *
 *----------------------------------------------------------------*
 * The type DenseMat is defined to be a pointer to a structure    *
 * with a size and a data field. The size field indicates the     *
 * number of columns (== number of rows) of a dense matrix, while *
 * the data field is a two dimensional array used for component   *
 * storage. The elements of a dense matrix are stored columnwise  *
 * (i.e columns are stored one on top of the other in memory). If *
 * A is of type DenseMat, then the (i,j)th element of A (with     *
 * 0 <= i,j <= size-1) is given by the expression (A->data)[j][i] *
 * or by the expression (A->data)[0][j*n+i]. The macros below     *
 * allow a user to access efficiently individual matrix           *
 * elements without writing out explicit data structure           *
 * references and without knowing too much about the underlying   *
 * element storage. The only storage assumption needed is that    *
 * elements are stored columnwise and that a pointer to the jth   *
 * column of elements can be obtained via the DENSE_COL macro.    *
 * Users should use these macros whenever possible.               *
 *                                                                *
 ******************************************************************/
typedef struct {
    int64 size;
    double **data;
} *DenseMat;

/* DenseMat accessor macros */

/******************************************************************
 *                                                                *
 * Macro : DENSE_ELEM                                             *
 * Usage : DENSE_ELEM(A,i,j) = a_ij;  OR                          *
 *         a_ij = DENSE_ELEM(A,i,j);                              *
 *----------------------------------------------------------------*
 * DENSE_ELEM(A,i,j) references the (i,j)th element of the N by N *
 * DenseMat A, 0 <= i,j <= N-1.                                   *
 *                                                                *
 ******************************************************************/
#define DENSE_ELEM(A, i, j) ((A->data)[j][i])

/******************************************************************
 *                                                                *
 * Macro : DENSE_COL                                              *
 * Usage : col_j = DENSE_COL(A,j);                                *
 *----------------------------------------------------------------*
 * DENSE_COL(A,j) references the jth column of the N by N         *
 * DenseMat A, 0 <= j <= N-1. The type of the expression          *
 * DENSE_COL(A,j) is double *. After the assignment in the usage    *
 * above, col_j may be treated as an array indexed from 0 to N-1. *
 * The (i,j)th element of A is referenced by col_j[i].            *
 *                                                                *
 ******************************************************************/
#define DENSE_COL(A, j) ((A->data)[j])

/* Functions that use the DenseMat representation for a dense matrix */

/******************************************************************
 *                                                                *
 * Function : dense_alloc_mat                                       *
 * Usage    : A = dense_alloc_mat(N);                               *
 *            if (A == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * dense_alloc_mat allocates memory for an N by N dense matrix and  *
 * returns the storage allocated (type DenseMat). dense_alloc_mat   *
 * returns NULL if the request for matrix storage cannot be       *
 * satisfied. See the above documentation for the type DenseMat   *
 * for matrix storage details.                                    *
 *                                                                *
 ******************************************************************/
DenseMat dense_alloc_mat(int64 N);

/******************************************************************
 *                                                                *
 * Function : DenseAllocPiv                                       *
 * Usage    : p = dense_alloc_piv(N);                               *
 *            if (p == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * DenseAllocPiv allocates memory for pivot information to be     *
 * filled in by the dense_factor routine during the factorization  *
 * of an N by N dense matrix. The underlying type for pivot       *
 * information is an array of N integers and this routine returns *
 * the pointer to the memory it allocates. If the request for     *
 * pivot storage cannot be satisfied, DenseAllocPiv returns NULL. *
 *                                                                *
 ******************************************************************/
int64 *dense_alloc_piv(int64 N);

/******************************************************************
 *                                                                *
 * Function : dense_factor                                         *
 * Usage    : ier = dense_factor(A, p);                            *
 *            if (ier != 0) ... A is singular                     *
 *----------------------------------------------------------------*
 * dense_factor performs the LU factorization of the N by N dense  *
 * matrix A. This is done using standard Gaussian elimination     *
 * with partial pivoting.                                         *
 *                                                                *
 * A successful LU factorization leaves the matrix A and the      *
 * pivot array p with the following information:                  *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., N-1.  *
 *                                                                *
 * (2) If the unique LU factorization of A is given by PA = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of A     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of A contains the multipliers, I-L.        *
 *                                                                *
 * dense_factor returns 0 if successful. Otherwise it encountered  *
 * a zero diagonal element during the factorization. In this case *
 * it returns the column index (numbered from one) at which       *
 * it encountered the zero.                                       *
 *                                                                *
 ******************************************************************/
int64 dense_factor(DenseMat A, int64 *p);

/******************************************************************
 *                                                                *
 * Function : dense_back_solve                                      *
 * Usage    : dense_back_solve(A, p, b);                            *
 *----------------------------------------------------------------*
 * dense_back_solve solves the N-dimensional system A x = b using   *
 * the LU factorization in A and the pivot information in p       *
 * computed in dense_factor. The solution x is returned in b. This *
 * routine cannot fail if the corresponding call to dense_factor   *
 * did not fail.                                                  *
 *                                                                *
 ******************************************************************/
void dense_back_solve(DenseMat A, int64 *p, Vector b);

/******************************************************************
 *                                                                *
 * Function : dense_zero                                           *
 * Usage    : dense_zero(A);                                       *
 *----------------------------------------------------------------*
 * dense_zero sets all the elements of the N by N matrix A to 0.0. *
 *                                                                *
 ******************************************************************/
void dense_zero(DenseMat A);

/******************************************************************
 *                                                                *
 * Function : dense_copy                                           *
 * Usage    : dense_copy(A, B);                                    *
 *----------------------------------------------------------------*
 * dense_copy copies the contents of the N by N matrix A into the  *
 * N by N matrix B.                                               *
 *                                                                *
 ******************************************************************/
void dense_copy(DenseMat A, DenseMat B);

/******************************************************************
 *                                                                *
 * Function: dense_scal                                           *
 * Usage   : dense_scal(c, A);                                    *
 *----------------------------------------------------------------*
 * dense_scal scales the elements of the N by N matrix A by the   *
 * constant c and stores the result back in A.                    *
 *                                                                *
 ******************************************************************/
void dense_scal(double c, DenseMat A);

/******************************************************************
 *                                                                *
 * Function : dense_add_i                                           *
 * Usage    : dense_add_i(A);                                       *
 *----------------------------------------------------------------*
 * dense_add_i adds the identity matrix to A and stores the result  *
 * back in A.                                                     *
 *                                                                *
 ******************************************************************/
void dense_add_i(DenseMat A);

/******************************************************************
 *                                                                *
 * Function : dense_free_mat                                        *
 * Usage    : dense_free_mat(A);                                    *
 *----------------------------------------------------------------*
 * dense_free_mat frees the memory allocated by dense_alloc_mat for   *
 * the N by N matrix A.                                           *
 *                                                                *
 ******************************************************************/
void dense_free_mat(DenseMat A);

/******************************************************************
 *                                                                *
 * Function : dense_free_piv                                        *
 * Usage    : dense_free_piv(p);                                    *
 *----------------------------------------------------------------*
 * dense_free_piv frees the memory allocated by DenseAllocPiv for   *
 * the pivot information array p.                                 *
 *                                                                *
 ******************************************************************/
void dense_free_piv(int64 *p);

/******************************************************************
 *                                                                *
 * Function : dense_print                                          *
 * Usage    : dense_print(A);                                      *
 *----------------------------------------------------------------*
 * This routine prints the N by N dense matrix A to standard      *
 * output as it would normally appear on paper. It is intended    *
 * as a debugging tool with small values of N. The elements are   *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 *                                                                *
 ******************************************************************/
void dense_print(DenseMat A);

/* Functions that use the double ** representation for a dense matrix */

/******************************************************************
 *                                                                *
 * Function : denalloc                                            *
 * Usage    : double **a;                                           *
 *            a = dense_alloc2(n);                                    *
 *            if (a == NULL) ... memory request failed            *
 *----------------------------------------------------------------*
 * dense_alloc2(n) allocates storage for an n by n dense matrix. It   *
 * returns a pointer to the newly allocated storage if            *
 * successful. If the memory request cannot be satisfied, then    *
 * denalloc returns NULL. The underlying type of the dense matrix *
 * returned is double **. If we allocate a dense matrix double **a by *
 * a = dense_alloc2(n), then a[j][i] references the (i,j)th element   *
 * of the matrix a, 0 <= i,j <= n-1, and a[j] is a pointer to the *
 * first element in the jth column of a. The location a[0]        *
 * contains a pointer to n^2 contiguous locations which contain   *
 * the elements of a.                                             *
 *                                                                *
 ******************************************************************/
double **dense_alloc2(int64 n);

/******************************************************************
 *                                                                *
 * Function : denallocpiv                                         *
 * Usage    : int64 *pivot;                                     *
 *            pivot = dense_alloc_piv2(n);                             *
 *            if (pivot == NULL) ... memory request failed        *
 *----------------------------------------------------------------*
 * dense_alloc_piv2(n) allocates an array of n integers. It returns a  *
 * pointer to the first element in the array if successful. It    *
 * returns NULL if the memory request could not be satisfied.     *
 *                                                                *
 ******************************************************************/
int64 *dense_alloc_piv2(int64 n);

/******************************************************************
 *                                                                *
 * Function : gefa                                                *
 * Usage    : int64 ier;                                        *
 *            ier = dense_gefa(a,n,p);                                  *
 *            if (ier > 0) ... zero element encountered during    *
 *                             the factorization                  *
 *----------------------------------------------------------------*
 * dense_gefa(a,n,p) factors the n by n dense matrix a. It overwrites   *
 * the elements of a with its LU factors and keeps track of the   *
 * pivot rows chosen in the pivot array p.                        *
 *                                                                *
 * A successful LU factorization leaves the matrix a and the      *
 * pivot array p with the following information:                  *
 *                                                                *
 * (1) p[k] contains the row number of the pivot element chosen   *
 *     at the beginning of elimination step k, k=0, 1, ..., n-1.  *
 *                                                                *
 * (2) If the unique LU factorization of a is given by Pa = LU,   *
 *     where P is a permutation matrix, L is a lower triangular   *
 *     matrix with all 1's on the diagonal, and U is an upper     *
 *     triangular matrix, then the upper triangular part of a     *
 *     (including its diagonal) contains U and the strictly lower *
 *     triangular part of a contains the multipliers, I-L.        *
 *                                                                *
 * gefa returns 0 if successful. Otherwise it encountered a zero  *
 * diagonal element during the factorization. In this case it     *
 * returns the column index (numbered from one) at which it       *
 * encountered the zero.                                          *
 *                                                                *
 ******************************************************************/
int64 dense_gefa(double **a, int64 n, int64 *p);

/******************************************************************
 *                                                                *
 * Function : gesl                                                *
 * Usage    : double *b;                                            *
 *            ier = dense_gefa(a,n,p);                                  *
 *            if (ier == 0) dense_gesl(a,n,p,b);                        *
 *----------------------------------------------------------------*
 * dense_gesl(a,n,p,b) solves the n by n linear system ax = b. It       *
 * assumes that a has been LU factored and the pivot array p has  *
 * been set by a successful call to dense_gefa(a,n,p). The solution x   *
 * is written into the b array.                                   *
 *                                                                *
 ******************************************************************/
void dense_gesl(double **a, int64 n, int64 *p, double *b);

/******************************************************************
 *                                                                *
 * Function : den_zero                                             *
 * Usage    : dense_zero2(a,n);                                       *
 *----------------------------------------------------------------*
 * dense_zero2(a,n) sets all the elements of the n by n dense matrix  *
 * a to be 0.0.                                                   *
 *                                                                *
 ******************************************************************/
void dense_zero2(double **a, int64 n);

/******************************************************************
 *                                                                *
 * Function : den_copy                                             *
 * Usage    : dense_copy2(a,b,n);                                     *
 *----------------------------------------------------------------*
 * dense_copy2(a,b,n) copies the n by n dense matrix a into the       *
 * n by n dense matrix b.                                         *
 *                                                                *
 ******************************************************************/
void dense_copy2(double **a, double **b, int64 n);

/******************************************************************
 *                                                                *
 * Function : den_scale                                            *
 * Usage    : dense_scale2(c,a,n);                                    *
 *----------------------------------------------------------------*
 * dense_scale2(c,a,n) scales every element in the n by n dense       *
 * matrix a by c.                                                 *
 *                                                                *
 ******************************************************************/
void dense_scale2(double c, double **a, int64 n);

/******************************************************************
 *                                                                *
 * Function : den_add_i                                             *
 * Usage    : dense_add_i2(a,n);                                       *
 *----------------------------------------------------------------*
 * dense_add_i2(a,n) increments the n by n dense matrix a by the       *
 * identity matrix.                                               *
 *                                                                *
 ******************************************************************/
void dense_add_i2(double **a, int64 n);

/******************************************************************
 *                                                                *
 * Function : den_free_piv                                          *
 * Usage    : dense_free_piv2(p);                                      *
 *----------------------------------------------------------------*
 * dense_free_piv2(p) frees the pivot array p allocated by             *
 * denallocpiv.                                                   *
 *                                                                *
 ******************************************************************/
void dense_free_piv2(int64 *p);

/******************************************************************
 *                                                                *
 * Function : den_free                                             *
 * Usage    : dense_free2(a);                                         *
 *----------------------------------------------------------------*
 * dense_free2(a) frees the dense matrix a allocated by denalloc.     *
 *                                                                *
 ******************************************************************/
void dense_free2(double **a);

/******************************************************************
 *                                                                *
 * Function : den_print                                            *
 * Usage    : dense_print2(a,n);                                      *
 *----------------------------------------------------------------*
 * dense_print2(a,n) prints the n by n dense matrix a to standard     *
 * output as it would normally appear on paper. It is intended as *
 * a debugging tool with small values of n. The elements are      *
 * printed using the %g option. A blank line is printed before    *
 * and after the matrix.                                          *
 *                                                                *
 ******************************************************************/
void dense_print2(double **a, int64 n);

#endif
