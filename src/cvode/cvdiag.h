/******************************************************************
 *                                                                *
 * File          : cvdiag.h                                       *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the header file for the CVODE diagonal linear solver,  *
 * CVDIAG.                                                        *
 *                                                                *
 * Note: The type int64 must be large enough to store the value *
 * of the linear system size N.                                   *
 *                                                                *
 ******************************************************************/

#ifndef _cvdiag_h
#define _cvdiag_h

#include <stdio.h>
#include "cvode.h"
#include "vector.h"

/******************************************************************
 *                                                                *
 * CVDIAG solver statistics indices                               *
 *----------------------------------------------------------------*
 * The following enumeration gives a symbolic name to each        *
 * CVDIAG statistic. The symbolic names are used as indices into  *
 * the iopt and ropt arrays passed to cvode_malloc.                *
 * The CVDIAG statistics are:                                     *
 *                                                                *
 * iopt[DIAG_LRW] : size (in double words) of double workspace        *
 *                  vectors used by this solver.                  *
 *                                                                *
 * iopt[DIAG_LIW] : size (in int64 words) of int64            *
 *                  workspace vectors used by this solver.        *
 *                                                                *
 * The number of diagonal approximate Jacobians formed is equal   *
 * to the number of cv_diag_setup calls. This number is available   *
 * in cv_iopt[NSETUPS].                                           *
 *                                                                *
 ******************************************************************/

enum {
    DIAG_LRW = CVODE_IOPT_SIZE,
    DIAG_LIW
};

/******************************************************************
 *                                                                *
 * Function : CVDiag                                              *
 *----------------------------------------------------------------*
 * A call to the CVDiag function links the main CVODE integrator  *
 * with the CVDIAG linear solver.                                 *
 *                                                                *
 * cvode_mem is the pointer to CVODE memory returned by           *
 *              cvode_malloc.                                      *
 *                                                                *
 ******************************************************************/

void cv_diag(void *cvode_mem);

#endif
