/******************************************************************
 *                                                                *
 * File          : llnlmath.h                                     *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the header file for a C math library. The routines     *
 * listed here work with the type double as defined in llnltyps.h.  *
 * To do single precision floating point arithmetic, set the type *
 * double to be float. To do double precision arithmetic, set the   *
 * type double to be double. The default implementations for        *
 * RPowerR and RSqrt call standard math library functions which   *
 * do double precision arithmetic. If this is unacceptable when   *
 * double is float, then the user should re-implement these two     *
 * routines by calling single precision routines available on     *
 * his/her machine.                                               *
 *                                                                *
 ******************************************************************/

#ifndef _llnlmath_h
#define _llnlmath_h

#include "integers.h"
#include "llnltyps.h"

/******************************************************************
 *                                                                *
 * Macros : MIN, MAX, ABS, SQR                                    *
 *----------------------------------------------------------------*
 * MIN(A, B) returns the minimum of A and B.                      *
 *                                                                *
 * MAX(A, B) returns the maximum of A and B.                      *
 *                                                                *
 * ABS(A) returns the absolute value of A.                        *
 *                                                                *
 * SQR(A) returns the square of A.                                *
 *                                                                *
 ******************************************************************/

#define MIN(A, B) ((A) < (B) ? (A) : (B))

#define MAX(A, B) ((A) > (B) ? (A) : (B))

#define ABS(A) ((A > 0) ? (A) : -(A))

#define SQR(A) ((A) * (A))

/******************************************************************
 *                                                                *
 * Function : UnitRoundoff                                        *
 * Usage    : double uround;                                        *
 *            uround = UnitRoundoff();                            *
 *----------------------------------------------------------------*
 * UnitRoundoff returns the unit roundoff u for double floating     *
 * point arithmetic, where u is defined to be the largest         *
 * positive double such that 1.0 + u != 1.0.                        *
 *                                                                *
 ******************************************************************/

double UnitRoundoff(void);

/******************************************************************
 *                                                                *
 * Function : RPowerI                                             *
 * Usage    : int32 exponent;                                       *
 *            double base, ans;                                     *
 *            ans = RPowerI(base,exponent);                       *
 *----------------------------------------------------------------*
 * RPowerI returns the value base^exponent, where base is a double  *
 * and exponent is an int32.                                        *
 *                                                                *
 ******************************************************************/

double RPowerI(double base, int32 exponent);

/******************************************************************
 *                                                                *
 * Function : RPowerR                                             *
 * Usage    : double base, exponent, ans;                           *
 *            ans = RPowerI(base,exponent);                       *
 *----------------------------------------------------------------*
 * RPowerR returns the value base^exponent, where both base and   *
 * exponent are reals. If base < 0.0, then RPowerR returns 0.0.   *
 *                                                                *
 ******************************************************************/

double RPowerR(double base, double exponent);

/******************************************************************
 *                                                                *
 * Function : RSqrt                                               *
 * Usage    : double sqrt_x;                                        *
 *            sqrt_x = RSqrt(x);                                  *
 *----------------------------------------------------------------*
 * RSqrt(x) returns the square root of x. If x < 0.0, then RSqrt  *
 * returns 0.0.                                                   *
 *                                                                *
 ******************************************************************/

double RSqrt(double x);

#endif
