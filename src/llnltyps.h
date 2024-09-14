/******************************************************************
 *                                                                *
 * File          : llnltyps.h                                     *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Last Modified : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This header file exports three types: double, int64, and bool  *
 * (short for boolean), as well as the constants TRUE and FALSE.  *
 *                                                                *
 * Users should #include "llnltyps.h" in any file that should     *
 * be easily modifiable to work with different double or int64    *
 * types and use the exported names double and int64 within such  *
 * a file. The types for double and int64 below have been set to  *
 * double and int32, respectively. A user should modify these       *
 * type declarations as he/she sees fit. For example, if a user   *
 * wants the work with type float because double precision        *
 * floating point arithmetic is too expensive on the user's       *
 * machine, then the definition below should be changed to:       *
 *                                                                *
 * typedef float double;                                            *
 *                                                                *
 * Similarly, if a user needs to work with extremely large        *
 * integers (see the system header file <limits.h> for the limits *
 * on type int32 and   int32 on your machine), then the user       *
 * should change the definition below to:                         *
 *                                                                *
 * typedef   int32 int64;                                      *
 *                                                                *
 * The constants FLOAT, DOUBLE, INT, LONG_INT indicate the        *
 * underlying types for double and int64. They should be set as   *
 * follows:                                                       *
 *                                                                *
 * (1) #define FLOAT 1                                            *
 *     #define DOUBLE 0     (double is float)                       *
 *                                                                *
 * (2) #define FLOAT 0                                            *
 *     #define DOUBLE 1     (double is double)                      *
 *                                                                *
 * (3) #define INT 1                                              *
 *     #define LONG_INT 0   (int64 is int32)                      *
 *                                                                *
 * (4) #define INT 0                                              *
 *     #define LONG_INT 1   (int64 is   int32)                 *
 *                                                                *
 * Thus the legal types for double are float and double, while      *
 * the legal types for int64 are int32 and   int32. The macro    *
 * RCONST gives a user a convenient way to define double            *
 * constants. To use the double constant 1.0, for example, the      *
 * user should write                                              *
 *                                                                *
 * #define ONE RCONST(1.0)                                        *
 *                                                                *
 * If double is double, then RCONST(1.0) expands to 1.0. If double is *
 * float, then RCONST(1.0) expands to 1.0F. There is never a      *
 * need to explicitely cast 1.0 to (double).                        *
 *                                                                *
 ******************************************************************/

#ifndef llnltyps_h
#define llnltyps_h

/******************************************************************
 *                                                                *
 * Types : double, int64                                          *
 *----------------------------------------------------------------*
 * The types double and int64 are currently set to double and     *
 * int32, respectively. See the documentation at the top for        *
 * usage details and a description of associated constants and    *
 * macros.                                                        *
 *                                                                *
 ******************************************************************/
#include <stdint.h>
#include "integers.h"

#define FLOAT 0
#define DOUBLE 1

#define INT 1
#define LONG_INT 0

#if FLOAT

#define RCONST(x) x##F

#elif DOUBLE

#define RCONST(x) x

#endif

/******************************************************************
 *                                                                *
 * Type : bool                                                    *
 * Constants : FALSE, TRUE                                        *
 *----------------------------------------------------------------*
 * ANSI C does not have a built-in boolean type. Below is the     *
 * definition for a new type bool. The advantage of using the     *
 * name bool (instead of int32) is an increase in code readability. *
 * It allows the programmer to make a distinction between int32 and *
 * boolean data. Variables of type bool are intended to have only *
 * the two values FALSE and TRUE which are defined below to be    *
 * equal to 0 and 1, respectively.                                *
 *                                                                *
 ******************************************************************/

#ifndef bool
#define bool int32
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#endif
