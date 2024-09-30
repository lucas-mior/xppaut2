/* f2c.h  --  Standard Fortran to C header file */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

        - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include "integers.h"
#include <complex_math.h>

#define min(a, b) ((a) <= (b) ? (a) : (b))
#define max(a, b) ((a) >= (b) ? (a) : (b))

#define ARRAY2D(array, i, j) array[(i) + (j)*array##_dim1]
#define ARRAY3D(array, i, j, k) array[(i) + ((j) + (k)*array##_dim2)*array##_dim1]

#endif
