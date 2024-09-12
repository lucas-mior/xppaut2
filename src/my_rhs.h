#ifndef _my_rhs_h_
#define _my_rhs_h_
#include "integers.h"

int32 MAIN__(void);
int32 main(int32 argc, char **argv);
void extra(double *y__y, double t, int32 nod, int32 neq);
void set_fix_rhs(double t, double *y);
int32 my_rhs(double t, double *y, double *ydot, int32 neq);
void update_based_on_current(void);
void fix_only(void);
void rhs_only(double *ydot);
void vec_rhs(double t, double *y, double *ydot, int32 neq);

#endif
