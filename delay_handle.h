
#ifndef _delay_handle_h_
#define _delay_handle_h_
#include "integers.h"

double delay_stab_eval(double delay, int32 var);
int32 alloc_delay(double big);
void free_delay(void);
void stor_delay(double *y);
double get_delay_old(int32 in, double tau);
void polint(double *xa, double *ya, int32 n, double x, double *y, double *dy);
double get_delay(int32 in, double tau);
int32 do_init_delay(double big);

#endif
