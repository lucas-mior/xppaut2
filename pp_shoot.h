#ifndef _pp_shoot_h_
#define _pp_shoot_h_
#include "integers.h"

#include <stdio.h>

void do_bc(double *y__0, double t0, double *y__1, double t1, double *f, int32 n);
void compile_bvp(void);
void reset_bvp(void);
void init_shoot_range(char *s);
void dump_shoot_range(FILE *fp, int32 f);
void bad_shoot(int32 iret);
void do_sh_range(double *ystart, double *yend);
int32 set_up_homoclinic(void);
int32 set_up_periodic(int32 *ipar, int32 *ivar, double *sect, int32 *ishow);
void find_bvp_com(int32 com);
void last_shot(int32 flag);
int32 set_up_sh_range(void);
void bvshoot(double *y, double *yend, double err, double eps, int32 maxit,
             int32 *iret, int32 n, int32 ishow, int32 iper, int32 ipar, int32 ivar,
             double sect);

#endif
