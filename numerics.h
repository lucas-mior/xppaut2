#ifndef _numerics_h_
#define _numerics_h_
#include "integers.h"

/*       Numerics.h   */

extern double DELTA_T, TEND, T0, TRANS, NULL_ERR, EVEC_ERR, NEWT_ERR;
extern double BOUND, DELAY, TOLER, HMIN, HMAX;
extern float *fft_data, *hist_data, color_scale, min_scale;
extern double POIPLN;

extern int32 NMESH, NJMP, METHOD, color_flag, NC_ITER;
extern int32 EVEC_ITER, FOREVER;

extern int32 POIMAP, POIVAR, POISGN, SOS;

extern int32 HIST, HVAR, hist_ind;

extern int32 XSHFT, YSHFT, ZSHFT;

void chk_volterra(void);
void check_pos(int32 *j);
void quick_num(int32 com);
void get_num_par(int32 ch);
void chk_delay(void);
void set_delay(void);
void ruelle(void);
void init_numerics(void);
void meth_dialog(void);
void get_pmap_pars_com(int32 l);
void get_method(void);
void set_col_par_com(int32 i);
void do_meth(void);
void set_total(double total);
void user_set_color_par(int32 flag, char *via, double lo, double hi);
void compute_one_period(double period, double *x, char *name);
#endif
