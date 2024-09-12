#ifndef _integrate_h_
#define _integrate_h_
#include "integers.h"

#include <stdio.h>

void init_ar_ic(void);
void dump_range(FILE *fp, int32 f);
void init_range(void);
int32 set_up_eq_range(void);
void cont_integ(void);
int32 range_item(void);
int32 range_item2(void);
int32 set_up_range(void);
int32 set_up_range2(void);
void init_monte_carlo(void);
void monte_carlo(void);
void do_monte_carlo_search(int32 append, int32 stuffbrowse, int32 ishoot);
void do_eq_range(double *x);
void swap_color(int32 *col, int32 rorw);
void set_cycle(int32 flag, int32 *icol);
int32 do_range(double *x, int32 flag);
void find_equilib_com(int32 com);
void batch_integrate(void);
void do_batch_dry_run(void);
void batch_integrate_once(void);
int32 write_this_run(char *file, int32 i);
void do_init_data(int32 com);
void run_from_x(double *x);
void run_now(void);
void do_start_flags(double *x, double *t);
void usual_integrate_stuff(double *x);
void do_new_array_ic(char *new, int32 j1, int32 j2);
void store_new_array_ic(char *new, int32 j1, int32 j2, char *formula);
void evaluate_ar_ic(char *v, char *f, int32 j1, int32 j2);
int32 extract_ic_data(char *big);
void arr_ic_start(void);
int32 set_array_ic(void);
int32 form_ic(void);
void get_ic(int32 it, double *x);
int32 ode_int(double *y, double *t, int32 *istart, int32 ishow);
int32 integrate(double *t, double *x, double tend, double dt, int32 count, int32 nout,
              int32 *start);
void send_halt(double *y, double t);
void send_output(double *y, double t);
void do_plot(float *oldxpl, float *oldypl, float *oldzpl, float *xpl,
             float *ypl, float *zpl);
void export_data(FILE *fp);
void plot_the_graphs(float *xv, float *xvold, int32 node, int32 neq, double ddt,
                     int32 *tc, int32 flag);
void plot_one_graph(float *xv, float *xvold, int32 node, int32 neq, double ddt,
                    int32 *tc);
void restore(int32 i1, int32 i2);
void comp_color(float *v1, float *v2, int32 n, double dt);
void shoot(double *x, double *xg, double *evec, int32 sgn);
void shoot_easy(double *x);
void stop_integration(void);
int32 stor_full(void);
int32 do_auto_range_go();

#endif
