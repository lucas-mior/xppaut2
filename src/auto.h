#ifndef _auto_h_
#define _auto_h_
#include "integers.h"

#include <X11/Xlib.h>
#include <stdio.h>

int32 get_auto_str(char *xlabel, char *ylabel);
int32 draw_ps_axes(void);
int32 draw_bif_axes(void);
int32 byeauto_(int32 *iflag);
int32 IXVal(double x);
int32 IYVal(double y);
int32 Circle(int32 x, int32 y, int32 r);
int32 XORCross(int32 x, int32 y);
int32 FillCircle(int32 x, int32 y, int32 r);
int32 LineWidth(int32 wid);
int32 renamef(char *old, char *new);
int32 copyf(char *old, char *new);
int32 appendf(char *old, char *new);
int32 deletef(char *old);
int32 close_auto(int32 flag);
int32 open_auto(int32 flag);
int32 do_auto(int32 iold, int32 isave);
int32 set_auto(void);
int32 auto_name_to_index(char *s);
int32 auto_par_to_name(int32 index, char *s);
int32 auto_per_par(void);
int32 auto_params(void);
int32 auto_num_par(void);
int32 auto_plot_par(void);
int32 auto_fit(void);
int32 auto_zoom(int32 i1, int32 j1, int32 i2, int32 j2);
int32 auto_xy_plot(double *x, double *y1, double *y2, double par1, double par2,
                   double per, double *uhigh, double *ulow, double *ubar,
                   double a);
int32 plot_point(int32 flag2, int32 icp1, int32 icp2);
int32 add_ps_point(double *par, double per, double *uhigh, double *ulow,
                   double *ubar, double a, int32 type, int32 flag, int32 lab,
                   int32 npar, int32 icp1, int32 icp2, int32 flag2, double *evr,
                   double *evi);
int32 add_point(double *par, double per, double *uhigh, double *ulow,
                double *ubar, double a, int32 type, int32 flag, int32 lab,
                int32 npar, int32 icp1, int32 icp2, int32 flag2, double *evr,
                double *evi);
int32 redraw_auto_menus(void);
int32 get_bif_sym(char *at, int32 itp);
int32 info_header(int32 icp1, int32 icp2);
void new_info(int32 ibr, int32 pt, char *ty, int32 lab, double *par,
              double norm, double u0, double per,
              int32 icp1, int32 icp2);
int32 traverse_diagram(void);
int32 clear_auto_plot(void);
int32 do_auto_win(void);
int32 load_last_plot(int32 flag);
int32 keep_last_plot(int32 flag);
int32 init_auto_win(void);
int32 plot_stab(double *evr, double *evi, int32 n);
int32 clr_stab(void);
int32 auto_motion(XEvent ev);
int32 display_auto(Window w);
Window lil_button(Window root, int32 x, int32 y, char *name);
int32 make_auto(char *wname, char *iname);
int32 yes_reset_auto(void);
int32 reset_auto(void);
int32 auto_grab(void);
int32 auto_start_diff_ss(void);
int32 auto_start_at_bvp(void);
int32 auto_start_at_per(void);
int32 get_start_period(double *p);
void get_start_orbit(double *u, double t, double p, int32 n);
int32 auto_new_ss(void);
int32 auto_new_discrete(void);
int32 auto_extend_ss(void);
int32 auto_start_choice(void);
int32 torus_choice(void);
int32 per_doub_choice(void);
int32 periodic_choice(void);
int32 hopf_choice(void);
int32 auto_new_per(void);
int32 auto_extend_bvp(void);
int32 auto_switch_per(void);
int32 auto_switch_bvp(void);
int32 auto_switch_ss(void);
int32 auto_2p_limit(int32 ips);
int32 auto_twopar_double(void);
int32 auto_torus(void);
int32 auto_2p_branch(void);
int32 auto_branch_choice(int32 ibr);
int32 auto_2p_fixper(void);
int32 auto_2p_hopf(void);
int32 auto_period_double(void);
int32 auto_err(char *s);
int32 auto_run(void);
int32 load_auto_orbit(void);
int32 save_auto(void);
int32 save_auto_numerics(FILE *fp);
int32 load_auto_numerics(FILE *fp);
int32 save_auto_graph(FILE *fp);
int32 load_auto_graph(FILE *fp);
int32 save_q_file(FILE *fp);
int32 make_q_file(FILE *fp);
int32 noinfo(char *s);
int32 load_auto(void);
int32 move_to_label(int32 mylab, int32 *nrow, int32 *ndim, FILE *fp);
int32 get_a_row(double *u, double *t, int32 n, FILE *fp);
int32 auto_file(void);
int32 a_msg(int32 i, int32 v);
int32 auto_enter(Window w, int32 v);
int32 auto_button(XEvent ev);
int32 auto_kill(void);
int32 auto_keypress(XEvent ev, int32 *used);

#endif
