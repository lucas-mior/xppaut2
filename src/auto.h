#ifndef auto_h
#define auto_h
#include "integers.h"

#include <X11/Xlib.h>
#include <stdio.h>

int32 auto_nox_get_str(char *xlabel, char *ylabel);
int32 auto_nox_draw_ps_axes(void);
int32 auto_nox_draw_bix_axes(void);
int32 auto_x11_bye(int32 *iflag);
int32 auto_nox_ix_val(double x);
int32 auto_nox_iy_val(double y);
int32 auto_x11_circle(int32 x, int32 y, int32 r);
int32 auto_x11_xor_cross(int32 x, int32 y);
int32 auto_x11_fill_circle(int32 x, int32 y, int32 r);
int32 auto_x11_line_width(int32 wid);
int32 auto_nox_renamef(char *old, char *new);
int32 auto_nox_copyf(char *old, char *new);
int32 auto_nox_deletef(char *old);
int32 auto_close(int32 flag);
int32 auto_nox_open(int32 flag);
int32 auto_nox_do(int32 iold, int32 isave);
int32 auto_nox_set_auto(void);
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
int32 auto_nox_plot_point(int32 flag2, int32 icp1, int32 icp2);
int32 auto_add_ps_point(double *par, double per, double *uhigh, double *ulow,
                        double *ubar, double a, int32 type, int32 flag,
                        int32 lab, int32 npar, int32 icp1, int32 icp2,
                        int32 flag2, double *evr, double *evi);
int32 auto_add_point(double *par, double per, double *uhigh, double *ulow,
                     double *ubar, double a, int32 type, int32 flag, int32 lab,
                     int32 npar, int32 icp1, int32 icp2, int32 flag2,
                     double *evr, double *evi);
int32 auto_x11_redraw_menus(void);
int32 auto_nox_get_bif_sym(char *at, int32 itp);
int32 auto_nox_info_header(int32 icp1, int32 icp2);
void auto_nox_new_info(int32 ibr, int32 pt, char *ty, int32 lab, double *par,
                       double norm, double u0, double per, int32 icp1,
                       int32 icp2);
int32 auto_x11_traverse_diagram(void);
int32 auto_x11_clear_plot(void);
int32 auto_nox_win(void);
int32 auto_nox_load_last_plot(int32 flag);
int32 auto_nox_keep_last_plot(int32 flag);
int32 auto_nox_init_win(void);
int32 auto_nox_plot_stab(double *evr, double *evi, int32 n);
int32 auto_x11_clr_stab(void);
int32 auto_x11_motion(XEvent event);
int32 auto_x11_display(Window window);
Window auto_x11_lil_button(Window root, int32 x, int32 y);
int32 auto_x11_make(char *wname, char *iname);
int32 auto_nox_yes_reset(void);
int32 auto_nox_reset(void);
int32 auto_grab(void);
int32 auto_start_diff_ss(void);
int32 auto_start_at_bvp(void);
int32 auto_start_at_per(void);
int32 auto_nox_get_start_period(double *p);
void auto_nox_get_start_orbit(double *u, double t, int32 n);
int32 auto_new_ss(void);
int32 auto_new_discrete(void);
int32 auto_extend_ss(void);
int32 auto_start_choice(void);
int32 auto_nox_torus_choice(void);
int32 auto_nox_per_doub_choice(void);
int32 auto_nox_periodic_choice(void);
int32 auto_nox_hopf_choice(void);
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
int32 auto_nox_load_orbit(void);
int32 auto_nox_save(void);
int32 auto_nox_save_numerics(FILE *fp);
int32 auto_nox_load_numerics(FILE *fp);
int32 auto_nox_save_graph(FILE *fp);
int32 auto_nox_load_graph(FILE *fp);
int32 auto_nox_save_q_file(FILE *fp);
int32 auto_nox_make_q_file(FILE *fp);
int32 auto_nox_no_info_noinfo(char *s);
int32 auto_nox_load(void);
int32 auto_nox_move_to_label(int32 mylab, int32 *nrow, int32 *ndim, FILE *fp);
int32 auto_nox_get_a_row(double *u, double *t, int32 n, FILE *fp);
int32 auto_file(void);
int32 auto_x11_enter(Window window, int32 v);
int32 auto_x11_button(XEvent event);
int32 auto_x11_kill(void);
int32 auto_x11_keypress(XEvent event, int32 *used);

#endif
