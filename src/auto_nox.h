#ifndef AUTO_NOX_H
#define AUTO_NOX_H
#include "integers.h"

#include <stdio.h>
#include "autlim.h"
#define MAX_AUT_PER 10

extern int32 HomoFlag;

typedef struct {
    int32 irot;
    int32 nrot[1000];
    double torper;
} Rotchk;

typedef struct {
    int32 exist;
    int32 ntst;
    int32 nmx;
    int32 npr;
    double ds, dsmax, dsmin, rl0, rl1, a0, a1;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double lastx;
    double lasty;
    int32 wid;
    int32 hgt;
    int32 x0;
    int32 y0;
    int32 st_wid;
    int32 nfpar;
    int32 nbc;
    int32 ips;
    int32 irs;
    int32 ilp;
    int32 isp;
    int32 isw;
    int32 itp;
    int32 plot;
    int32 var;
    int32 icp1;
    int32 icp2;
    int32 icp3;
    int32 icp4;
    int32 icp5;
    int32 nper;
    char hinttxt[256];
    double period[MAX_AUT_PER];
    int64 uzrpar[MAX_AUT_PER];
    double epsl;
    double epsu;
    double epss;
    int32 ncol;
} Bifurcation;

typedef struct {
    int32 iad;
    int32 mxbf;
    int32 iid;
    int32 itmx;
    int32 itnw;
    int32 nwtn;
    int32 iads;
} AdvAuto;

typedef struct {
    int32 package;
    int32 ibr;
    int32 ntot;
    int32 itp;
    int32 lab;
    double norm, uhi[NAUTO], ulo[NAUTO], u0[NAUTO], ubar[NAUTO];
    double par[20], per, torper;
    int32 index;
    int32 nfpar;
    int32 icp1;
    int32 icp2;
    int32 icp3;
    int32 icp4;
    int32 icp5;
    int32 flag;
} GrabPoint;

typedef struct diagram {
    int32 package;
    int32 ibr;
    int32 ntot;
    int32 itp;
    int32 lab;
    int32 calc;
    double norm, *uhi, *ulo, *u0, *ubar, *evr, *evi;
    double par[20], per, torper;
    int32 index;
    int32 nfpar;
    int32 icp1;
    int32 icp2;
    int32 icp3;
    int32 icp4;
    int32 icp5;
    int32 flag2;
    struct diagram *prev;
    struct diagram *next;
} Diagram;

typedef struct {
    int32 plot;
    int32 var;
    int32 icp1;
    int32 icp2;
    int32 icp3;
    int32 icp4;
    int32 icp5;
    double xmin, ymin, xmax, ymax;
} AutoAX;

extern GrabPoint grabpt;
extern Rotchk blrtn;
extern Bifurcation Auto;
extern AdvAuto aauto;

extern int32 NewPeriodFlag;

extern double homo_l[100];
extern double homo_r[100];

extern int32 TypeOfCalc;
extern int32 SEc;
extern int32 UEc;
extern int32 SPc;
extern int32 UPc;

extern int32 RestartLabel;

extern int32 auto_ntst;
extern int32 auto_nmx;
extern int32 auto_npr;
extern int32 auto_ncol;
extern double auto_ds;
extern double auto_dsmax;
extern double auto_dsmin;
extern double auto_rl0;
extern double auto_rl1;
extern double auto_a0;
extern double auto_a1;
extern double auto_xmax;
extern double auto_xmin;
extern double auto_ymax;
extern double auto_ymin;
extern double auto_epsl;
extern double auto_epsu;
extern double auto_epss;
extern int32 auto_var;

extern int32 Auto_index_to_array[8];
extern int32 AutoPar[8];
extern double outperiod[20];
extern int64 UzrPar[20];
extern int32 NAutoUzr;

extern char fort3[200];
extern char fort7[200];
extern char fort8[200];
extern char fort9[200];

extern int32 load_all_labeled_orbits;

extern int32 AutoTwoParam;
extern int32 NAutoPar;

void auto_nox_colset(int32 type);
void auto_nox_colset2(int32 flag2);
void auto_nox_get_str(char *xlabel, char *ylabel);
void auto_nox_draw_ps_axes(void);
void auto_nox_draw_svg_axes(void);
void auto_nox_draw_bix_axes(void);
int32 auto_nox_ix_val(double x);
int32 auto_nox_iy_val(double y);
int32 auto_nox_check_bnds(int32 ix, int32 iy);
void auto_nox_renamef(char *old, char *new);
void auto_nox_copyf(char *old, char *new);
void auto_nox_appendf(char *old, char *new);
void auto_nox_deletef(char *old);
void auto_nox_close(int32 flag);
void auto_nox_open(int32 flag);
void auto_nox_do(int32 iold, int32 isave);
void auto_nox_set_auto(void);
int32 auto_nox_name_to_index(char *s);
int32 auto_nox_par_to_name(int64 index, char *s);
void auto_nox_per_par(void);
void auto_nox_params(void);
void auto_nox_num_par(void);
void auto_nox_plot_par(void);
void auto_fit(void);
void auto_default(void);
void auto_nox_zoom_in(int32 i1, int32 j1, int32 i2, int32 j2);
void auto_nox_zoom_out(int32 i1, int32 j1, int32 i2, int32 j2);
void auto_nox_xy_plot(double *x, double *y1, double *y2, double par1,
                      double par2, double per, double *uhigh, double *ulow,
                      double *ubar, double a);
int32 auto_nox_plot_point(int32 flag2, int32 icp1, int32 icp2);
void auto_nox_add_ps_point(double *par, double per, double *uhigh, double *ulow,
                           double *ubar, double a, int32 type, int32 flag,
                           int32 icp1, int32 icp2, int32 flag2);
void auto_nox_line(double x1i, double y1i, double x2i, double y2i);
void auto_nox_add_point(double *par, double per, double *uhigh, double *ulow,
                        double *ubar, double a, int32 type, int32 flg,
                        int32 lab, int32 icp1, int32 icp2, int32 flag2,
                        double *evr, double *evi);
void auto_nox_get_bif_sym(char *at, int32 itp);
void auto_nox_info_header(int32 icp1, int32 icp2);
void auto_nox_new_info(int32 ibr, int32 pt, char *ty, int32 lab, double *par,
                       double norm, double u0, double per, int32 icp1,
                       int32 icp2);
void auto_nox_traverse_out(Diagram *d, int32 *ix, int32 *iy, int32 dodraw);
void auto_nox_win(void);
void auto_nox_load_last_plot(int32 flag);
void auto_nox_keep_last_plot(int32 flag);
void auto_nox_init_win(void);
void auto_nox_plot_stab(double *evr, double *evi, int32 n);
int32 auto_nox_yes_reset(void);
int32 auto_nox_reset(void);
void auto_nox_grab(void);

void auto_start_diff_ss(void);
void auto_start_at_bvp(void);
void auto_start_at_per(void);
void auto_nox_find_best_homo_shift(int32 n);
void auto_nox_get_start_period(double *p);
void auto_nox_get_start_orbit(double *u, double t, int32 n);
void auto_nox_get_shifted_orbit(double *u, double t, double p, int32 n);
void auto_nox_new_ss(void);
void auto_new_discrete(void);
void auto_nox_extend_ss(void);
void auto_start_choice(void);
void auto_nox_torus_choice(void);
void auto_nox_per_doub_choice(void);
void auto_nox_periodic_choice(void);
void auto_nox_hopf_choice(void);
void auto_nox_new_per(void);
int32 auto_nox_get_homo_info(int32 *nun, int32 *nst, double *ul, double *ur);
void auto_extend_homoclinic(void);
void auto_nox_2p_limit(int32 ips);
void auto_nox_2p_branch(int32 ips);
void auto_nox_branch_choice(int32 ibr, int32 ips);
void auto_nox_homo_choice(int32 itp);
void auto_nox_2p_fixper(void);
void auto_nox_err(char *s);
void auto_nox_run(void);
void auto_nox_load_orbit(void);
void auto_nox_save(void);
void auto_nox_save_numerics(FILE *fp);
void auto_nox_load_numerics(FILE *fp);
void auto_nox_save_graph(FILE *fp);
void auto_nox_load_graph(FILE *fp);
void auto_nox_save_q_file(FILE *fp);
void auto_nox_make_q_file(FILE *fp);
int32 auto_nox_no_info_noinfo(char *s);
void auto_nox_load(void);
int32 auto_nox_move_to_label(int32 mylab, int32 *nrow, int32 *ndim, FILE *fp);
void auto_nox_get_a_row(double *u, double *t, int32 n, FILE *fp);
void auto_nox_file(void);
int32 auto_nox_check_plot_type(int32 flag2, int32 icp1, int32 icp2);
void auto_nox_store_point(double x, double y);

#endif
