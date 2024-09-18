#ifndef auto_nox_h
#define auto_nox_h
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
    int32 ntst, nmx, npr;
    double ds, dsmax, dsmin, rl0, rl1, a0, a1;
    double xmin, xmax, ymin, ymax;
    double lastx, lasty;
    int32 wid, hgt, x0, y0, st_wid;
    int32 nfpar, nbc;
    int32 ips, irs, ilp, isp, isw, itp;
    int32 plot, var;
    int32 icp1, icp2, icp3, icp4, icp5;
    int32 nper;
    char hinttxt[256];
    double period[MAX_AUT_PER];
    int64 uzrpar[MAX_AUT_PER];
    double epsl, epsu, epss;
    int32 ncol;
} BIFUR;

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
    int32 ibr, ntot, itp, lab;
    double norm, uhi[NAUTO], ulo[NAUTO], u0[NAUTO], ubar[NAUTO];
    double par[20], per, torper;
    int32 index, nfpar, icp1, icp2, icp3, icp4, icp5;
    int32 flag;
} GrabPoint;

typedef struct diagram {
    int32 package;
    int32 ibr, ntot, itp, lab, calc;
    double norm, *uhi, *ulo, *u0, *ubar, *evr, *evi;
    double par[20], per, torper;
    int32 index, nfpar;
    int32 icp1, icp2, icp3, icp4, icp5, flag2;
    struct diagram *prev;
    struct diagram *next;
} Diagram;

typedef struct {
    int32 plot, var, icp1, icp2, icp3, icp4, icp5;
    double xmin, ymin, xmax, ymax;
} AUTOAX;

extern GrabPoint grabpt;
extern Rotchk blrtn;

extern int32 TypeOfCalc;
extern int32 SEc;
extern int32 UEc;
extern int32 SPc;
extern int32 UPc;

extern int32 RestartLabel;

extern int32 auto_ntst, auto_nmx, auto_npr, auto_ncol;
extern double auto_ds, auto_dsmax, auto_dsmin;
extern double auto_rl0, auto_rl1, auto_a0, auto_a1;
extern double auto_xmax, auto_xmin, auto_ymax, auto_ymin;
extern double auto_epsl, auto_epsu, auto_epss;
extern int32 auto_var;

extern int32 Auto_index_to_array[8];
extern int32 AutoPar[8];
extern double outperiod[20];

extern int32 load_all_labeled_orbits;

extern int32 AutoTwoParam;
extern int32 NAutoPar;

void colset(int32 type);
void pscolset2(int32 flag2);
void colset2(int32 flag2);
void get_auto_str(char *xlabel, char *ylabel);
void draw_ps_axes(void);
void draw_svg_axes(void);
void draw_bif_axes(void);
int32 IXVal(double x);
int32 IYVal(double y);
int32 chk_auto_bnds(int32 ix, int32 iy);
void renamef(char *old, char *new);
void copyf(char *old, char *new);
void appendf(char *old, char *new);
void deletef(char *old);
void close_auto(int32 flag);
void open_auto(int32 flag);
void do_auto(int32 iold, int32 isave);
void set_auto(void);
int32 auto_name_to_index(char *s);
int32 auto_par_to_name(int64 index, char *s);
void auto_per_par(void);
void auto_params(void);
void auto_num_par(void);
void auto_plot_par(void);
void auto_fit(void);
void auto_default(void);
void auto_zoom_in(int32 i1, int32 j1, int32 i2, int32 j2);
void auto_zoom_out(int32 i1, int32 j1, int32 i2, int32 j2);
void auto_xy_plot(double *x, double *y1, double *y2, double par1, double par2,
                  double per, double *uhigh, double *ulow, double *ubar,
                  double a);
int32 plot_point(int32 flag2, int32 icp1, int32 icp2);
void add_ps_point(double *par, double per, double *uhigh, double *ulow,
                  double *ubar, double a, int32 type, int32 flag, int32 icp1,
                  int32 icp2, int32 flag2);
void auto_line(double x1i, double y1i, double x2i, double y2i);
void add_point(double *par, double per, double *uhigh, double *ulow,
               double *ubar, double a, int32 type, int32 flg, int32 lab,
               int32 icp1, int32 icp2, int32 flag2, double *evr, double *evi);
void get_bif_sym(char *at, int32 itp);
void info_header(int32 icp1, int32 icp2);
void new_info(int32 ibr, int32 pt, char *ty, int32 lab, double *par,
              double norm, double u0, double per, int32 icp1, int32 icp2);
void traverse_out(Diagram *d, int32 *ix, int32 *iy, int32 dodraw);
void do_auto_win(void);
void load_last_plot(int32 flag);
void keep_last_plot(int32 flag);
void init_auto_win(void);
void plot_stab(double *evr, double *evi, int32 n);
int32 yes_reset_auto(void);
int32 reset_auto(void);
void auto_grab(void);

void auto_start_diff_ss(void);
void auto_start_at_bvp(void);
void auto_start_at_per(void);
void find_best_homo_shift(int32 n);
void get_start_period(double *p);
void get_start_orbit(double *u, double t, int32 n);
void get_shifted_orbit(double *u, double t, double p, int32 n);
void auto_new_ss(void);
void auto_new_discrete(void);
void auto_extend_ss(void);
void auto_start_choice(void);
void torus_choice(void);
void per_doub_choice(void);
void periodic_choice(void);
void hopf_choice(void);
void auto_new_per(void);
void auto_start_at_homoclinic(void);
int32 get_homo_info(int32 *nun, int32 *nst, double *ul, double *ur);
void auto_extend_homoclinic(void);
void auto_extend_bvp(void);
void auto_switch_per(void);
void auto_switch_bvp(void);
void auto_switch_ss(void);
void auto_2p_limit(int32 ips);
void auto_twopar_double(void);
void auto_torus(void);
void auto_2p_branch(int32 ips);
void auto_branch_choice(int32 ibr, int32 ips);
void auto_homo_choice(int32 itp);
void auto_2p_fixper(void);
void auto_2p_hopf(void);
void auto_period_double(void);
void auto_err(char *s);
void auto_run(void);
void load_auto_orbit(void);
void load_auto_orbitx(int32 ibr, int32 flag, int32 lab, double per);
void save_auto(void);
void save_auto_numerics(FILE *fp);
void load_auto_numerics(FILE *fp);
void save_auto_graph(FILE *fp);
void load_auto_graph(FILE *fp);
void save_q_file(FILE *fp);
void make_q_file(FILE *fp);
int32 noinfo(char *s);
void load_auto(void);
int32 move_to_label(int32 mylab, int32 *nrow, int32 *ndim, FILE *fp);
void get_a_row(double *u, double *t, int32 n, FILE *fp);
void auto_file(void);
int32 check_plot_type(int32 flag2, int32 icp1, int32 icp2);
void storeautopoint(double x, double y);

#endif
