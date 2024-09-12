#ifndef _diagram_h_
#define _diagram_h_
#include "integers.h"

#include <stdio.h>
#include "auto_nox.h"

void start_diagram(int32 n);
int32 find_diagram(int32 irs, int32 n, int32 *index, int32 *ibr, int32 *ntot, int32 *itp,
                 int32 *nfpar, double *a, double *uhi, double *ulo, double *u0,
                 double *par, double *per, int32 *icp1, int32 *icp2, int32 *icp3,
                 int32 *icp4);
void edit_start(int32 ibr, int32 ntot, int32 itp, int32 lab, int32 nfpar, double a,
                double *uhi, double *ulo, double *u0, double *ubar, double *par,
                double per, int32 n, int32 icp1, int32 icp2, int32 icp3, int32 icp4,
                double *evr, double *evi);
void edit_diagram(DIAGRAM *d, int32 ibr, int32 ntot, int32 itp, int32 lab, int32 nfpar,
                  double a, double *uhi, double *ulo, double *u0, double *ubar,
                  double *par, double per, int32 n, int32 icp1, int32 icp2, int32 icp3,
                  int32 icp4, int32 flag2, double *evr, double *evi, double tp);
void add_diagram(int32 ibr, int32 ntot, int32 itp, int32 lab, int32 nfpar, double a,
                 double *uhi, double *ulo, double *u0, double *ubar,
                 double *par, double per, int32 n, int32 icp1, int32 icp2, int32 icp3,
                 int32 icp4, int32 flag2, double *evr, double *evi);
void kill_diagrams(void);
void redraw_diagram(void);
void write_info_out(void);
void write_init_data_file(void);
void write_pts(void);
void post_auto(void);
void svg_auto(void);
void bound_diagram(double *xlo, double *xhi, double *ylo, double *yhi);
int32 save_diagram(FILE *fp, int32 n);
int32 load_diagram(FILE *fp, int32 node);

#endif
