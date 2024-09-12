#ifndef _adj2_h_
#define _adj2_h_
#include "integers.h"

#include <stdio.h>

void init_trans(void);
void dump_transpose_info(FILE *fp, int32 f);
int32 do_transpose(void);
int32 create_transpose(void);
void alloc_h_stuff(void);
void data_back(void);
void adj_back(void);
void h_back(void);
void make_adj_com(int32 com);
void adjoint_parameters(void);
void new_h_fun(int32 silent);
void dump_h_stuff(FILE *fp, int32 f);
int32 make_h(float **orb, float **adj, int32 nt, int32 node, int32 silent);
void new_adjoint(void);
void test_test(void);
void compute_one_orbit(double *ic, double per);
int32 adjoint(float **orbit, float **adjnt, int32 nt, double dt, double eps,
              double minerr, int32 maxit, int32 node);
void eval_rhs(double **jac, int32 k1, int32 k2, double t, double *y, double *yp,
              int32 node);
int32 rk_interp(double **jac, int32 k1, int32 k2, double *y, double *work,
                int32 neq, double del, int32 nstep);
int32 step_eul(double **jac, int32 k, int32 k2, double *yold, double *work,
               int32 node, double dt);
void do_liapunov(void);
void alloc_liap(int32 n);
void do_this_liaprun(int32 i, double p);
void norm_vec(double *v, double *mu, int32 n);
int32 hrw_liapunov(double *liap, int32 batch, double eps);
void new_adjoint(void);

#endif
