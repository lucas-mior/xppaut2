#ifndef _markov_h_
#define _markov_h_
#include "integers.h"

#include <stdio.h>

void add_wiener(int32 index);
void set_wieners(double dt, double *x, double t);
void add_markov(int32 nstate, char *name);
int32 build_markov(char **ma, char *name);
int32 old_build_markov(FILE *fptr, char *name);
void extract_expr(char *source, char *dest, int32 *i0);
void create_markov(int32 nstates, double *st, int32 type, char *name);
void add_markov_entry(int32 index, int32 j, int32 k, char *expr);
void compile_all_markov(void);
int32 compile_markov(int32 index, int32 j, int32 k);
void update_markov(double *x, double t, double dt);
double new_state(double old, int32 index, double dt);
void make_gill_nu(double *nu, int32 n, int32 m, double *v);
void one_gill_step(int32 meth, int32 nrxn, int32 *rxn, double *v);
void do_stochast_com(int32 i);
void mean_back(void);
void variance_back(void);
void compute_em(void);
void free_stoch(void);
void init_stoch(int32 len);
void append_stoch(int32 first, int32 length);
void do_stats(int32 ierr);
double gammln(double xx);
double poidev(double xm);
double ndrand48(void);
void nsrand48(int32 seed);
double ran1(long *idum);

#endif
