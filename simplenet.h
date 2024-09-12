#ifndef _simplenet_h_
#define _simplenet_h_
#include "integers.h"

int32 get_vector_info();
double net_interp(double x, int32 i);
double network_value(double x, int32 i);
void init_net(double *v, int32 n);
int32 add_spec_fun(char *name, char *rhs);
void add_special_name(char *name, char *rhs);
int32 add_vectorizer(char *name, char *rhs);
void add_vectorizer_name(char *name, char *rhs);
int32 is_network(char *s);
void eval_all_nets(void);
void evaluate_network(int32 ind);
void update_all_ffts(void);
void update_fft(int32 ind);
void fft_conv(int32 it, int32 n, double *values, double *yy, double *fftr,
              double *ffti, double *dr, double *di);
int32 gilparse(char *s, int32 *ind, int32 *nn);
int32 g_namelist(char *s, char *root, int32 *flag, int32 *i1, int32 *i2);

#endif
