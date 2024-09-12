#ifndef _dae_fun_h_
#define _dae_fun_h_
#include "integers.h"

int32 add_svar(char *name, char *rhs);
int32 add_svar_names(void);
int32 add_aeqn(char *rhs);
int32 compile_svars(void);
void reset_dae(void);
void set_init_guess(void);
void err_dae(void);
void init_dae_work(void);
void get_dae_fun(double *y, double *f);
void do_daes(void);
int32 solve_dae(void);
void get_new_guesses(void);

#endif
