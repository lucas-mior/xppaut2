#ifndef _tabular_h_
#define _tabular_h_
#include "integers.h"

void set_auto_eval_flags(int32 f);
void set_table_name(char *name, int32 index);
void view_table(int32 index);
void new_lookup_com(int32 i);
void new_lookup_ok(void);
double lookupxy(double x, int32 n, double *xv, double *yv);
double tab_interp(double xlo, double h, double x, double *y, int32 n, int32 i);
double lookup(double x, int32 index);
void init_table(void);
void redo_all_fun_tables(void);
int32 eval_fun_table(int32 n, double xlo, double xhi, char *formula, double *y);
int32 create_fun_table(int32 npts, double xlo, double xhi, char *formula,
                       int32 index);
int32 load_table(char *filename, int32 index);
int32 get_lookup_len(int32 i);

#endif
