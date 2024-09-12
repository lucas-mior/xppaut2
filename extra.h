#ifndef _extra_h_
#define _extra_h_

#define MAXW 50
#include "integers.h"

void load_new_dll(void);
int32 my_fun(double *in, double *out, int32 nin, int32 nout, double *v,
             double *c);
void auto_load_dll(void);
void do_in_out(void);
void add_export_list(char *in, char *out);
void check_inout(void);
int32 get_export_count(char *s);
void do_export_list(void);
void parse_inout(char *l, int32 flag);
void get_import_values(int32 n, double *ydot, char *soname, char *sofun,
                       int32 ivar, double *wgt[MAXW], double *var, double *con);

#endif
