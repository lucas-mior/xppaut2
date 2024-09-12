#ifndef _lunch_new_h_
#define _lunch_new_h_
#include "integers.h"

#include <stdio.h>
#include "form_ode.h"

void file_inf(void);
void ps_write_pars(FILE *fp);
void do_info(FILE *fp);
int32 read_lunch(FILE *fp);
void do_lunch(int32 f);
void dump_eqn(FILE *fp);
void io_numerics(int32 f, FILE *fp);
void io_parameter_file(char *fn, int32 flag);
void io_ic_file(char *fn, int32 flag);
void io_parameters(int32 f, FILE *fp);
void io_exprs(int32 f, FILE *fp);
void io_graph(int32 f, FILE *fp);
void io_int(int32 *i, FILE *fp, int32 f, char *ss);
void io_double(double *z, FILE *fp, int32 f, char *ss);
void io_float(float *z, FILE *fp, int32 f, char *ss);
void io_string(char *s, int32 len, FILE *fp, int32 f);

#endif
