#ifndef _extra_h_
#define _extra_h_

#define MAXW 50

void load_new_dll(void);
int my_fun(double *in, double *out, int nin, int nout, double *v, double *c);
void auto_load_dll(void);
void do_in_out(void);
void add_export_list(char *in, char *out);
void check_inout(void);
int get_export_count(char *s);
void do_export_list(void);
void parse_inout(char *l, int flag);
void get_import_values(int n, double *ydot, char *soname, char *sofun, int ivar,
                       double *wgt[MAXW], double *var, double *con);

#endif
