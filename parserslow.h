#ifndef _parserslow_h_
#define _parserslow_h_
#include "integers.h"

#include "volterra.h"
#include "xpplim.h"

#define FUN1TYPE 9
#define FUN2TYPE 1
#define VARTYPE 3 /* standard variable */
#define CONTYPE 2 /* standard parameter */
#define UFUNTYPE 24
#define SVARTYPE 4  /* shifted variable */
#define SCONTYPE 32 /* shifted constant  */
#define NETTYPE 6
#define TABTYPE 7
#define USTACKTYPE 8
#define KERTYPE 10
#define VECTYPE 13       /* for vectorized stuff */
#define EVECTYPE 14      /* treat vector like a function */
#define MAXTYPE 20000000 /* this is the maximum number of named stuff */

#define COM(a, b) ((a) * MAXTYPE + (b))

#define MAXARG 20
#define NEGATE 9
#define MINUS 4
#define LPAREN 0
#define RPAREN 1
#define COMMA 2
#define STARTTOK 10
#define ENDTOK 11

#define ENDEXP 999
#define ENDFUN 998
#define STARTDELAY 980
#define DELSYM 42
#define ENDDELAY 996
#define MYIF 995
#define MYELSE 993
#define MYTHEN 994
#define SUMSYM 990
#define ENDSUM 991
#define SHIFTSYM 64
#define ISHIFTSYM 67
#define ENDSHIFT 988
#define SUMINDEX 989
#define LASTTOK MAX_SYMBS - 2
#define NUMSYM 987
#define NUMTOK 59
#define CONV 2
#define FIRST_ARG 73
#define ENDDELSHFT 986
#define DELSHFTSYM 65
#define ENDISHIFT 985
#define SETSYM 72
#define ENDSET 981
#define INDX 68
#define INDXVAR 984

/*#define STDSYM 95
 */
#define STDSYM 96

#define INDXCOM 922
#define STARTINDX 70
#define ENDINDX 69

/* #define MXLEN 32 */
#define MXLEN 10

typedef struct {
    char name[MXLEN + 1];
    int32 len;
    int32 com;
    int32 arg;
    int32 pri;
} SYMBOL;

typedef struct {
    int32 narg;
    char args[MAXARG][14];
} UFUN_ARG;

#define VECT_ROOT 500

void init_rpn(void);
void free_ufuns(void);
int32 duplicate_name(char *junk);
int32 add_constant(char *junk);
int32 get_var_index(char *name);
int32 get_type(int32 index);
int32 add_con(char *name, double value);
int32 add_kernel(char *name, double mu, char *expr);
int32 add_var(char *junk, double value);
int32 add_expr(char *expr, int32 *command, int32 *length);
int32 add_net_name(int32 index, char *name);
int32 add_vector_name(int32 index, char *name);
int32 add_2d_table(char *name, char *file);
int32 add_file_table(int32 index, char *file);
int32 add_table_name(int32 index, char *name);
int32 add_form_table(int32 index, int32 nn, double xlo, double xhi,
                     char *formula);
void set_old_arg_names(int32 narg);
void set_new_arg_names(int32 narg, char args[10][14]);
int32 add_ufun_name(char *name, int32 index, int32 narg);
void fixup_endfun(int32 *u, int32 l, int32 narg);
int32 add_ufun_new(int32 index, int32 narg, char *rhs, char args[20][14]);
int32 add_ufun(char *junk, char *expr, int32 narg);
int32 check_num(int32 *tok, double value);
int32 is_ufun(int32 x);
int32 is_ucon(int32 x);
int32 is_uvar(int32 x);
int32 isvar(int32 y);
int32 iscnst(int32 y);
int32 isker(int32 y);
int32 is_kernel(int32 x);
int32 is_lookup(int32 x);
int32 find_lookup(char *name);
void find_name(char *string, int32 *index);
int32 get_param_index(char *name);
int32 get_val(char *name, double *value);
int32 set_val(char *name, double value);
void set_ivar(int32 i, double value);
double get_ivar(int32 i);
int32 alg_to_rpn(int32 *toklist, int32 *command);
void pr_command(int32 *command);
void show_where(char *string, int32 index);
int32 function_sym(int32 token);
int32 unary_sym(int32 token);
int32 binary_sym(int32 token);
int32 pure_number(int32 token);
int32 gives_number(int32 token);
int32 check_syntax(int32 oldtoken, int32 newtoken);
int32 make_toks(char *dest, int32 *my_token);
void tokeninfo(int32 tok);
int32 do_num(char *source, char *num, double *value, int32 *ind);
void convert(char *source, char *dest);
void find_tok(char *source, int32 *index, int32 *tok);
double pmod(double x, double y);
void two_args(void);
double bessel_j(double x, double y);
double bessel_y(double x, double y);
double bessi(double nn, double x);
double bessi0(double x);
double bessi1(double x);
double bessis(double nn, double x);
double bessis0(double x);
double bessis1(double x);
char *com_name(int32 com);
double do_shift(double shift, double variable);
double do_ishift(double shift, double variable);
double do_delay_shift(double delay, double shift, double variable);
double do_delay(double delay, double i);
void one_arg(void);
double normal(double mean, double std);
double max(double x, double y);
double min(double x, double y);

double neg(double z);
double recip(double z);
double heaviside(double z);
double rndom(double z);
double signum(double z);
double dnot(double x);
double dand(double x, double y);
double dor(double x, double y);
double dge(double x, double y);
double dle(double x, double y);
double deq(double x, double y);
double dne(double x, double y);
double dgt(double x, double y);
double dlt(double x, double y);
double evaluate(int32 *equat);
double eval_rpn(int32 *equat);

/*  STRING STUFF  */
#ifndef STRUPR
void strupr(char *s);
void strlwr(char *s);
#endif

/*****************************************************/

#endif
