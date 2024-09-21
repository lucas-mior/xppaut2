#ifndef parser_h
#define parser_h
#include "integers.h"

#include "volterra.h"
#include "xpplim.h"
#define MAX_SYMBS 10000
#define MAXARG 20
#define NEGATE 9
#define MINUS 4
#define LPAREN 0
#define RPAREN 1
#define COMMA 2
#define STARTTOK 10
#define ENDTOK 11
#define UFUN 24
#define MAX_UFUN 50
#define ENDEXP 999
#define ENDFUN 998
#define STARTDELAY 980
#define FUNTYPE 24
#define VARTYPE 2
#define CONTYPE 3
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
#define ENDISHIFT 985
#define DELSHFTSYM 65
#define SETSYM 72
#define ENDSET 981
/* #define EQSYM 68 */
#define INDXVAR 984
#define INDXCOM 922
#define INDX 68
#define STARTINDX 70
#define ENDINDX 69
#define MAX_TAB 50
/* #define MXLEN 32 */
#define MXLEN 10

/* number of standard symbols */

#define STDSYM 95

#define NEXT1VAR 26
#define NEXT2VAR 120

typedef struct {
    char name[MXLEN + 1];
    int32 len;
    int32 com;
    int32 arg;
    int32 pri;
    int32 id;
} SYMBOL;

SYMBOL my_symb[MAX_SYMBS] = {
    "(",
    1,
    999,
    0,
    1,
    0, /*  0   */
    ")",
    1,
    999,
    0,
    2,
    0,
    ",",
    1,
    999,
    0,
    3,
    0,
    "+",
    1,
    100,
    0,
    4,
    0,
    "-",
    1,
    101,
    0,
    4,
    0,
    "*",
    1,
    102,
    0,
    6,
    0,
    "/",
    1,
    103,
    0,
    6,
    0,
    "^",
    1,
    105,
    0,
    7,
    0,
    "**",
    2,
    105,
    0,
    7,
    0,
    "~",
    1,
    14,
    0,
    6,
    0,
    "START",
    5,
    -1,
    0,
    0,
    0, /* 10  */
    "END",
    3,
    999,
    0,
    -1,
    0,
    "ATAN2",
    5,
    104,
    2,
    10,
    0,
    "MAX",
    3,
    106,
    2,
    10,
    0,
    "MIN",
    3,
    107,
    2,
    10,
    0,
    "SIN",
    3,
    0,
    0,
    10,
    0,
    "COS",
    3,
    1,
    0,
    10,
    0,
    "TAN",
    3,
    2,
    0,
    10,
    0,
    "ASIN",
    4,
    3,
    0,
    10,
    0,
    "ACOS",
    4,
    4,
    0,
    10,
    0,
    "ATAN",
    4,
    5,
    0,
    10,
    0, /* 20  */
    "SINH",
    4,
    6,
    0,
    10,
    0,
    "TANH",
    4,
    7,
    0,
    10,
    0,
    "COSH",
    4,
    8,
    0,
    10,
    0,
    "ABS",
    3,
    9,
    0,
    10,
    0,
    "EXP",
    3,
    10,
    0,
    10,
    0,
    "LN",
    2,
    11,
    0,
    10,
    0,
    "LOG",
    3,
    11,
    0,
    10,
    0,
    "LOG10",
    5,
    12,
    0,
    10,
    0,
    "SQRT",
    4,
    13,
    0,
    10,
    0,
    "HEAV",
    4,
    16,
    0,
    10,
    0, /*  30 */
    "SIGN",
    4,
    17,
    0,
    10,
    0,
    "#$%1",
    4,
    800,
    0,
    10,
    0,
    "#$%2",
    4,
    801,
    0,
    10,
    0,
    "#$%3",
    4,
    802,
    0,
    10,
    0,
    "#$%4",
    4,
    803,
    0,
    10,
    0,
    "#$%5",
    4,
    804,
    0,
    10,
    0,
    "#$%6",
    4,
    805,
    0,
    10,
    0,
    "#$%7",
    4,
    806,
    0,
    10,
    0,
    "#$%8",
    4,
    807,
    0,
    10,
    0,
    "FLR",
    3,
    18,
    0,
    10,
    0, /*  40 */
    "MOD",
    3,
    108,
    2,
    10,
    0, /*  41 */
    "DELAY",
    5,
    ENDDELAY,
    2,
    10,
    0,
    /*  42 */ /*  Delay symbol */
    "RAN",
    3,
    19,
    1,
    10,
    0, /* 43 */
    "&",
    1,
    109,
    0,
    6,
    0, /* logical stuff  */
    "|",
    1,
    110,
    0,
    4,
    0,
    ">",
    1,
    111,
    0,
    7,
    0,
    "<",
    1,
    112,
    0,
    7,
    0,
    "==",
    2,
    113,
    0,
    7,
    0,
    ">=",
    2,
    114,
    0,
    7,
    0,
    "<=",
    2,
    115,
    0,
    7,
    0, /*50 */
    "IF",
    2,
    995,
    1,
    10,
    0,
    "THEN",
    4,
    994,
    1,
    10,
    0,
    "ELSE",
    4,
    993,
    1,
    10,
    0,
    "!=",
    2,
    116,
    0,
    7,
    0,
    "NOT",
    3,
    20,
    0,
    6,
    0,
    "NORMAL",
    6,
    117,
    2,
    10,
    0, /* returns normally dist number */
    "BESSELJ",
    7,
    118,
    2,
    10,
    0, /* Bessel J   */
    "BESSELY",
    7,
    119,
    2,
    10,
    0, /* Bessel Y */
    "NXXQQ",
    5,
    NUMSYM,
    0,
    10,
    0,
    "ERF",
    3,
    21,
    0,
    10,
    0, /* 60 */
    "ERFC",
    4,
    22,
    0,
    10,
    0,
    "SUM",
    3,
    SUMSYM,
    2,
    10,
    0,
    "OF",
    2,
    ENDSUM,
    0,
    10,
    0,
    "SHIFT",
    5,
    ENDSHIFT,
    2,
    10,
    0,
    "DEL_SHFT",
    8,
    ENDDELSHFT,
    3,
    10,
    0,
    "HOM_BCS",
    7,
    23,
    0,
    10,
    0,
    "ISHIFT",
    6,
    ENDISHIFT,
    2,
    10,
    0,
    "$",
    1,
    INDXCOM,
    0,
    10,
    0, /*68 */
    "]",
    1,
    ENDSHIFT,
    0,
    10,
    0,
    "[",
    1,
    ENDSHIFT,
    0,
    10,
    0, /*70 */
    "POISSON",
    7,
    24,
    0,
    10,
    0, /* 71 */
    "SET",
    3,
    ENDSET,
    3,
    10,
    0, /* 72 */
    "ARG1",
    4,
    800,
    0,
    10,
    0,
    "ARG2",
    4,
    801,
    0,
    10,
    0,
    "ARG3",
    4,
    802,
    0,
    10,
    0,
    "ARG4",
    4,
    803,
    0,
    10,
    0,
    "ARG5",
    4,
    804,
    0,
    10,
    0,
    "ARG6",
    4,
    805,
    0,
    10,
    0,
    "ARG7",
    4,
    806,
    0,
    10,
    0,
    "ARG8",
    4,
    807,
    0,
    10,
    0,
    "ARG9",
    4,
    808,
    0,
    10,
    0,
    "ARG10",
    5,
    809,
    0,
    10,
    0,
    "ARG11",
    5,
    810,
    0,
    10,
    0,
    "ARG12",
    5,
    811,
    0,
    10,
    0,
    "ARG13",
    5,
    812,
    0,
    10,
    0,
    "ARG14",
    5,
    813,
    0,
    10,
    0,
    "ARG15",
    5,
    814,
    0,
    10,
    0,
    "ARG16",
    5,
    815,
    0,
    10,
    0,
    "ARG17",
    5,
    816,
    0,
    10,
    0,
    "ARG18",
    5,
    817,
    0,
    10,
    0,
    "ARG19",
    5,
    818,
    0,
    10,
    0,
    "ARG20",
    5,
    819,
    0,
    10,
    0,
    "BESSELI",
    7,
    120,
    2,
    10,
    0, /* Bessel I  # 93 */
    "LGAMMA",
    6,
    25,
    1,
    10,
    0, /* log-gamma #94  */
};

#define IFUN1 0
#define IFUN2 1
#define IVAR 2
#define ICON 3
#define ILOOKUP 4
#define IUPTR 5
#define IKERN 6
#define ISHIFT 7
#define IUFUN 8
#define IDIV 9
#define ISUBTR 10
#define IMULT 11
#define IADD 12
#define INUMSYM 13
#define IENDFUN 14
#define IMYIF 15
#define IMYTHEN 16
#define IMYELSE 17
#define IENDDELAY 18
#define IENDSHIFT 19
#define ISUMSYM 20
#define IENDSUM 21
#define IENDEXPR 22
#define INETWORK 23
#define IENDDELSHFT 24
#define ISETEQUAL 25
#define IENDISHIFT 26
#define I_INDX 27
#define VECT_ROOT 500

int32 NSYM = STDSYM, NCON = 0, NVAR = 0, NFUN = 0;

/*     pointers to functions    */

void (*fun1[50])();
void (*fun2[50])();

/*****************************************************/

int32 stack_pointer, uptr;
double constants[MAX_PAR];
double variables[MAX_ODE1];
int32 *ufun[MAX_UFUN];
char *ufun_def[MAX_UFUN];
char ufun_names[MAX_UFUN][12];
int32 narg_fun[MAX_UFUN];
double stack[200], ustack[200];

KERNEL kernel[MAX_KER];
int32 NKernel;
int32 MaxPoints;
double *Memory[MAX_KER];
int32 NTable;

typedef struct {
    int32 narg;
    char args[MAXARG][11];
} UFUN_ARG;

UFUN_ARG ufun_arg[MAX_UFUN];

void zz_pmod(void);
void zz_atan2(void);
void zz_pow(void);
void zz_max(void);
void zz_min(void);
void zz_dand(void);
void zz_dor(void);
void zz_dnot(void);
void zz_dge(void);
void zz_dle(void);
void zz_dlt(void);
void zz_dgt(void);
void zz_deq(void);
void zz_dne(void);
void zz_normal(void);
void zz_bessel_j(void);
void zz_bessel_i(void);
void zz_bessel_y(void);
void zz_add(void);
void zz_subt(void);
void zz_div(void);
void zz_mult(void);
void zz_sin(void);
void zz_cos(void);
void zz_tan(void);
void zz_pois(void);
void zz_asin(void);
void zz_acos(void);
void zz_atan(void);
void zz_sinh(void);
void zz_cosh(void);
void zz_tanh(void);
void zz_fabs(void);
void zz_exp(void);
void zz_log(void);
void zz_log10(void);
void zz_sqrt(void);
void zz_neg(void);
void zz_recip(void);
void zz_heaviside(void);
void zz_signum(void);
void zz_floor(void);
void zz_erf(void);
void zz_erfc(void);
void zz_hom_bcs(void);
void zz_rndom(void);
void zz_lgamma(void);
int32 init_rpn(void);
int32 free_ufuns(void);
int32 duplicate_name(char *junk);
int32 add_constant(char *junk);
int32 get_type(int32 index);
int32 add_con(char *name, double value);
int32 add_kernel(char *name, double mu, char *expr);
int32 add_var(char *junk, double value);
int32 add_expr(char *expr, int32 *command, int32 *length);
int32 add_expr_no3(char *expr, int32 *command, int32 *length);
int32 add_vect_name(int32 index, char *name);
int32 add_net_name(int32 index, char *name);
int32 add_file_table(int32 index, char *file);
int32 add_table_name(int32 index, char *name);
int32 add_form_table(int32 index, int32 nn, double xlo, double xhi,
                     char *formula);
int32 set_old_arg_names(int32 narg);
int32 set_new_arg_names(int32 narg, char args[10][11]);
int32 fixup_endfun(int32 *u, int32 l, int32 narg);
int32 add_ufun_name(char *name, int32 index, int32 narg);
int32 add_ufun_new(int32 index, int32 narg, char *rhs, char args[20][11]);
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
int32 find_name(char *string, int32 *index);
int32 get_param_index(char *name);
int32 get_val(char *name, double *value);
int32 get_var_index(char *name);
int32 set_val(char *name, double value);
int32 set_ivar(int32 i, double value);
double get_ivar(int32 i);
int32 alg_to_rpn(int32 *toklist, int32 *command);
int32 pr_command(int32 *command);
int32 fpr_command(int32 *command);
int32 show_where(char *string, int32 index);
int32 function_sym(int32 token);
int32 unary_sym(int32 token);
int32 binary_sym(int32 token);
int32 pure_number(int32 token);
int32 gives_number(int32 token);
int32 check_syntax(int32 oldtoken, int32 newtoken);
int32 make_toks(char *dest, int32 *my_token);
int32 tokeninfo(int32 tok);
int32 do_num(char *source, char *num, double *value, int32 *ind);
int32 convert(char *source, char *dest);
int32 find_tok(char *source, int32 *index, int32 *tok);
double pmod(double x, double y);
int32 two_args(void);
double bessel_j(double x, double y);
double bessel_y(double x, double y);
double bessi(int32 n, double x);
double bessi0(double x);
double bessi1(double x);
char *com_name(int32 com);
double do_set_shift(double value, double shift, double variable);
double do_ishift(double shift, double variable);
double do_shift(double shift, double variable);
double do_delay_shift(double delay, double shift, double variable);
double do_delay(double delay, double i);
int32 one_arg(void);
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
int32 pass3(int32 *com1, int32 *com2, int32 *len);
double feval_rpn(int32 *comz);
int32 strupr(char *s);
int32 strlwr(char *s);

#endif
