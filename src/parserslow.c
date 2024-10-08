#include "parserslow.h"
#include "integers.h"
#include "xmalloc.h"
#include <stdbool.h>

#include <time.h>
#include "functions.h"

#include <stdlib.h>

#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#include <math.h>
/* #include <xmalloc.h> */
#include <stdio.h>
#include <string.h>

#include "xpplim.h"

#define MAXEXPLEN 1024
#define DOUB_EPS 2.23E-15
#define POP stack[--stack_pointer]
static double zippy;
#define PUSH(a)                                                                                    \
    do {                                                                                           \
        zippy = (a);                                                                               \
        stack[stack_pointer++] = zippy;                                                            \
    } while (0)

/* #define COM(a) my_symb[toklist[(a)]].com */
int32 errout;
int32 NDELAYS = 0;
static double hom_bcs(int32);
static double box_muller;
static int32 box_muller_flag = 0;
int32 rand_seed = 12345678;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

static double current_index = 0;
static int32 SumIndex = 1;

/* FIXXX */
static int32 stack_pointer;
static int32 uptr;
double constants[MAX_PAR];
double variables[MAX_ODE1];
int32 *ufun[MAX_UFUN];
char *ufun_def[MAX_UFUN];
char ufun_names[MAX_UFUN][12];
int32 narg_fun[MAX_UFUN];
static double stack[200];
static double ustack[200];

Kernel kernel[MAX_KER];
int32 nkernel;
int32 max_points;
int32 NTable;

UFUN_ARG ufun_arg[MAX_UFUN];

static struct Symbol {
    char name[MXLEN + 1];
    int32 len;
    int32 com;
    int32 arg;
    int32 pri;
} my_symb[MAX_SYMBS] = {
    {"(", 1, 999, 0, 1},  //  0
    {")", 1, 999, 0, 2},
    {",", 1, 999, 0, 3},
    {"+", 1, COM(FUN2TYPE, 0), 0, 4},
    {"-", 1, COM(FUN2TYPE, 1), 0, 4},
    {"*", 1, COM(FUN2TYPE, 2), 0, 6},
    {"/", 1, COM(FUN2TYPE, 3), 0, 6},
    {"^", 1, COM(FUN2TYPE, 5), 0, 7},
    {"**", 2, COM(FUN2TYPE, 5), 0, 7},
    {"~", 1, COM(FUN1TYPE, 14), 0, 6},
    {"START", 5, -1, 0, 0},  // 10
    {"END", 3, 999, 0, -1},
    {"ATAN2", 5, COM(FUN2TYPE, 4), 2, 10},
    {"MAX", 3, COM(FUN2TYPE, 6), 2, 10},
    {"MIN", 3, COM(FUN2TYPE, 7), 2, 10},
    {"SIN", 3, COM(FUN1TYPE, 0), 0, 10},
    {"COS", 3, COM(FUN1TYPE, 1), 0, 10},
    {"TAN", 3, COM(FUN1TYPE, 2), 0, 10},
    {"ASIN", 4, COM(FUN1TYPE, 3), 0, 10},
    {"ACOS", 4, COM(FUN1TYPE, 4), 0, 10},
    {"ATAN", 4, COM(FUN1TYPE, 5), 0, 10},  // 20
    {"SINH", 4, COM(FUN1TYPE, 6), 0, 10},
    {"TANH", 4, COM(FUN1TYPE, 7), 0, 10},
    {"COSH", 4, COM(FUN1TYPE, 8), 0, 10},
    {"ABS", 3, COM(FUN1TYPE, 9), 0, 10},
    {"EXP", 3, COM(FUN1TYPE, 10), 0, 10},
    {"LN", 2, COM(FUN1TYPE, 11), 0, 10},
    {"LOG", 3, COM(FUN1TYPE, 11), 0, 10},
    {"LOG10", 5, COM(FUN1TYPE, 12), 0, 10},
    {"SQRT", 4, COM(FUN1TYPE, 13), 0, 10},
    {"HEAV", 4, COM(FUN1TYPE, 16), 0, 10},  //  30
    {"SIGN", 4, COM(FUN1TYPE, 17), 0, 10},
    {"#$%1", 4, COM(USTACKTYPE, 0), 0, 10},
    {"#$%2", 4, COM(USTACKTYPE, 1), 0, 10},
    {"#$%3", 4, COM(USTACKTYPE, 2), 0, 10},
    {"#$%4", 4, COM(USTACKTYPE, 3), 0, 10},
    {"#$%5", 4, COM(USTACKTYPE, 4), 0, 10},
    {"#$%6", 4, COM(USTACKTYPE, 5), 0, 10},
    {"#$%7", 4, COM(USTACKTYPE, 6), 0, 10},
    {"#$%8", 4, COM(USTACKTYPE, 7), 0, 10},
    {"FLR", 3, COM(FUN1TYPE, 18), 0, 10},  //  40
    {"MOD", 3, COM(FUN2TYPE, 8), 2, 10},   //  41
    {"delay2", 5, ENDDELAY, 2, 10},
    /*  42 */                              //  delay2 symbol
    {"RAN", 3, COM(FUN1TYPE, 19), 1, 10},  // 43
    {"&", 1, COM(FUN2TYPE, 9), 0, 6},      // logical stuff
    {"|", 1, COM(FUN2TYPE, 10), 0, 4},
    {">", 1, COM(FUN2TYPE, 11), 0, 7},
    {"<", 1, COM(FUN2TYPE, 12), 0, 7},
    {"==", 2, COM(FUN2TYPE, 13), 0, 7},
    {">=", 2, COM(FUN2TYPE, 14), 0, 7},
    {"<=", 2, COM(FUN2TYPE, 15), 0, 7},  // 50
    {"IF", 2, 995, 1, 10},
    {"THEN", 4, 994, 1, 10},
    {"ELSE", 4, 993, 1, 10},
    {"!=", 2, COM(FUN2TYPE, 16), 0, 7},
    {"NOT", 3, COM(FUN1TYPE, 20), 0, 6},
    {"NORMAL", 6, COM(FUN2TYPE, 17), 2, 10},   // returns normally dist number
    {"BESSELJ", 7, COM(FUN2TYPE, 18), 2, 10},  // Bessel J
    {"BESSELY", 7, COM(FUN2TYPE, 19), 2, 10},  // Bessel Y
    {"NXXQQ", 5, NUMSYM, 0, 10},
    {"ERF", 3, COM(FUN1TYPE, 21), 0, 10},  // 60
    {"ERFC", 4, COM(FUN1TYPE, 22), 0, 10},
    {"SUM", 3, SUMSYM, 2, 10},
    {"OF", 2, ENDSUM, 0, 10},
    {"SHIFT", 5, ENDSHIFT, 2, 10},
    {"DEL_SHFT", 8, ENDDELSHFT, 3, 10},  // 65
    {"HOM_BCS", 7, COM(FUN1TYPE, 23), 0, 10},
    {"ISHIFT", 6, ENDISHIFT, 2, 10},  // 67
    {"@", 1, INDXCOM, 0, 10},         // 68
    {"]", 1, ENDSHIFT, 0, 10},
    {"[", 1, ENDSHIFT, 0, 10},                 // 70
    {"POISSON", 7, COM(FUN1TYPE, 24), 0, 10},  // 71
    {"SET", 3, ENDSET, 3, 10},                 // 72
    {"ARG1", 4, COM(USTACKTYPE, 0), 0, 10},    //  FIXXX ????
    {"ARG2", 4, COM(USTACKTYPE, 1), 0, 10},
    {"ARG3", 4, COM(USTACKTYPE, 2), 0, 10},
    {"ARG4", 4, COM(USTACKTYPE, 3), 0, 10},
    {"ARG5", 4, COM(USTACKTYPE, 4), 0, 10},
    {"ARG6", 4, COM(USTACKTYPE, 5), 0, 10},
    {"ARG7", 4, COM(USTACKTYPE, 6), 0, 10},
    {"ARG8", 4, COM(USTACKTYPE, 7), 0, 10},
    {"ARG9", 4, COM(USTACKTYPE, 8), 0, 10},
    {"ARG10", 5, COM(USTACKTYPE, 9), 0, 10},
    {"ARG11", 5, COM(USTACKTYPE, 10), 0, 10},
    {"ARG12", 5, COM(USTACKTYPE, 11), 0, 10},
    {"ARG13", 5, COM(USTACKTYPE, 12), 0, 10},
    {"ARG14", 5, COM(USTACKTYPE, 13), 0, 10},
    {"ARG15", 5, COM(USTACKTYPE, 14), 0, 10},
    {"ARG16", 5, COM(USTACKTYPE, 15), 0, 10},
    {"ARG17", 5, COM(USTACKTYPE, 16), 0, 10},
    {"ARG18", 5, COM(USTACKTYPE, 17), 0, 10},
    {"ARG19", 5, COM(USTACKTYPE, 18), 0, 10},
    {"ARG20", 5, COM(USTACKTYPE, 19), 0, 10},
    {"BESSELI", 7, COM(FUN2TYPE, 20), 2, 10},   // Bessel I  # 93
    {"LGAMMA", 6, COM(FUN1TYPE, 25), 1, 10},    // Log Gamma  #94
    {"BESSELIS", 8, COM(FUN2TYPE, 21), 2, 10},  // Bessel I Scaled  # 95
};

int32 NCON = 0;
int32 NVAR = 0;
int32 nfun = 0;
int32 NSYM = STDSYM;

typedef double (*FunctionInt32)(int32);
typedef double (*FunctionInt32Int32)(int32, int32);
typedef double (*FunctionDouble)(double);
typedef double (*FunctionDoubleDouble)(double, double);

static void *fun1[50];
static void *fun2[50];

/*************************
  RPN COMPILER           *
**************************/

/*****************************
 *      PARSER.C              *
 *
 *
 *     parses any algebraic expression
 *     and converts to an int64 array
 *     to be interpreted by the rpe_val
 *     function.
 *
 *     the main data structure is a contiguous
 *     list of symbols with their priorities
 *     and their symbol value
 *
 *     on the first pass, the expression is converted to
 *     a list of integers without any checking except for
 *     valid symbols and for numbers
 *
 *     next this list of integers is converted to an RPN expression
 *     for evaluation.
 *
 *
 *  6/95  stuff added to add names to namelist without compilation
 *************************************************************/

void
init_rpn(void) {
    errout = 1;
    NCON = 0;
    nfun = 0;
    NVAR = 0;
    nkernel = 0;

    max_points = 4000;
    NSYM = STDSYM;

    fun2[4] = (void *)atan2;
    fun2[5] = (void *)pow;
    fun2[6] = (void *)max;
    fun2[7] = (void *)min;
    //  fun2[8]= (void*)mod;
    fun2[8] = (void *)pmod;  // This always gives an answer in [0,y) for mod(x,y)
    fun2[9] = (void *)dand;
    fun2[10] = (void *)dor;
    fun2[11] = (void *)dgt;
    fun2[12] = (void *)dlt;
    fun2[13] = (void *)deq;
    fun2[14] = (void *)dge;
    fun2[15] = (void *)dle;
    fun2[16] = (void *)dne;
    fun2[17] = (void *)normal;
    fun2[18] = (void *)bessel_j;
    fun2[19] = (void *)bessel_y;
    fun2[20] = (void *)bessi;
    fun2[21] = (void *)bessis;

    fun1[0] = (void *)sin;
    fun1[1] = (void *)cos;
    fun1[2] = (void *)tan;
    fun1[3] = (void *)asin;
    fun1[4] = (void *)acos;
    fun1[5] = (void *)atan;
    fun1[6] = (void *)sinh;
    fun1[7] = (void *)tanh;
    fun1[8] = (void *)cosh;
    fun1[9] = (void *)fabs;
    fun1[10] = (void *)exp;
    fun1[11] = (void *)log;
    fun1[12] = (void *)log10;
    fun1[13] = (void *)sqrt;
    fun1[14] = (void *)neg;
    fun1[15] = (void *)recip;
    fun1[16] = (void *)heaviside;
    fun1[17] = (void *)signum;
    fun1[18] = (void *)floor;
    fun1[19] = (void *)rndom;
    fun1[20] = (void *)dnot;
    fun1[21] = (void *)erf;
    fun1[22] = (void *)erfc;
    fun1[23] = (void *)hom_bcs;
    fun1[24] = (void *)markov_poidev;
    fun1[25] = (void *)lgamma;

    parserslow_add_con("PI", M_PI);

    parserslow_add_con("I'", 0.0);
    /*   This is going to be for interacting with the
         animator */
    SumIndex = NCON - 1;
    parserslow_add_con("mouse_x", 0.0);
    parserslow_add_con("mouse_y", 0.0);
    parserslow_add_con("mouse_vx", 0.0);
    parserslow_add_con("mouse_vy", 0.0);

    // end animator stuff

    tabular_init_table();
    if (newseed == 1) {
        rand_seed = (int32)time(0);
    }
    markov_nsrand48(rand_seed);
    return;
}

int32
duplicate_name(char *junk) {
    int32 i;
    find_name(junk, &i);
    if (i >= 0) {
        if (errout) {
            printf("%s is a duplicate name\n", junk);
        }
        return 1;
    }
    return 0;
}

/*  ADD_CONSTANT   */

int32
parserslow_add_constant(char *junk) {
    int32 len;
    char string[100];
    if (duplicate_name(junk) == 1) {
        return 1;
    }

    if (NCON >= MAX_PAR) {
        if (errout) {
            printf("too many constants !!\n");
        }
        return 1;
    }
    convert(junk, string);
    len = (int32)strlen(string);
    if (len < 1) {
        ggets_plintf("Empty parameter - remove spaces\n");
        return 1;
    }
    if (len > MXLEN) {
        len = MXLEN;
    }
    strncpy(my_symb[NSYM].name, string, (usize)len);
    my_symb[NSYM].name[len] = '\0';
    my_symb[NSYM].len = len;
    my_symb[NSYM].pri = 10;
    my_symb[NSYM].arg = 0;
    my_symb[NSYM].com = COM(CONTYPE, NCON - 1);
    NSYM++;
    return 0;
}

int32
get_var_index(char *name) {
    int32 type;
    int32 com;
    find_name(name, &type);
    if (type < 0) {
        return -1;
    }
    com = my_symb[type].com;
    if (is_uvar(com)) {
        return com % MAXTYPE;
    }
    return -1;
}

/* GET_TYPE   */

int32
get_type(int32 index) {
    return my_symb[index].com;
}

/*   ADD_CON      */

int32
parserslow_add_con(char *name, double value) {
    if (NCON >= MAX_PAR) {
        if (errout) {
            printf("too many constants !!\n");
        }
        return 1;
    }
    constants[NCON] = value;
    NCON++;
    return parserslow_add_constant(name);
}

int32
parserslow_add_kernel(char *name, double mu, char *expr) {
    char string[100];

    int32 len;
    int32 in = -1;
    if (duplicate_name(name) == 1) {
        return 1;
    }
    if (nkernel == MAX_KER) {
        ggets_plintf("Too many kernels..\n");
        return 1;
    }
    if (mu < 0 || mu >= 1.0) {
        ggets_plintf(" mu must lie in [0,1.0) \n");
        return 1;
    }
    convert(name, string);
    len = (int32)strlen(string);
    if (len > MXLEN) {
        len = MXLEN;
    }
    strncpy(my_symb[NSYM].name, string, (usize)len);
    my_symb[NSYM].name[len] = '\0';
    my_symb[NSYM].len = len;
    my_symb[NSYM].pri = 10;
    my_symb[NSYM].arg = 0;
    my_symb[NSYM].com = COM(KERTYPE, nkernel);
    kernel[nkernel].k_n1 = 0.0;
    kernel[nkernel].mu = mu;
    kernel[nkernel].k_n = 0.0;
    kernel[nkernel].k_n1 = 0.0;
    kernel[nkernel].flag = 0;
    for (int32 i = 0; i < (int32)strlen(expr); i++) {
        if (expr[i] == '#') {
            in = i;
        }
    }
    if (in == 0 || in == ((int32)strlen(expr) - 1)) {
        ggets_plintf("Illegal use of convolution...\n");
        return 1;
    }
    if (in > 0) {
        kernel[nkernel].flag = CONV;
        kernel[nkernel].expr = xmalloc(strlen(expr) + 2 - (usize)in);
        kernel[nkernel].kerexpr = xmalloc((usize)in + 1);
        for (int32 i = 0; i < in; i++) {
            kernel[nkernel].kerexpr[i] = expr[i];
        }
        kernel[nkernel].kerexpr[in] = 0;
        for (int32 i = in + 1; i < (int32)strlen(expr); i++) {
            kernel[nkernel].expr[i - in - 1] = expr[i];
        }
        kernel[nkernel].expr[strlen(expr) - (usize)in - 1] = 0;
        ggets_plintf("Convolving %s with %s\n", kernel[nkernel].kerexpr, kernel[nkernel].expr);
    } else {
        kernel[nkernel].expr = xmalloc(strlen(expr) + 2);
        strcpy(kernel[nkernel].expr, expr);
    }
    NSYM++;
    nkernel++;
    return 0;
}

/*  ADD_VAR          */

int32
parserslow_add_var(char *junk, double value) {
    char string[100];
    int32 len;
    if (duplicate_name(junk) == 1) {
        return 1;
    }
    if (NVAR >= MAX_ODE1) {
        if (errout) {
            printf("too many variables !!\n");
        }
        return 1;
    }
    convert(junk, string);
    len = (int32)strlen(string);
    if (len > MXLEN) {
        len = MXLEN;
    }
    strncpy(my_symb[NSYM].name, string, (usize)len);
    my_symb[NSYM].name[len] = '\0';
    my_symb[NSYM].len = len;
    my_symb[NSYM].pri = 10;
    my_symb[NSYM].arg = 0;
    my_symb[NSYM].com = COM(VARTYPE, NVAR);
    NSYM++;
    variables[NVAR] = value;
    NVAR++;
    return 0;
}

/* ADD_EXPR   */

int32
parserslow_add_expr(char *expr, int32 *command, int32 *length) {
    char dest[1024];
    int32 my_token[1024];
    int32 err;
    int32 i;
    convert(expr, dest);
    err = make_toks(dest, my_token);

    if (err != 0) {
        return 1;
    }
    err = parserslow_alg_to_rpn(my_token, command);
    if (err != 0) {
        return 1;
    }
    i = 0;
    while (command[i] != ENDEXP) {
        i++;
    }
    *length = i + 1;
    return 0;
}

int32
parserslow_add_vector_name(int32 index, char *name) {
    char string[50];
    int32 len = (int32)strlen(name);
    ggets_plintf(" Adding vectorizer %s %d \n", name, index);
    if (duplicate_name(name) == 1) {
        return 1;
    }
    convert(name, string);
    printf(" 1\n");
    if (len > MXLEN) {
        len = MXLEN;
    }
    strncpy(my_symb[NSYM].name, string, (usize)len);
    my_symb[NSYM].name[len] = '\0';
    my_symb[NSYM].len = len;
    my_symb[NSYM].pri = 10;
    my_symb[NSYM].arg = 1;
    my_symb[NSYM].com = COM(VECTYPE, index);

    NSYM++;
    return 0;
}

int32
parserslow_add_net_name(int32 index, char *name) {
    char string[50];
    int32 len = (int32)strlen(name);
    ggets_plintf(" Adding net %s %d \n", name, index);
    if (duplicate_name(name) == 1) {
        return 1;
    }
    convert(name, string);
    if (len > MXLEN) {
        len = MXLEN;
    }
    strncpy(my_symb[NSYM].name, string, (usize)len);
    my_symb[NSYM].name[len] = '\0';
    my_symb[NSYM].len = len;
    my_symb[NSYM].pri = 10;
    my_symb[NSYM].arg = 1;
    my_symb[NSYM].com = COM(NETTYPE, index);
    NSYM++;
    return 0;
}

/* ADD LOOKUP TABLE   */

int32
parserslow_add_file_table(int32 index, char *file) {
    char file2[1000];
    int32 i2 = 0;
    int32 i1 = 0;
    int32 n;
    char ch;
    n = (int32)strlen(file);
    for (i1 = 0; i1 < n; i1++) {
        ch = file[i1];
        if (((int32)ch > 31) && ((int32)ch < 127)) {
            file2[i2] = ch;
            i2++;
        }
    }
    file2[i2] = 0;
    if (tabular_load_table(file2, index) == 0) {
        if (errout) {
            printf("Problem with creating table !!\n");
        }
        return 1;
    }

    return 0;
}

int32
parserslow_add_table_name(int32 index, char *name) {
    char string[50];
    int32 len = (int32)strlen(name);
    if (duplicate_name(name) == 1) {
        return 1;
    }
    convert(name, string);
    if (len > MXLEN) {
        len = MXLEN;
    }
    strncpy(my_symb[NSYM].name, string, (usize)len);
    my_symb[NSYM].name[len] = '\0';
    my_symb[NSYM].len = len;
    my_symb[NSYM].pri = 10;
    my_symb[NSYM].arg = 1;
    my_symb[NSYM].com = COM(TABTYPE, index);
    tabular_set_table_name(name, index);
    NSYM++;
    return 0;
}
/* ADD LOOKUP TABLE   */

int32
parserslow_add_form_table(int32 index, int32 nn, double xlo, double xhi, char *formula) {
    if (tabular_create_fun(nn, xlo, xhi, formula, index) == 0) {
        if (errout) {
            printf("Problem with creating table !!\n");
        }
        return 1;
    }
    return 0;
}

void
set_old_arg_names(int32 narg) {
    for (int32 i = 0; i < narg; i++) {
        snprintf(my_symb[FIRST_ARG + i].name, sizeof(my_symb[FIRST_ARG + i].name), "ARG%d", i + 1);
        my_symb[FIRST_ARG + i].len = 4;
    }
    return;
}

void
set_new_arg_names(int32 narg, char args[10][14]) {
    for (int32 i = 0; i < narg; i++) {
        strcpy(my_symb[FIRST_ARG + i].name, args[i]);
        my_symb[FIRST_ARG + i].len = (int32)strlen(args[i]);
    }
    return;
}

/* NEW ADD_FUN for new form_ode code  */

int32
parserslow_add_ufun_name(char *name, int32 index, int32 narg) {
    char string[50];
    int32 len = (int32)strlen(name);
    if (duplicate_name(name) == 1) {
        return 1;
    }
    if (index >= MAX_UFUN) {
        if (errout) {
            printf("too many functions !!\n");
        }
        return 1;
    }
    ggets_plintf(" Added user fun %s \n", name);
    convert(name, string);
    if (len > MXLEN) {
        len = MXLEN;
    }
    strncpy(my_symb[NSYM].name, string, (usize)len);
    my_symb[NSYM].name[len] = '\0';
    my_symb[NSYM].len = len;
    my_symb[NSYM].pri = 10;
    my_symb[NSYM].arg = narg;
    my_symb[NSYM].com = COM(UFUNTYPE, index);
    NSYM++;
    strcpy(ufun_names[index], name);
    return 0;
}

void
fixup_endfun(int32 *u, int32 l, int32 narg) {
    u[l - 1] = ENDFUN;
    u[l] = narg;
    u[l + 1] = ENDEXP;
    return;
}

int32
parserslow_add_ufun_new(int32 index, int32 narg, char *rhs, char args[MAXARG][14]) {
    int32 l;
    int32 end;
    if (narg > MAXARG) {
        ggets_plintf("Maximal arguments exceeded \n");
        return 1;
    }
    if ((ufun[index] = xmalloc(1024)) == NULL) {
        if (errout) {
            printf("not enough memory!!\n");
        }
        return 1;
    }
    if ((ufun_def[index] = xmalloc(MAXEXPLEN)) == NULL) {
        if (errout) {
            printf("not enough memory!!\n");
        }
        return 1;
    }
    ufun_arg[index].narg = narg;
    for (int32 i = 0; i < narg; i++) {
        strcpy(ufun_arg[index].args[i], args[i]);
    }
    set_new_arg_names(narg, args);
    if (parserslow_add_expr(rhs, ufun[index], &end) == 0) {
        ufun[index][end - 1] = ENDFUN;
        ufun[index][end] = narg;
        ufun[index][end + 1] = ENDEXP;
        strcpy(ufun_def[index], rhs);
        l = (int32)strlen(ufun_def[index]);
        ufun_def[index][l] = 0;
        narg_fun[index] = narg;
        set_old_arg_names(narg);
        return 0;
    }

    set_old_arg_names(narg);
    if (errout) {
        printf(" ERROR IN FUNCTION DEFINITION\n");
    }
    return 1;
}

/* ADD_UFUN   */

int32
parserslow_add_ufun(char *junk, char *expr, int32 narg) {
    char string[50];
    int32 l;
    int32 end;
    int32 len = (int32)strlen(junk);

    if (duplicate_name(junk) == 1) {
        return 1;
    }
    if (nfun >= MAX_UFUN) {
        if (errout) {
            printf("too many functions !!\n");
        }
        return 1;
    }
    if ((ufun[nfun] = xmalloc(1024)) == NULL) {
        if (errout) {
            printf("not enough memory!!\n");
        }
        return 1;
    }
    if ((ufun_def[nfun] = xmalloc(MAXEXPLEN)) == NULL) {
        if (errout) {
            printf("not enough memory!!\n");
        }
        return 1;
    }

    convert(junk, string);
    if (parserslow_add_expr(expr, ufun[nfun], &end) == 0) {
        if (len > MXLEN) {
            len = MXLEN;
        }
        strncpy(my_symb[NSYM].name, string, (usize)len);
        my_symb[NSYM].name[len] = '\0';
        my_symb[NSYM].len = len;
        my_symb[NSYM].pri = 10;
        my_symb[NSYM].arg = narg;
        my_symb[NSYM].com = COM(UFUNTYPE, nfun);
        NSYM++;
        ufun[nfun][end - 1] = ENDFUN;
        ufun[nfun][end] = narg;
        ufun[nfun][end + 1] = ENDEXP;
        strcpy(ufun_def[nfun], expr);
        l = (int32)strlen(ufun_def[nfun]);
        ufun_def[nfun][l - 1] = 0;
        strcpy(ufun_names[nfun], junk);
        narg_fun[nfun] = narg;
        for (int32 i = 0; i < narg; i++) {
            snprintf(ufun_arg[nfun].args[i], sizeof(ufun_arg[nfun].args[i]), "ARG%d", i + 1);
        }
        nfun++;
        return 0;
    }
    if (errout) {
        printf(" ERROR IN FUNCTION DEFINITION\n");
    }
    return 1;
}

int32
check_num(int32 *tok, double value) {
    int32 bob;
    int32 in;
    for (int32 i = 0; i < NSYM; i++) {
        if (strncmp(my_symb[i].name, "NUM##", 5) == 0) {
            bob = my_symb[i].com;
            in = bob % MAXTYPE;
            if (constants[in] == value) {
                *tok = i;
                return 1;
            }
        }
    }
    return 0;
}

/* is_ufun         */

int32
is_ufun(int32 x) {
    if ((x / MAXTYPE) == UFUNTYPE) {
        return 1;
    } else {
        return 0;
    }
}

/* IS_UCON        */

int32
is_ucon(int32 x) {
    if (x / MAXTYPE == CONTYPE) {
        return 1;
    } else {
        return 0;
    }
}

/* IS_UVAR       */

int32
is_uvar(int32 x) {
    if (x / MAXTYPE == VARTYPE) {
        return 1;
    } else {
        return 0;
    }
}

int32
isvar(int32 y) {
    return y == VARTYPE;
}

int32
iscnst(int32 y) {
    return y == CONTYPE;
}

int32
isker(int32 y) {
    return y == KERTYPE;
}

int32
is_kernel(int32 x) {
    if ((x / MAXTYPE) == KERTYPE) {
        return 1;
    } else {
        return 0;
    }
}

int32
is_lookup(int32 x) {
    if ((x / MAXTYPE) == TABTYPE) {
        return 1;
    } else {
        return 0;
    }
}

int32
find_lookup(char *name) {
    int32 index;
    int32 com;
    find_name(name, &index);
    if (index == -1) {
        return -1;
    }
    com = my_symb[index].com;
    if (is_lookup(com)) {
        return com % MAXTYPE;
    }
    return -1;
}

/* FIND_NAME    */

void
find_name(char *string, int32 *index) {
    char junk[100];
    int32 i;
    int32 len;
    convert(string, junk);
    len = (int32)strlen(junk);
    for (i = 0; i < NSYM; i++) {
        if (len == my_symb[i].len) {
            if (strncmp(my_symb[i].name, junk, (usize)len) == 0) {
                break;
            }
        }
    }
    if (i < NSYM) {
        *index = i;
    } else {
        *index = -1;
    }
    return;
}

int32
get_param_index(char *name) {
    int32 type;
    int32 com;
    find_name(name, &type);
    if (type < 0) {
        return -1;
    }
    com = my_symb[type].com;
    if (is_ucon(com)) {
        return com % MAXTYPE;
    }
    return -1;
}

/* GET_VAL   */

int32
get_val(char *name, double *value) {
    int32 type;
    int32 com;
    *value = 0.0;
    find_name(name, &type);
    if (type < 0) {
        return 0;
    }
    com = my_symb[type].com;
    if (is_ucon(com)) {
        *value = constants[com % MAXTYPE];
        return 1;
    }
    if (is_uvar(com)) {
        *value = variables[com % MAXTYPE];
        return 1;
    }
    return 0;
}

/* SET_VAL         */

int32
set_val(char *name, double value) {
    int32 type;
    int32 com;
    find_name(name, &type);
    if (type < 0) {
        return 0;
    }
    com = my_symb[type].com;
    if (is_ucon(com)) {
        constants[com % MAXTYPE] = value;

        return 1;
    }
    if (is_uvar(com)) {
        variables[com % MAXTYPE] = value;
        return 1;
    }
    return 0;
}

void
set_ivar(int32 i, double value) {
    SETVAR(i, value);
    return;
}

double
get_ivar(int32 i) {
    return GETVAR(i);
}

int32
parserslow_alg_to_rpn(int32 *toklist, int32 *command) {
    int32 tokstak[500];
    int32 comptr = 0;
    int32 tokptr = 0;
    int32 lstptr = 0;
    int32 temp;
    int32 ncomma = 0;
    int32 loopstk[100];
    int32 lptr = 0;
    int32 nif = 0;
    int32 nthen = 0;
    int32 nelse = 0;
    int32 newtok;
    int32 oldtok;
    int32 my_com;
    int32 my_arg;
    int32 jmp;

    tokstak[0] = STARTTOK;
    tokptr = 1;
    oldtok = STARTTOK;
    while (true) {
    getnew:
        newtok = toklist[lstptr++];
        //        check for delay2 symbol
        if (newtok == DELSYM) {
            temp = my_symb[toklist[lstptr + 1]].com;
            /* !! */ if (is_uvar(temp)) {
                /* ram -- is this right? not sure I understand what was
                 * happening here
                 */
                my_symb[LASTTOK].com = COM(SVARTYPE, temp % MAXTYPE);  // create a temporary sybol
                NDELAYS++;
                toklist[lstptr + 1] = LASTTOK;

                my_symb[LASTTOK].pri = 10;

            } else {
                printf("Illegal use of delay2 \n");
                return 1;
            }
        }

        //        check for delshft symbol
        if (newtok == DELSHFTSYM) {
            temp = my_symb[toklist[lstptr + 1]].com;
            /* !! */ if (is_uvar(temp)) {
                // ram -- same issue
                my_symb[LASTTOK].com = COM(SVARTYPE, temp % MAXTYPE);  // create a temporary sybol
                NDELAYS++;
                toklist[lstptr + 1] = LASTTOK;

                my_symb[LASTTOK].pri = 10;

            } else {
                printf("Illegal use of delay2 Shift \n");
                return 1;
            }
        }

        if (newtok == SETSYM) {
            temp = my_symb[toklist[lstptr + 1]].com;
            if (is_uvar(temp)) {
                // ram -- same issue
                my_symb[LASTTOK].com = COM(SVARTYPE, temp % MAXTYPE);  // create a temporary sybol
                toklist[lstptr + 1] = LASTTOK;

                my_symb[LASTTOK].pri = 10;
            } else {
                printf("Illegal use of set - variables only\n");
                return 1;
            }
        }

        // check for shift
        if (newtok == SHIFTSYM || newtok == ISHIFTSYM) {
            temp = my_symb[toklist[lstptr + 1]].com;
            /* !! */ if (is_uvar(temp) || is_ucon(temp)) {
                // ram -- same issue
                if (is_uvar(temp)) {
                    my_symb[LASTTOK].com = COM(SVARTYPE, temp % MAXTYPE);
                }
                if (is_ucon(temp)) {
                    my_symb[LASTTOK].com = COM(SCONTYPE, temp % MAXTYPE);
                }
                // create a temporary sybol

                toklist[lstptr + 1] = LASTTOK;

                my_symb[LASTTOK].pri = 10;

            } else {
                printf("Illegal use of SHIFT \n");
                return 1;
            }
        }

    next:
        if ((newtok == ENDTOK) && (oldtok == STARTTOK)) {
            break;
        }

        if (newtok == LPAREN) {
            tokstak[tokptr] = LPAREN;
            tokptr++;
            oldtok = LPAREN;
            goto getnew;
        }
        if (newtok == RPAREN) {
            switch (oldtok) {
            case LPAREN:
                tokptr--;
                oldtok = tokstak[tokptr - 1];
                goto getnew;
            case COMMA:
                tokptr--;
                ncomma++;
                oldtok = tokstak[tokptr - 1];
                goto next;
            default:
                break;
            }
        }
        if ((newtok == COMMA) && (oldtok == COMMA)) {
            tokstak[tokptr] = COMMA;
            tokptr++;
            goto getnew;
        }
        // ram -- the THOUS problem
        /*           if(my_symb[oldtok%THOUS].pri>=my_symb[newtok%THOUS].pri)
        {
         command[comptr]=my_symb[oldtok%THOUS].com;
         if((my_symb[oldtok%THOUS].arg==2)&&
         (my_symb[oldtok%THOUS].com/MAXTYPE==FUN2TYPE)) */

        if (my_symb[oldtok].pri >= my_symb[newtok].pri) {
            command[comptr] = my_symb[oldtok].com;
            if ((my_symb[oldtok].arg == 2) && (my_symb[oldtok].com / MAXTYPE == FUN2TYPE)) {
                ncomma--;
            }
            my_com = command[comptr];
            comptr++;
            //   New code   3/95
            if (my_com == NUMSYM) {
                tokptr--;
                command[comptr] = tokstak[tokptr - 1];
                comptr++;
                tokptr--;
                command[comptr] = tokstak[tokptr - 1];
                comptr++;
            }
            //   end new code    3/95
            if (my_com == SUMSYM) {
                loopstk[lptr] = comptr;
                comptr++;
                lptr++;
                ncomma -= 1;
            }
            if (my_com == ENDSUM) {
                lptr--;
                jmp = comptr - loopstk[lptr] - 1;
                command[loopstk[lptr]] = jmp;
            }
            if (my_com == MYIF) {
                loopstk[lptr] = comptr;  // add some space for jump
                comptr++;
                lptr++;
                nif++;
            }
            if (my_com == MYTHEN) {
                // First resolve the if jump
                lptr--;
                jmp = comptr - loopstk[lptr];  // -1 is old
                command[loopstk[lptr]] = jmp;
                // Then set up for the then jump
                loopstk[lptr] = comptr;
                lptr++;
                comptr++;
                nthen++;
            }
            if (my_com == MYELSE) {
                lptr--;
                jmp = comptr - loopstk[lptr] - 1;
                command[loopstk[lptr]] = jmp;
                nelse++;
            }

            if (my_com == ENDDELAY || my_com == ENDSHIFT || my_com == ENDISHIFT) {
                ncomma -= 1;
            }
            if (my_com == ENDDELSHFT || my_com == ENDSET) {
                ncomma -= 2;
            }
            /*  if(my_com==CONV||my_com==DCONV){
               ncomma-=1;
              }  */

            //    CHECK FOR USER FUNCTION
            if (is_ufun(my_com)) {
                my_arg = my_symb[oldtok].arg;
                command[comptr] = my_arg;
                comptr++;
                ncomma = ncomma + 1 - my_arg;
            }
            //      USER FUNCTION OKAY
            tokptr--;
            oldtok = tokstak[tokptr - 1];
            goto next;
        }
        //    NEW code       3/95
        if (newtok == NUMTOK) {
            tokstak[tokptr++] = toklist[lstptr++];
            tokstak[tokptr++] = toklist[lstptr++];
        }
        //  end  3/95
        tokstak[tokptr] = newtok;
        oldtok = newtok;
        tokptr++;
        goto getnew;
    }
    if (ncomma != 0) {
        ggets_plintf("Illegal number of arguments\n");
        return 1;
    }
    if ((nif != nelse) || (nif != nthen)) {
        ggets_plintf("If statement missing ELSE or THEN \n");
        return 1;
    }
    command[comptr] = my_symb[ENDTOK].com;

    return 0;
}

void
pr_command(int32 *command) {
    int32 i = 0;
    int32 token;
    while (true) {
        token = command[i];
        ggets_plintf("%d %d \n", i, token);
        if (token == ENDEXP) {
            break;
        }
        i++;
    }
    return;
}

void
show_where(char *string, int32 index) {
    char junk[MAXEXPLEN];
    for (int32 i = 0; i < index; i++) {
        junk[i] = ' ';
    }
    junk[index] = '^';
    junk[index + 1] = 0;
    ggets_plintf("%s\n%s\n", string, junk);
    return;
}

int32
function_sym(int32 token) {
    // functions should have ( after them
    int32 com = my_symb[token].com;
    int32 i1 = com / MAXTYPE;

    if (i1 == FUN1TYPE && !unary_sym(token)) {
        return 1;  // single variable functions
    }
    if (i1 == FUN2TYPE && !binary_sym(token)) {
        return 1;  // two-variable function
    }
    /* ram this was: if(i1==UFUN||i1==7||i1==6||i1==5)return 1; recall: 5 was
     * bad
     */
    if (i1 == UFUNTYPE || i1 == TABTYPE || i1 == VECTYPE || i1 == NETTYPE) {
        return 1;
    }
    if (token == DELSHFTSYM || token == SETSYM || token == DELSYM || token == SHIFTSYM ||
        token == ISHIFTSYM || com == MYIF || com == MYTHEN || com == MYELSE || com == SUMSYM ||
        com == ENDSUM) {
        return 1;
    }
    return 0;
}

int32
unary_sym(int32 token) {
    // ram: these are tokens not byte code, so no change here?
    if (token == 9 || token == 55) {
        return 1;
    }
    return 0;
}

int32
binary_sym(int32 token) {
    // ram: these are tokens not byte code, so no change here?
    if (token > 2 && token < 9) {
        return 1;
    }
    if (token > 43 && token < 51) {
        return 1;
    }
    if (token == 54) {
        return 1;
    }
    return 0;
}

int32
pure_number(int32 token) {
    int32 com = my_symb[token].com;
    int32 i1 = com / MAXTYPE;
    /* !! */ if (token == NUMTOK || isvar(i1) || iscnst(i1) || isker(i1) || i1 == USTACKTYPE ||
                 token == INDX) {
        return 1;
    }
    return 0;
}

int32
gives_number(int32 token) {
    int32 com = my_symb[token].com;
    int32 i1 = com / MAXTYPE;
    if (token == INDX) {
        return 1;
    }
    if (token == NUMTOK) {
        return 1;
    }
    if (i1 == FUN1TYPE && !unary_sym(token)) {
        return 1;  // single variable functions
    }
    if (i1 == FUN2TYPE && !binary_sym(token)) {
        return 1;  // two-variable function
    }
    // !!
    /* ram: 5 issue; was
     * if(i1==8||isvar(i1)||iscnst(i1)||i1==7||i1==6||i1==5||isker(i1)||i1==UFUN)return
     * 1;
     */
    if (i1 == USTACKTYPE || isvar(i1) || iscnst(i1) || i1 == TABTYPE || i1 == VECTYPE ||
        i1 == NETTYPE || isker(i1) || i1 == UFUNTYPE) {
        return 1;
    }
    if (com == MYIF || token == DELSHFTSYM || token == SETSYM || token == DELSYM ||
        token == SHIFTSYM || token == ISHIFTSYM || com == SUMSYM) {
        return 1;
    }
    return 0;
}

int32
check_syntax(  // 1 is BAD!
    int32 oldtoken, int32 newtoken) {
    int32 com2 = my_symb[newtoken].com;

    /* if the first symbol or (  or binary symbol then must be unary symbol or
       something that returns a number or another (
    */

    if (unary_sym(oldtoken) || oldtoken == COMMA || oldtoken == STARTTOK || oldtoken == LPAREN ||
        binary_sym(oldtoken)) {
        if (unary_sym(newtoken) || gives_number(newtoken) || newtoken == LPAREN) {
            return 0;
        }
        return 1;
    }

    /* if this is a regular function, then better have (
     */

    if (function_sym(oldtoken)) {
        if (newtoken == LPAREN) {
            return 0;
        }
        return 1;
    }

    /* if we have a constant or variable or ) or kernel then better
       have binary symbol or "then" or "else" as next symbol
    */

    if (pure_number(oldtoken)) {
        if (binary_sym(newtoken) || newtoken == RPAREN || newtoken == COMMA || newtoken == ENDTOK) {
            return 0;
        }

        return 1;
    }

    if (oldtoken == RPAREN) {
        if (binary_sym(newtoken) || newtoken == RPAREN || newtoken == COMMA || newtoken == ENDTOK) {
            return 0;
        }
        if (com2 == MYELSE || com2 == MYTHEN || com2 == ENDSUM) {
            return 0;
        }

        return 1;
    }

    ggets_plintf("Bad token %d \n", oldtoken);
    return 1;
}

/******************************
 *    PARSER                   *
 ******************************/

int32
make_toks(char *dest, int32 *my_token) {
    char num[40];
    double value;
    int32 old_tok = STARTTOK;
    int32 tok_in = 0;
    int32 index = 0;
    int32 token;
    int32 nparen = 0;
    int32 lastindex = 0;
    union  //  WARNING  -- ASSUMES 32 bit int32  and 64 bit double
    {
        struct {
            int32 int1;
            int32 int2;
        } pieces;
        struct {
            double z;
        } num;
    } encoder;

    while (dest[index] != '\0') {
        lastindex = index;
        find_tok(dest, &index, &token);
        if ((token == MINUS) &&
            ((old_tok == STARTTOK) || (old_tok == COMMA) || (old_tok == LPAREN))) {
            token = NEGATE;
        }
        if (token == LPAREN) {
            ++nparen;
        }
        if (token == RPAREN) {
            --nparen;
        }

        if (token == NSYM) {
            if (do_num(dest, num, &value, &index)) {
                show_where(dest, index);
                return 1;
            }
            //    new code        3/95
            encoder.num.z = value;
            my_token[tok_in++] = NUMTOK;
            my_token[tok_in++] = encoder.pieces.int1;
            my_token[tok_in++] = encoder.pieces.int2;
            if (check_syntax(old_tok, NUMTOK) == 1) {
                ggets_plintf("Illegal syntax \n");
                show_where(dest, lastindex);
                return 1;
            }
            old_tok = NUMTOK;

        }

        else {
            my_token[tok_in++] = token;
            if (check_syntax(old_tok, token) == 1) {
                ggets_plintf("Illegal syntax (Ref:%d %d) \n", old_tok, token);
                show_where(dest, lastindex);
                tokeninfo(old_tok);
                tokeninfo(token);
                return 1;
            }

            old_tok = token;
        }
    }

    my_token[tok_in++] = ENDTOK;
    if (check_syntax(old_tok, ENDTOK) == 1) {
        ggets_plintf("Premature end of expression \n");
        show_where(dest, lastindex);
        return 1;
    }
    if (nparen != 0) {
        if (errout) {
            printf(" parentheses don't match\n");
        }
        return 1;
    }
    return 0;
}

void
tokeninfo(int32 tok) {
    ggets_plintf(" %s %d %d %d %d \n", my_symb[tok].name, my_symb[tok].len, my_symb[tok].com,
                 my_symb[tok].arg, my_symb[tok].pri);
    return;
}

int32
do_num(char *source, char *num, double *value, int32 *ind) {
    int32 j = 0;
    int32 i = *ind;
    int32 error = 0;
    int32 ndec = 0;
    int32 nexp = 0;
    int32 ndig = 0;
    char ch;
    char oldch;
    oldch = '\0';
    *value = 0.0;
    while (true) {
        ch = source[i];
        if (((ch == '+') || (ch == '-')) && (oldch != 'E')) {
            break;
        }
        if ((ch == '*') || (ch == '^') || (ch == '/') || (ch == ',') || (ch == ')') ||
            (ch == '\0') || (ch == '|') || (ch == '>') || (ch == '<') || (ch == '&') ||
            (ch == '=')) {
            break;
        }
        if ((ch == 'E') || (ch == '.') || (ch == '+') || (ch == '-') || isdigit(ch)) {
            if (isdigit(ch)) {
                ndig++;
            }
            switch (ch) {
            case 'E':
                nexp++;
                if ((nexp == 2) || (ndig == 0)) {
                    goto err;
                }
                break;
            case '.':
                ndec++;
                if ((ndec == 2) || (nexp == 1)) {
                    goto err;
                }
                break;
            default:
                break;
            }
            num[j] = ch;
            j++;
            i++;
            oldch = ch;
        } else {
        err:
            num[j] = ch;
            j++;
            error = 1;
            break;
        }
    }
    num[j] = '\0';
    if (error == 0) {
        *value = atof(num);
    } else if (errout) {
        printf(" illegal expression: %s\n", num);
    }
    *ind = i;
    return error;
}

void
convert(char *source, char *dest) {
    char ch;
    int32 i = 0;
    int32 j = 0;
    while (true) {
        ch = source[i];
        if (!isspace(ch)) {
            dest[j++] = ch;
        }
        i++;
        if (ch == '\0') {
            break;
        }
    }
    strupr(dest);
    return;
}

void
find_tok(char *source, int32 *index, int32 *tok) {
    int32 i = *index;
    int32 maxlen = 0;
    int32 symlen;
    int32 my_tok;
    int32 match;
    my_tok = NSYM;
    for (int32 k = 0; k < NSYM; k++) {
        symlen = my_symb[k].len;
        if (symlen <= maxlen) {
            continue;
        }

        match = 1;
        for (int32 j = 0; j < symlen; j++) {
            if (source[i + j] != my_symb[k].name[j]) {
                match = 0;
                break;
            }
        }
        if (match != 0) {
            my_tok = k;
            maxlen = symlen;
        }
    }
    *index = *index + maxlen;
    *tok = my_tok;
    return;
}

double
pmod(double x, double y) {
    double z = fmod(x, y);
    if (z < 0) {
        z += y;
    }
    return z;
}

/*   These are the Bessel Functions; if you dont have them then
     return some sort of dummy value or else write a program
     to compute them
*/

double
bessel_j(double x, double y) {
    int32 n = (int32)x;
    return jn(n, y);
}

double
bessel_y(double x, double y) {
    int32 n = (int32)x;
    return yn(n, y);
}

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

double
bessi(double nn, double x) {
    int32 n;
    double bi;
    double bim;
    double bip;
    double tox;
    double ans;
    n = (int32)nn;
    if (n == 0) {
        return bessi0(x);
    }
    if (n == 1) {
        return bessi1(x);
    }
    if (x == 0.0) {
        return 0.0;
    } else {
        tox = 2.0 / fabs(x);
        bip = ans = 0.0;
        bi = 1.0;
        for (int32 j = 2*(n + (int32)(sqrt(ACC*n))); j > 0; j--) {
            bim = bip + j*tox*bi;
            bip = bi;
            bi = bim;
            if (fabs(bi) > BIGNO) {
                ans *= BIGNI;
                bi *= BIGNI;
                bip *= BIGNI;
            }
            if (j == n) {
                ans = bip;
            }
        }
        ans *= bessi0(x) / bi;
        return x < 0.0 && (n & 1) ? -ans : ans;
    }
}

double
bessi0(double x) {
    double ax;
    double ans;
    double y;

    if ((ax = fabs(x)) < 3.75) {
        y = x / 3.75;
        y *= y;
        ans = 1.0 +
              y*(3.5156229 +
                   y*(3.0899424 +
                        y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
    } else {
        y = 3.75 / ax;
        ans =
            (exp(ax) / sqrt(ax)) *
            (0.39894228 +
             y*(0.1328592e-1 +
                  y*(0.225319e-2 + y*(-0.157565e-2 +
                                          y*(0.916281e-2 + y*(-0.2057706e-1 +
                                                                  y*(0.2635537e-1 +
                                                                       y*(-0.1647633e-1 +
                                                                            y*0.392377e-2))))))));
    }
    return ans;
}

double
bessi1(double x) {
    double ax;
    double ans;
    double y;

    if ((ax = fabs(x)) < 3.75) {
        y = x / 3.75;
        y *= y;
        ans = ax*(0.5 + y*(0.87890594 +
                               y*(0.51498869 +
                                    y*(0.15084934 + y*(0.2658733e-1 +
                                                           y*(0.301532e-2 + y*0.32411e-3))))));
    } else {
        y = 3.75 / ax;
        ans = 0.2282967e-1 + y*(-0.2895312e-1 + y*(0.1787654e-1 - y*0.420059e-2));
        ans = 0.39894228 +
              y*(-0.3988024e-1 +
                   y*(-0.362018e-2 + y*(0.163801e-2 + y*(-0.1031555e-1 + y*ans))));
        ans *= (exp(ax) / sqrt(ax));
    }
    return x < 0.0 ? -ans : ans;
}

double
bessis(double nn, double x) {
    int32 n;
    double bi;
    double bim;
    double bip;
    double tox;
    double ans;
    n = (int32)nn;
    if (n == 0) {
        return bessis0(x);
    }
    if (n == 1) {
        return bessis1(x);
    }
    if (x == 0.0) {
        return 0.0;
    } else {
        tox = 2.0 / fabs(x);
        bip = ans = 0.0;
        bi = 1.0;
        for (int32 j = 2*(n + (int32)(sqrt(ACC*n))); j > 0; j--) {
            bim = bip + j*tox*bi;
            bip = bi;
            bi = bim;
            if (fabs(bi) > BIGNO) {
                ans *= BIGNI;
                bi *= BIGNI;
                bip *= BIGNI;
            }
            if (j == n) {
                ans = bip;
            }
        }
        ans *= bessis0(x) / bi;
        return x < 0.0 && (n & 1) ? -ans : ans;
    }
}

double
bessis0(double x) {
    double ax;
    double ans;
    double y;

    if ((ax = fabs(x)) < 3.75) {
        y = x / 3.75;
        y *= y;
        ans = (1.0 +
               y*(3.5156229 +
                    y*(3.0899424 +
                         y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))))) *
              exp(-ax);
    } else {
        y = 3.75 / ax;
        ans =
            (1.0 / sqrt(ax)) *
            (0.39894228 +
             y*(0.1328592e-1 +
                  y*(0.225319e-2 + y*(-0.157565e-2 +
                                          y*(0.916281e-2 + y*(-0.2057706e-1 +
                                                                  y*(0.2635537e-1 +
                                                                       y*(-0.1647633e-1 +
                                                                            y*0.392377e-2))))))));
    }
    return ans;
}

double
bessis1(double x) {
    double ax;
    double ans;
    double y;

    if ((ax = fabs(x)) < 3.75) {
        y = x / 3.75;
        y *= y;
        ans = exp(-ax)*ax *
              (0.5 + y*(0.87890594 +
                          y*(0.51498869 +
                               y*(0.15084934 +
                                    y*(0.2658733e-1 + y*(0.301532e-2 + y*0.32411e-3))))));
    } else {
        y = 3.75 / ax;
        ans = 0.2282967e-1 + y*(-0.2895312e-1 + y*(0.1787654e-1 - y*0.420059e-2));
        ans = 0.39894228 +
              y*(-0.3988024e-1 +
                   y*(-0.362018e-2 + y*(0.163801e-2 + y*(-0.1031555e-1 + y*ans))));
        ans *= (1. / sqrt(ax));
    }
    return x < 0.0 ? -ans : ans;
}

#undef ACC
#undef BIGNO
#undef BIGNI

/*********************************************
          FANCY delay2 HERE                   *-------------------------<<<
*********************************************/

char *
com_name(int32 com) {
    int32 i;
    for (i = 0; i < NSYM; i++) {
        if (my_symb[i].com == com) {
            break;
        }
    }
    if (i < NSYM) {
        return my_symb[i].name;
    } else {
        return "";
    }
}

double
do_shift(double shift, double variable) {
    int32 it;
    int32 in;
    int32 i = (int32)(variable);
    int32 ish = (int32)shift;

    if (i < 0) {
        return 0.0;
    }
    it = i / MAXTYPE;
    in = (i % MAXTYPE) + ish;
    switch (it) {
    case CONTYPE:
        if (in > NCON) {
            return 0.0;
        } else {
            return constants[in];
        }
    case VARTYPE:
        if (in > MAX_ODE) {
            return 0.0;
        } else {
            return variables[in];
        }
    default:
        ggets_plintf("This can't happen: Invalid symbol index for SHIFT: i = %d\n", i);
        return 0.0;
    }
}

double
do_ishift(double shift, double variable) {
    return variable + shift;
}

double
do_delay_shift(double delay2, double shift, double variable) {
    int32 in;
    int32 i = (int32)(variable);
    int32 ish = (int32)shift;
    if (i < 0) {
        return 0.0;
    }
    in = (i % MAXTYPE) + ish;

    if (in > MAX_ODE) {
        return 0.0;
    }

    if (del_stab_flag > 0) {
        if (delay_flag && delay2 > 0.0) {
            return delay_handle_get_delay(in - 1, delay2);
        }
        return variables[in];
    }

    return delay_handle_stab_eval(delay2, in);
}

double
do_delay(double delay2, double i) {
    int32 variable;
    /* ram - this was a little weird, since i is a double... except I think it's
     * secretely an int64 */
    variable = ((int32)i) % MAXTYPE;

    if (del_stab_flag > 0) {
        if (delay_flag && delay2 > 0.0) {
            return delay_handle_get_delay(variable - 1, delay2);
        }
        return variables[variable];
    }

    return delay_handle_stab_eval(delay2, (int32)variable);
}

double
hom_bcs(int32 i) {
    (void)i;
    return 0.0;  // this is deprecated so no longer used
}

double
normal(double mean, double std) {
    double fac;
    double r;
    double v1;
    double v2;
    if (box_muller_flag == 0) {
        do {
            v1 = 2.0*markov_ndrand48() - 1.0;
            v2 = 2.0*markov_ndrand48() - 1.0;
            r = v1*v1 + v2*v2;
        } while (r >= 1.0);
        fac = sqrt(-2.0*log(r) / r);
        box_muller = v1*fac;
        box_muller_flag = 1;
        return v2*fac*std + mean;
    } else {
        box_muller_flag = 0;
        return box_muller*std + mean;
    }
}

double
max(double x, double y) {
    return ((x > y) ? x : y);
}

double
min(double x, double y) {
    return ((x < y) ? x : y);
}

double
neg(double z) {
    return -z;
}

double
recip(double z) {
    return 1.00 / z;
}

double
heaviside(double z) {
    double w = 1.0;
    if (z < 0) {
        w = 0.0;
    }
    return w;
}

double
rndom(double z) {
    return z*markov_ndrand48();
}

double
signum(double z) {
    if (z < 0.0) {
        return -1.0;
    }
    if (z > 0.0) {
        return 1.0;
    }
    return 0.0;
}

/*  logical stuff  */

double
dnot(double x) {
    return (double)(x == 0.0);
}

double
dand(double x, double y) {
    return (double)(fabs(x) >= 1e-13 && fabs(y) >= 1e-13);
}

double
dor(double x, double y) {
    return (double)(fabs(x) >= 1e-13 || fabs(y) >= 1e-13);
}

double
dge(double x, double y) {
    return (double)(x >= y);
}

double
dle(double x, double y) {
    return (double)(x <= y);
}

double
deq(double x, double y) {
    return (double)(x == y);
}

double
dne(double x, double y) {
    return (double)(x != y);
}

double
dgt(double x, double y) {
    return (double)(x > y);
}

double
dlt(double x, double y) {
    return (double)(x < y);
}

/*              end of logical stuff    */

double
evaluate(int32 *equat) {
    uptr = 0;
    stack_pointer = 0;
    return eval_rpn(equat);
}

double
eval_rpn(int32 *equat) {
    int32 i, it, in, j, *tmpeq;
    int32 is;

    int32 low;
    int32 high;
    int32 ijmp;
    int32 iv;
    double temx;
    double temy;
    double temz;
    double sum;
    union  //  WARNING  -- ASSUMES 32 bit int32  and 64 bit double
    {
        struct {
            int32 int1;
            int32 int2;
        } pieces;
        struct {
            double z;
        } num;
    } encoder;

    while ((i = *equat++) != ENDEXP) {
        switch (i) {
        case NUMSYM:
            encoder.pieces.int2 = *equat++;
            encoder.pieces.int1 = *equat++;
            PUSH(encoder.num.z);
            break;
        case ENDFUN:
            i = *equat++;

            uptr -= i;

            break;

        case MYIF:
            temx = POP;
            ijmp = *equat++;
            if (temx == 0.0) {
                equat += ijmp;
            }
            break;
        case MYTHEN:
            ijmp = *equat++;
            equat += ijmp;
            break;
        case MYELSE:
            break;

        case ENDDELSHFT:
            temx = POP;
            temy = POP;
            temz = POP;
            PUSH(do_delay_shift(temx, temy, temz));
            break;
        case ENDSET: /* indirectly set a variable + shift to a value
                        SET(name,shift,value)  */
            temx = POP;
            temy = POP;
            temz = POP;
            iv = (int32)temy + (((int32)temz) % MAXTYPE);
            variables[iv] = temx;
            PUSH(temx);
            break;
        case ENDDELAY:
            temx = POP;
            temy = POP;

            PUSH(do_delay(temx, temy));
            break;

        case ENDSHIFT:
            temx = POP;
            temy = POP;
            PUSH(do_shift(temx, temy));
            break;
        case ENDISHIFT:
            temx = POP;
            temy = POP;
            PUSH(do_ishift(temx, temy));
            break;
        case SUMSYM:
            temx = POP;
            high = (int32)temx;
            temx = POP;
            low = (int32)temx;
            ijmp = *equat++;
            sum = 0.0;
            if (low <= high) {
                for (is = low; is <= high; is++) {
                    tmpeq = equat;
                    constants[SumIndex] = (double)is;
                    sum += eval_rpn(tmpeq);
                }
            }
            equat += ijmp;
            PUSH(sum);
            break;

        case ENDSUM:
            return POP;
        case INDXCOM:
            PUSH(current_index);
            break;
        default: {
            it = i / MAXTYPE;
            in = i % MAXTYPE;
            switch (it) {
            case FUN1TYPE:
                PUSH(((FunctionDouble)fun1[in])(POP));
                break;
            case FUN2TYPE: {
                if (in == 0) {
                    temx = POP;
                    temy = POP;
                    PUSH(temx + temy);
                    goto bye;
                }
                if (in == 2) {
                    temx = POP;
                    temy = POP;
                    PUSH(temx*temy);
                    goto bye;
                }
                if (in == 1) {
                    temx = POP;
                    temy = POP;
                    PUSH(temy - temx);
                    goto bye;
                }
                if (in == 3) {
                    temx = POP;
                    if (temx == 0.0) {
                        temx = DOUB_EPS;
                    }
                    temy = POP;
                    PUSH(temy / temx);
                    goto bye;
                }
                temx = POP;
                temy = POP;
                PUSH(((FunctionDoubleDouble)fun2[in])(temy, temx));
                break;
            }
            case CONTYPE:
                PUSH(constants[in]);
                break;
            case VECTYPE:
                PUSH(simplenet_vector_value(POP, in));
                break;
            case NETTYPE:
                PUSH(simplenet_network_value(POP, in));
                break;
            case TABTYPE:
                PUSH(tabular_lookup(POP, in));
                break;

            case USTACKTYPE:
                /* ram: so this means ustacks really do need to be of USTACKTYPE
                 */
                PUSH(ustack[uptr - 1 - in]);
                break;
            case KERTYPE:
                PUSH(volterra_ker_val(in));
                break;
            case VARTYPE:
                PUSH(variables[in]);
                break;

                // indexes for shift and delay2 operators...
            case SCONTYPE:

                PUSH((double)(COM(CONTYPE, in)));
                break;
            case SVARTYPE:

                PUSH((double)(COM(VARTYPE, in)));
                break;

            case UFUNTYPE:
                i = *equat++;

                for (j = 0; j < i; j++) {
                    ustack[uptr] = POP;

                    uptr++;
                }
                PUSH(eval_rpn(ufun[in]));
                break;
            default:
                break;
            }
        bye:
            j = 0;
        }
        }
    }
    return POP;
}

/* code for log-gamma if you dont have it */

/*  STRING STUFF  */
#ifndef STRUPR
void
strupr(char *s) {
    int32 i = 0;
    while (s[i]) {
        if (islower(s[i])) {
            s[i] -= 32;
        }
        i++;
    }
    return;
}

void
strlwr(char *s) {
    int32 i = 0;
    while (s[i]) {
        if (isupper(s[i])) {
            s[i] += 32;
        }
        i++;
    }
    return;
}
#endif
