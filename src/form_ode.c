#include "functions.h"
#include "integers.h"
#include <stdbool.h>

#include "parserslow.h"
#include "read_dir.h"
#include <unistd.h>

#include "newpars.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#include "xpplim.h"

#include "newpars.h"

#define MAXONLY 1000

static int32 IN_INCLUDED_FILE = 0;

char uvar_names[MAX_ODE][MAX_ODE_NAME_LENGTH];
char upar_names[MAX_PAR][MAX_ODE_NAME_LENGTH];
char *ode_names[MAX_ODE];
char *save_eqn[MAXLINES];
double default_val[MAX_PAR];

typedef struct VarInfo {
    char lhs[MAXEXPLEN];
    char rhs[MAXEXPLEN];
    char args[MAXARG][NAMLEN + 1];
    int32 type;
    int32 nargs;
    double value;
    struct VarInfo *next, *prev;
} VarInfo;

static VarInfo *my_varinfo;
static int32 start_var_info = 0;

int32 *my_ode[MAX_ODE];

int32 leng[MAX_ODE];

static char *onlylist[MAXONLY];
int32 *plotlist;
static int32 N_only = 0;
int32 N_plist;

Action comments[MAXCOMMENTS];
static int32 is_a_map = 0;
int32 n_comments = 0;
BcStruct my_bc[MAX_ODE];

double default_ic[MAX_ODE];
int32 NODE;
int32 NUPAR;
int32 NLINES;
int32 PrimeStart;
int32 NCON_START;
int32 NSYM_START;
int32 BVP_N;

int32 ConvertStyle = 0;
FILE *convertf;
static int32 OldStyle = 1;
int32 IN_VARS;
int32 NMarkov;

int32 FIX_VAR;

int32 EqType[MAX_ODE];
static int32 Naux = 0;
static char aux_names[MAX_ODE][12];

FixInfo fixinfo[MAX_ODE];

static void form_ode_strncpy_trim(char *dest, char *source, int32 n);
static void form_ode_strcpy_trim(char *dest, char *source);
static void form_ode_advance_past_first_word(char **sptr);
static void form_ode_add_comment(char *s);
static int32 form_ode_is_comment(char *s);
static int32 form_ode_not_ker(char *s, int32 i);
static int32 form_ode_check_if_ic(char *big);
static void form_ode_read_a_line(FILE *fp, char *s);
static void form_ode_remove_blanks(char *s1);
static int32 form_ode_next_nonspace(char *s1, int32 i0, int32 *i1);
static int32 form_ode_strparse(char *s1, char *s2, int32 i0, int32 *i1);
static int32 form_ode_extract(char *s1, int32 *ie, int32 i1);
static void form_ode_init_varinfo(void);
static int32 form_ode_parse_a_string(char *s1, VarInfo *v);
static void form_ode_strpiece(char *dest, char *src, int32 i0, int32 ie);
static int32 form_ode_formula_or_number(char *expr, double *z);
static int32 form_ode_find_the_name(char list[MAX_ODE1][MAXVNAM], int32 n,
                                    char *name);
static void form_ode_break_up_list(char *rhs);
static void form_ode_add_only(char *s);
static int32 form_ode_do_new_parser(FILE *fp, char *first, int32 nnn);
static int32 form_ode_if_end_include(char *old);
static int32 form_ode_if_include_file(char *old, char *nf);
static void form_ode_clrscr(void);
static void form_ode_find_ker(char *string, int32 *alt);
static void form_ode_take_apart(char *bob, double *value, char *name);
static void form_ode_list_em(char *wild);
static int32 form_ode_get_a_filename(char *filename, char *wild);
static void form_ode_format_list(char **s, int32 n);

int32
form_ode_make_eqn(void) {
    FILE *fptr;
    char wild[256];
    char string[256];
    int32 okay;
    NEQ = 2;
    FIX_VAR = 0;
    NMarkov = 0;

    okay = 0;
    snprintf(wild, sizeof(wild), "*.ode");
    form_ode_get_a_filename(string, wild);
    if ((fptr = fopen(string, "r")) == NULL) {
        ggets_plintf("\n Cannot open %s \n", string);
        return 0;
    }
    strcpy(this_file, string);
    form_ode_clrscr();
    okay = form_ode_get_eqn(fptr);
    fclose(fptr);
    return okay;
}

void
form_ode_strip_saveqn(void) {
    for (int32 i = 0; i < NLINES; i++) {
        for (int32 j = 0; j < (int32)strlen(save_eqn[i]); j++) {
            if (save_eqn[i][j] < 32) {
                save_eqn[i][j] = 32;
            }
        }
    }
    return;
}

int32
form_ode_idsc(char *string) {
    char c;
    int32 i = 0;
    int32 l = (int32)strlen(string);
    int32 j = 0;
    int32 flag = 0;
    char end[256];
    if (is_a_map == 1) {
        return 1;
    }
    while (i < l) {
        c = string[i];
        if (flag == 1) {
            end[j] = c;
            j++;
        }
        if (c == '.') {
            flag = 1;
        }
        i++;
    }
    end[j] = 0;

    if (strcmp(end, "dis") == 0 || strcmp(end, "dif") == 0) {
        return 1;
    }
    return 0;
}

void
form_ode_format_list(char **s, int32 n) {
    int32 ip;
    int32 ncol;
    int32 k;
    int32 j;
    char fmat[30];
    int32 lmax = 0;
    int32 l = 0;
    for (int32 i = 0; i < n; i++) {
        l = (int32)strlen(s[i]);
        if (lmax < l) {
            lmax = l;
        }
    }
    ncol = 80 / (lmax + 2);
    if (ncol < 1) {
        ncol = 1;
    }
    if (ncol > 8) {
        ncol = 8;
    }
    k = n / ncol;
    j = n - ncol*k;
    snprintf(fmat, sizeof(fmat), "%s%d%s", "%", lmax + 2, "s");
    for (ip = 0; ip < k; ip++) {
        for (int32 i = 0; i < ncol; i++) {
            ggets_plintf(fmat, s[ip*ncol + i]);
        }
        ggets_plintf("\n");
    }
    for (int32 i = 0; i < j; i++) {
        ggets_plintf(fmat, s[k*ncol + i]);
    }
    ggets_plintf("\n");
    return;
}

int32
form_ode_get_a_filename(char *filename, char *wild) {
    if (xpp_batch) {
        char string[MAXEXPLEN];
        form_ode_list_em(wild);
        while (true) {
            ggets_plintf("(r)un (c)d (l)ist ");
            scanf("%s", string);
            if (string[0] == 'r') {
                ggets_plintf("Run file: ");
                scanf("%s", filename);
                ggets_plintf("Loading %s\n ", filename);
                return 1;
            } else {
                if (string[0] == 'l') {
                    ggets_plintf("List files of type: ");
                    scanf("%s", wild);
                    form_ode_list_em(wild);
                } else {
                    if (string[0] == 'c') {
                        ggets_plintf("Change to directory: ");
                        scanf("%s", string);
                        read_dir_change_dir(string);
                        form_ode_list_em(wild);
                    }
                }
            }
        }
    } else {
        int32 status;
        int32 m;
        read_dir_get_directory(filename);
        m = (int32)strlen(filename);
        if (filename[m - 1] != '/') {
            strcat(filename, "/");
        }
        status = init_conds_file_selector("Select an ODE file", filename, wild);
        if (status == 0) {
            main_bye_bye();
        } else {
            return 1;
        }
    }
}

void
form_ode_list_em(char *wild) {
    read_dir_get_directory(cur_dir);
    ggets_plintf("%s: \n", cur_dir);
    read_dir_get_fileinfo(wild, cur_dir, &my_ff);
    ggets_plintf("DIRECTORIES:\n");
    form_ode_format_list(my_ff.dirnames, my_ff.ndirs);
    ggets_plintf("FILES OF TYPE %s:\n", wild);
    form_ode_format_list(my_ff.filenames, my_ff.nfiles);

    read_dir_free_finfo(&my_ff);
    return;
}

int32
form_ode_get_eqn(FILE *fptr) {
    char bob[MAXEXPLEN];
    char filename[XPP_MAX_NAME + 4];
    int32 done = 1;
    int32 nn;
    int32 i;
    int32 flag;
    init_rpn();
    NLINES = 0;
    IN_VARS = 0;
    NODE = 0;
    BVP_N = 0;
    NUPAR = 0;
    NWiener = 0;
    /* load_eqn_check_for_xpprc();  This is now done just once and in
     * main_do_vis_env() */
    strcpy(options, "default.opt");
    parserslow_add_var("t", 0.0);
    fgets(bob, MAXEXPLEN, fptr);
    nn = (int32)strlen(bob) + 1;
    if (NLINES > MAXLINES) {
        fprintf(stderr, "whoops! NLINES>MAXLINES in form_ode.c ...\n");
        exit(1);
    }
    if ((save_eqn[NLINES] = xmalloc((usize)nn)) == NULL) {
        ggets_plintf("Out of memory...");
        exit(0);
    }

    strncpy(save_eqn[NLINES++], bob, (usize)nn);
    i = atoi(bob);
    if (i <= 0) { /* New parser ---   */

        OldStyle = 0;
        ConvertStyle = 0;
        flag = form_ode_do_new_parser(fptr, bob, 0);
        if (flag < 0) {
            exit(0);
        }
    } else {
        OldStyle = 1;
        NEQ = i;
        ggets_plintf("NEQ=%d\n", NEQ);
        if (ConvertStyle) {
            if (strlen(this_file) == 0) {
                strcpy(filename, "convert.ode");
            } else {
                snprintf(filename, sizeof(filename), "%s.new", this_file);
            }
            if ((convertf = fopen(filename, "w")) == NULL) {
                printf(" Cannot open %s - no conversion done \n", filename);
                ConvertStyle = 0;
            }
            fprintf(convertf, "# converted %s \n", this_file);
        }
        while (done) {
            fgets(bob, MAXEXPLEN, fptr);
            nn = (int32)strlen(bob) + 1;
            if ((save_eqn[NLINES] = xmalloc((usize)nn)) == NULL) {
                exit(0);
            }
            strncpy(save_eqn[NLINES++], bob, (usize)nn);
            done = form_ode_compiler(bob, fptr);
        }
        if (ConvertStyle) {
            fprintf(convertf, "done\n");
            fclose(convertf);
        }
    }
    if ((NODE + NMarkov) == 0) {
        ggets_plintf(
            " Must have at least one equation! \n Probably not an ODE file.\n");
        exit(0);
    }
    if (BVP_N > IN_VARS) {
        ggets_plintf("Too many boundary conditions\n");
        exit(0);
    }
    if (BVP_N < IN_VARS) {
        if (BVP_N > 0) {
            printf("Warning: Too few boundary conditions\n");
        }
        for (i = BVP_N; i < IN_VARS; i++) {
            my_bc[i].com = xmalloc(200*sizeof(*(my_bc[i].com)));
            my_bc[i].string = xmalloc(256);
            my_bc[i].name = xmalloc(10);
            my_bc[i].side = 0;
            strcpy(my_bc[i].string, "0");
            strcpy(my_bc[i].name, "0=");
        }
    }
    BVP_FLAG = 1;

    if (NODE != NEQ + FIX_VAR - NMarkov) {
        ggets_plintf(" Too many/few equations\n");
        exit(0);
    }
    if (IN_VARS > NEQ) {
        ggets_plintf(" Too many variables\n");
        exit(0);
    }
    NODE = IN_VARS;

    for (i = 0; i < Naux; i++) {
        strcpy(uvar_names[i + NODE + NMarkov], aux_names[i]);
    }

    for (i = NODE + NMarkov + Naux; i < NEQ; i++) {
        snprintf(uvar_names[i], sizeof(uvar_names[i]), "AUX%d",
                 i - NODE - NMarkov + 1);
    }

    for (i = 0; i < NEQ; i++) {
        strupr(uvar_names[i]);
        strupr(ode_names[i]);
        ani_de_space(ode_names[i]);
    }
    /* add primed variables */
    PrimeStart = NVAR;
    if (NVAR < MAX_PRIME_VAR) {
        parserslow_add_var("t'", 0.0);
        for (i = 0; i < NODE; i++) {
            char prim[sizeof(uvar_names[i]) + 1];
            snprintf(prim, sizeof(prim), "%s'", uvar_names[i]);
            parserslow_add_var(prim, 0.0);
        }
    } else {
        ggets_plintf(
            " Warning: primed variables not added must have < %d variables\n",
            MAX_PRIME_VAR);
        ggets_plintf(" Averaging and boundary value problems cannot be done\n");
    }
    if (NMarkov > 0) {
        markov_compile_all();
    }
    if (flags_compile() == 1) {
        ggets_plintf(" Error in compiling a flag \n");
        exit(0);
    }
    flags_show();
    /* add auxiliary variables */
    for (i = NODE + NMarkov; i < NEQ; i++) {
        parserslow_add_var(uvar_names[i], 0.0);
    }
    NCON_START = NCON;
    NSYM_START = NSYM;
    xpp_version_maj = (double)MAJOR_VERSION;
    xpp_version_min = (double)MINOR_VERSION;
    ggets_plintf("Used %d constants and %d symbols \n", NCON, NSYM);
    ggets_plintf("XPPAUT %g.%g Copyright (C) 2002-now  Bard Ermentrout \n",
                 xpp_version_maj, xpp_version_min);
    return 1;
}

int32
form_ode_compiler(char *bob, FILE *fptr) {
    double value;
    double xlo;
    double xhi;
    int32 narg;
    int32 done;
    int32 nn;
    int32 iflg = 0;
    int32 VFlag = 0;
    int32 nstates;
    int32 alt;
    int32 index;
    int32 sign;
    char *ptr, *my_string, *command;
    char name[20];
    char formula[MAXEXPLEN];
    char condition[MAXEXPLEN];
    char fixname[MAX_ODE1][12];
    int32 nlin;
    int32 i;
    ptr = bob;
    done = 1;
    if (bob[0] == '@') {
        load_eqn_stor_internopts(bob);
        if (ConvertStyle) {
            fprintf(convertf, "%s\n", bob);
        }
        return done;
    }
    command = form_ode_get_first(ptr, " ,");
    strlwr(command);
    switch (*command) {
    case 'd':
        done = 0;
        break;
    case 's':
        /* form ode show syms */
        ggets_plintf("(    ,    )    +    -      *    ^    **    / \n");
        ggets_plintf("sin  cos  tan  atan  atan2 acos asin\n");
        ggets_plintf("exp  ln   log  log10 tanh  cosh sinh \n");
        ggets_plintf("max  min  heav flr   mod   sign sqrt \n");
        ggets_plintf("t    pi   ran  \n");
        break;
    case 'h':
        /* form ode welcome */
        ggets_plintf("\n The commands are: \n");
        ggets_plintf(" P(arameter) -- declare parameters "
                     "<name1>=<value1>,<name2>=<value2>,...\n");
        ggets_plintf(" F(ixed)     -- declare fixed variables\n");
        ggets_plintf(" V(ariables) -- declare ode variables \n");
        ggets_plintf(" U(ser)      -- declare user functions <name> <nargs> "
                     "<formula>\n");
        ggets_plintf(" C(hange)    -- change option file   <filename>\n");
        ggets_plintf(" O(de)       -- declare RHS for equations\n");
        ggets_plintf(" D(one)      -- finished compiling formula\n");
        ggets_plintf(" H(elp)      -- this menu                 \n");
        ggets_plintf(" S(ymbols)   -- Valid functions and symbols\n");
        ggets_plintf(" I(ntegral)  -- rhs for integral eqn\n");
        ggets_plintf(" K(ernel)    -- declare kernel for integral eqns\n");
        ggets_plintf(" T(able)     -- lookup table\n");
        ggets_plintf(" A(ux)       -- name auxiliary variable\n");
        ggets_plintf(" N(umbers)   --  hidden parameters\n");
        ggets_plintf(" M(arkov)    --  Markov variables \n");
        ggets_plintf(" W(iener)    -- Wiener parameter \n");
        ggets_plintf("_________________________________________________________"
                     "____________"
                     "____\n");
        break;
    case 'x':
        my_string = form_ode_do_fit_get_next("{ ");
        strcpy(condition, my_string);

        my_string = form_ode_do_fit_get_next("}\n");
        strcpy(formula, my_string);
        load_eqn_add_intern_set(condition, formula);
        break;
    case 'w': /*  Make a Wiener (heh heh) constants  */
        ggets_plintf("Wiener constants\n");
        if (ConvertStyle) {
            fprintf(convertf, "wiener ");
        }
        form_ode_advance_past_first_word(&ptr);
        while ((my_string = form_ode_get_next2(&ptr)) != NULL) {
            form_ode_take_apart(my_string, &value, name);
            free(my_string);
            ggets_plintf("|%s|=%f ", name, value);
            if (ConvertStyle) {
                fprintf(convertf, "%s  ", name);
            }
            if (parserslow_add_con(name, value)) {
                ggets_plintf("ERROR at line %d\n", NLINES);
                exit(0);
            }
            markov_add_wiener(NCON - 1);
        }
        if (ConvertStyle) {
            fprintf(convertf, "\n");
        }
        ggets_plintf("\n");
        break;
    case 'n':
        ggets_plintf(" Hidden params:\n");
        if (ConvertStyle) {
            fprintf(convertf, "number ");
        }

        form_ode_advance_past_first_word(&ptr);
        while ((my_string = form_ode_get_next2(&ptr)) != NULL) {
            form_ode_take_apart(my_string, &value, name);
            free(my_string);
            if (ConvertStyle) {
                fprintf(convertf, "%s=%g  ", name, value);
            }

            ggets_plintf("|%s|=%f ", name, value);
            if (parserslow_add_con(name, value)) {
                ggets_plintf("ERROR at line %d\n", NLINES);
                exit(0);
            }
        }
        if (ConvertStyle) {
            fprintf(convertf, "\n");
        }
        ggets_plintf("\n");
        break;
    case 'g': /* global */
        my_string = form_ode_do_fit_get_next("{ ");
        sign = atoi(my_string);
        ggets_plintf(" GLOBAL: sign =%d \n", sign);
        my_string = form_ode_do_fit_get_next("{}");
        strcpy(condition, my_string);
        ggets_plintf(" condition = %s \n", condition);
        my_string = form_ode_do_fit_get_next("\n");
        strcpy(formula, my_string);
        ggets_plintf(" events=%s \n", formula);
        if (flags_add_global(condition, sign, formula)) {
            printf("Bad global !! \n");
            exit(0);
        }
        if (ConvertStyle) {
            fprintf(convertf, "global %d {%s} %s\n", sign, condition, formula);
        }
        break;
    case 'p':
        ggets_plintf("Parameters:\n");
        if (ConvertStyle) {
            fprintf(convertf, "par ");
        }

        form_ode_advance_past_first_word(&ptr);

        while ((my_string = form_ode_get_next2(&ptr)) != NULL) {
            form_ode_take_apart(my_string, &value, name);
            free(my_string);
            default_val[NUPAR] = value;
            strcpy(upar_names[NUPAR++], name);
            if (ConvertStyle) {
                fprintf(convertf, "%s=%g  ", name, value);
            }
            ggets_plintf("|%s|=%f ", name, value);
            if (parserslow_add_con(name, value)) {
                ggets_plintf("ERROR at line %d\n", NLINES);
                exit(0);
            }
        }
        if (ConvertStyle) {
            fprintf(convertf, "\n");
        }
        ggets_plintf("\n");
        break;
    case 'c':
        my_string = form_ode_do_fit_get_next(" \n");
        strcpy(options, my_string);
        ggets_plintf(" Loading new options file:<%s>\n", my_string);
        if (ConvertStyle) {
            fprintf(convertf, "option %s\n", options);
        }
        break;
    case 'f':
        iflg = 0;
        ggets_plintf("\nFixed variables:\n");
        goto vrs;
    case 'm': /* Markov variable  */
        my_string = form_ode_do_fit_get_next(" ");
        strcpy(name, my_string);
        my_string = form_ode_do_fit_get_next(" ");
        value = atof(my_string);
        my_string = form_ode_do_fit_get_next(" \n");
        nstates = atoi(my_string);
        parserslow_add_var(name, value);
        strcpy(uvar_names[IN_VARS + NMarkov], name);
        last_ic[IN_VARS + NMarkov] = value;
        default_ic[IN_VARS + NMarkov] = value;
        ggets_plintf(" Markov variable %s=%f has %d states \n", name, value,
                     nstates);
        if (OldStyle) {
            markov_add(nstates, name);
        }
        if (ConvertStyle) {
            fprintf(convertf, "%s(0)=%g\n", name, value);
        }
        break;
    case 'r': /* state table for Markov variables  */
        my_string = form_ode_do_fit_get_next("\n");
        strcpy(name, my_string);
        nlin = NLINES;
        index = markov_old_build(fptr, name);
        nn = (int32)strlen(save_eqn[nlin]);
        if ((ode_names[IN_VARS + index] = xmalloc((usize)nn + 10)) == NULL) {
            exit(0);
        }
        strcpy(formula, save_eqn[nlin]);
        sprintf(ode_names[IN_VARS + index], "{ %s ... }", formula);
        break;
    case 'v':
        iflg = 1;
        ggets_plintf("\nVariables:\n");
        if (ConvertStyle) {
            fprintf(convertf, "init ");
        }
    vrs:
        if (NMarkov > 0 && OldStyle) {
            printf(" Error at line %d \n Must declare Markov variables after "
                   "fixed "
                   "and regular variables\n",
                   NLINES);
            exit(0);
        }
        form_ode_advance_past_first_word(&ptr);
        while ((my_string = form_ode_get_next2(&ptr)) != NULL) {
            if ((IN_VARS > NEQ) || (IN_VARS == MAX_ODE)) {
                ggets_plintf(" too many variables at line %d\n", NLINES);
                exit(0);
            }
            form_ode_take_apart(my_string, &value, name);
            free(my_string);
            if (parserslow_add_var(name, value)) {
                ggets_plintf("ERROR at line %d\n", NLINES);
                exit(0);
            }
            if (iflg) {
                strcpy(uvar_names[IN_VARS], name);
                last_ic[IN_VARS] = value;
                default_ic[IN_VARS] = value;
                IN_VARS++;
                if (ConvertStyle) {
                    fprintf(convertf, "%s=%g  ", name, value);
                }
            } else {
                if (ConvertStyle) {
                    strcpy(fixname[FIX_VAR], name);
                }
                FIX_VAR++;
            }
            ggets_plintf("|%s| ", name);
        }
        ggets_plintf(" \n");
        if (iflg && ConvertStyle) {
            fprintf(convertf, "\n");
        }
        break;
    case 'b':
        my_string = form_ode_do_fit_get_next("\n");
        my_bc[BVP_N].com = xmalloc(200*sizeof(*(my_bc[BVP_N].com)));
        my_bc[BVP_N].string = xmalloc(256);
        my_bc[BVP_N].name = xmalloc(10);
        strcpy(my_bc[BVP_N].string, my_string);
        strcpy(my_bc[BVP_N].name, "0=");
        if (ConvertStyle) {
            fprintf(convertf, "bndry %s\n", my_bc[BVP_N].string);
        }

        ggets_plintf("|%s| |%s| \n", my_bc[BVP_N].name, my_bc[BVP_N].string);
        BVP_N++;
        break;
    case 'k':
        if (ConvertStyle) {
            printf(" Warning  kernel declaration cannot be converted \n");
        }
        my_string = form_ode_do_fit_get_next(" ");
        strcpy(name, my_string);
        my_string = form_ode_do_fit_get_next(" ");
        value = atof(my_string);
        my_string = form_ode_do_fit_get_next("$");
        strcpy(formula, my_string);
        ggets_plintf("Kernel mu=%f %s = %s \n", value, name, formula);
        if (parserslow_add_kernel(name, value, formula)) {
            printf("ERROR at line %d\n", NLINES);
            exit(0);
        }
        break;
    case 't':
        if (NTable >= MAX_TAB) {
            if (ERROUT) {
                printf("too many tables !!\n");
            }
            exit(0);
        }
        my_string = form_ode_do_fit_get_next(" ");
        strcpy(name, my_string);
        my_string = form_ode_do_fit_get_next(" \n");
        if (my_string[0] == '%') {
            printf(" Function form of table....\n");
            my_string = form_ode_do_fit_get_next(" ");
            nn = atoi(my_string);
            my_string = form_ode_do_fit_get_next(" ");
            xlo = atof(my_string);
            my_string = form_ode_do_fit_get_next(" ");
            xhi = atof(my_string);
            my_string = form_ode_do_fit_get_next("\n");
            strcpy(formula, my_string);
            printf(" %s has %d pts from %f to %f = %s\n", name, nn, xlo, xhi,
                   formula);
            parserslow_add_table_name(NTable, name);

            if (parserslow_add_form_table(NTable, nn, xlo, xhi, formula)) {
                ggets_plintf("ERROR at line %d\n", NLINES);
                exit(0);
            }

            if (ConvertStyle) {
                fprintf(convertf, "table %s %% %d %g %g %s\n", name, nn, xlo,
                        xhi, formula);
            }
            NTable++;
            printf(" NTable = %d \n", NTable);

        } else if (my_string[0] == '@') {
            ggets_plintf(" Two-dimensional array: \n ");
            my_string = form_ode_do_fit_get_next(" ");
            strcpy(formula, my_string);
            ggets_plintf(" %s = %s \n", name, formula);
        } else {
            strcpy(formula, my_string);
            ggets_plintf("Lookup table %s = %s \n", name, formula);
            parserslow_add_table_name(NTable, name);
            if (parserslow_add_file_table(NTable, formula)) {
                ggets_plintf("ERROR at line %d\n", NLINES);
                exit(0);
            }
            if (ConvertStyle) {
                fprintf(convertf, "table %s %s\n", name, formula);
            }
            NTable++;
        }
        break;

    case 'u':
        my_string = form_ode_do_fit_get_next(" ");
        strcpy(name, my_string);
        my_string = form_ode_do_fit_get_next(" ");
        narg = atoi(my_string);
        my_string = form_ode_do_fit_get_next("$");
        strcpy(formula, my_string);
        ggets_plintf("%s %d :\n", name, narg);
        if (ConvertStyle) {
            fprintf(convertf, "%s(", name);
            for (i = 0; i < narg; i++) {
                fprintf(convertf, "arg%d", i + 1);
                if (i < (narg - 1)) {
                    fprintf(convertf, ",");
                }
            }
            fprintf(convertf, ")=%s", formula);
        }
        if (parserslow_add_ufun(name, formula, narg)) {
            printf("ERROR at line %d\n", NLINES);
            exit(0);
        }

        ggets_plintf("user %s = %s\n", name, formula);
        break;
    case 'i':
        VFlag = 1;
        __attribute__((fallthrough));
    case 'o':
        if (NODE >= (NEQ + FIX_VAR - NMarkov)) {
            done = 0;
            break;
        }
        my_string = form_ode_do_fit_get_next("\n");
        strcpy(formula, my_string);
        nn = (int32)strlen(formula) + 1;
        if ((my_ode[NODE] = xmalloc(MAXEXPLEN*sizeof(int32))) == NULL) {
            printf("Out of memory at line %d\n", NLINES);
            exit(0);
        }

        if (NODE < IN_VARS) {
            if ((ode_names[NODE] = xmalloc((usize)nn + 5)) == NULL) {
                ggets_plintf("Out of memory at line %d\n", NLINES);
                exit(0);
            }
            strcpy(ode_names[NODE], formula);
            if (ConvertStyle) {
                if (VFlag) {
                    fprintf(convertf, "volt %s=%s\n", uvar_names[NODE],
                            formula);
                } else {
                    fprintf(convertf, "%s'=%s\n", uvar_names[NODE], formula);
                }
            }
            form_ode_find_ker(formula, &alt);

            EqType[NODE] = VFlag;
            VFlag = 0;
        }
        if (NODE >= IN_VARS && NODE < (IN_VARS + FIX_VAR)) {
            if (ConvertStyle) {
                fprintf(convertf, "%s=%s\n", fixname[NODE - IN_VARS], formula);
            }
            form_ode_find_ker(formula, &alt);
        }

        if (NODE >= (IN_VARS + FIX_VAR)) {
            i = NODE - (IN_VARS + FIX_VAR);
            if ((ode_names[NODE - FIX_VAR + NMarkov] = xmalloc((usize)nn)) ==
                NULL) {
                ggets_plintf("Out of memory at line %d\n", NLINES);
                exit(0);
            }
            strcpy(ode_names[NODE - FIX_VAR + NMarkov], formula);
            if (ConvertStyle) {
                if (i < Naux) {
                    fprintf(convertf, "aux %s=%s\n", aux_names[i], formula);
                } else {
                    fprintf(convertf, "aux aux%d=%s\n", i + 1, formula);
                }
            }
        }
        ggets_plintf("RHS(%d)=%s\n", NODE, formula);
        if (parserslow_add_expr(formula, my_ode[NODE], &leng[NODE])) {
            printf("ERROR at line %d\n", NLINES);
            exit(0);
        }
        NODE++;
        break;

    case 'a': /* name auxiliary variables */
        ggets_plintf("Auxiliary variables:\n");
        while ((my_string = form_ode_do_fit_get_next(" ,\n")) != NULL) {
            strcpy(aux_names[Naux], my_string);
            ggets_plintf("|%s| ", aux_names[Naux]);
            Naux++;
        }
        ggets_plintf("\n");
        break;

    default:
        if (ConvertStyle) {
            my_string = form_ode_do_fit_get_next("\n");
            fprintf(convertf, "%s %s\n", command, my_string);
        }
        break;
    }

    return done;
}

/* ram: do I need to strip the name of any whitespace? */
void
form_ode_take_apart(char *bob, double *value, char *name) {
    int32 k;
    int32 l;
    char number[40];
    l = (int32)strlen(bob);
    k = (int32)strcspn(bob, "=");
    if (k == l) {
        *value = 0.0;
        form_ode_strcpy_trim(name, bob);
    } else {
        form_ode_strncpy_trim(name, bob, k);
        name[k] = '\0';
        for (int32 i = k + 1; i < l; i++) {
            number[i - k - 1] = bob[i];
        }
        number[l - k - 1] = '\0';
        *value = atof(number);
    }
    return;
}

char *
form_ode_get_first(char *string, char *src) {
    char *ptr;
    ptr = strtok(string, src);
    return ptr;
}

char *
form_ode_do_fit_get_next(char *src) {
    char *ptr;
    ptr = strtok(NULL, src);
    return ptr;
}

void
form_ode_find_ker(char *string, int32 *alt) {
    /* this extracts the integral operators from the string */
    char new[MAXEXPLEN];
    char form[MAXEXPLEN];
    char num[MAXEXPLEN];
    double mu = 0.0;
    int32 fflag = 0;
    int32 in = 0;
    int32 i = 0;
    int32 ifr = 0;
    int32 inum = 0;
    int32 n = (int32)strlen(string);
    char name[20];
    char ch;
    *alt = 0;
    while (i < n) {
        ch = string[i];
        if (ch == '[') {
            in = in - 3;
            inum = 0;
            i++;
            while ((ch = string[i]) != ']') {
                num[inum] = ch;
                inum++;
                i++;
            }
            mu = atof(num);
            fflag = 1;
            *alt = 1;
            ifr = 0;
            i += 2;
            continue;
        }
        if (ch == '{') {
            in = in - 3;
            fflag = 1;
            ifr = 0;
            *alt = 1;
            i++;
            continue;
        }
        if (ch == '}') {
            form[ifr] = 0;
            snprintf(name, sizeof(name), "K##%d", NKernel);
            ggets_plintf("Kernel mu=%f %s = %s \n", mu, name, form);
            if (parserslow_add_kernel(name, mu, form)) {
                exit(0);
            }
            for (usize j = 0; j < strlen(name); j++) {
                new[in] = name[j];
                in++;
            }
            mu = 0.0;
            ifr = 0;
            fflag = 0;
            i++;
            continue;
        }
        if (fflag) {
            form[ifr] = ch;
            ifr++;
        } else {
            new[in] = ch;
            in++;
        }
        i++;
    }
    new[in] = 0;
    strcpy(string, new);
    return;
}

void
form_ode_clrscr(void) {
    system("clear");
    return;
}

/***   remove this for full PP   ***/

/*   This is the new improved parser for input files.
     It is much more natural.  The format is as follows:

 * # comments
 * par  name=val, ....
 * init name=val,...
 * number name=value, ...
 * wiener name,..
 * table name ...
 * markov name #states (replaces m r)
 * { }  ..... { }
 * .
 * .
 * { }  ..... { }
 * options filename
 * aux name = expression
 * bndry ....
 * global ...
 * special name=conv(....)
 * special name=sparse(...)
 *
 * u' = expression    \
 *                     ----   Differential equations (replaces o v)
 * du/dt = expression /
 *
 * u(t+1) = expression >--- Difference equation   (replace o v)
 *
 * u(t) = expression with int32{} or int32[] <-- volterra eq (replaces i v)
 *
 * f(x,y,...) = expression >----   function (replaces u)
 *
 * u = expression>---  fixed  (replaces f o)
 *
 * u(0) = value >---  initial data (replaces v, init is also OK ) */

/*
 * XPP INTERNALS DEMAND THE FOLLOWING ORDER CONVENTION:

 * external names :  ODES Markov AUXILLIARY (uvar_names)
 * internal names :  ODES FIXED Markov  (variables)
 * internal formula: ODES FIXED AUXILLIARY (my_ode)
 * external formula: odes markov auxilliary (ode_names)

 * NODE = #ode variables
 * NMarkov = # Markov variables
 * NAux = # named auxiliary variables
 * NEQ = ode+naux   --> plotted quantities

 * my_ode[] <---  formulas
 * ode_names[] <---- "rhs"
 * uvar_names[] <----\
 * aux_names[]  <----/ external names

 * New parser reads in each line storing it in the VarInfo structure
 * if it is a markov (the only truly multiline command) then it
 * ** immediately ** reads in the markov stuff

 * It makes free use of "compiler"  in the old parser by
 * sending it new strings

 * On the first pass it does nothing except markov stuff
 * On the second pass it imitates an ode file doing things in the
 * "correct" order

 *  Only functions have changed syntax ...  */

int32
form_ode_if_include_file(char *old, char *nf) {
    int32 i = 0;
    int32 n = (int32)strlen(old);
    char c;

    if (strncmp(old, "#include", 8) == 0) {
        while (true) {
            c = old[i];
            if (c == ' ') {
                break;
            }
            i++;
            if (i == n) {
                return 0;
            }
        }
        for (int32 j = i + 1; j < n; j++) {
            nf[j - i - 1] = old[j];
        }
        nf[n - i - 1] = 0;
        ani_de_space(nf);
        return 1;
    }
    return 0;
}

int32
form_ode_if_end_include(char *old) {
    if (IN_INCLUDED_FILE > 0) {
        if (strncmp(old, "#done", 5) == 0) {
            return 1;
        }
        if (strncmp(old, "done", 4) == 0) {
            return 1;
        }
        /*Note that feof termination of an included file
         is also possible but that condition is checked
        elsewhere (currently near the bottom of do_new_parser)
        */
    }
    return 0;
}

int32
form_ode_do_new_parser(FILE *fp, char *first, int32 nnn) {
    VarInfo v;
    char **markovarrays = NULL;
    char *strings[256];
    int32 nstrings = 0;
    int32 ns;
    char **markovarrays2 = NULL;
    int32 done = 0;
    int32 start = 0;
    int32 i0;
    int32 i1;
    int32 i2;
    int32 istates;
    int32 jj1 = 0;
    int32 jj2 = 0;
    int32 jj;
    int32 notdone = 1;
    int32 jjsgn = 1;
    char name[20];
    int32 nstates = 0;
    /*char newfile[256];*/
    char newfile[XPP_MAX_NAME];
    FILE *fnew;
    /*int32 nlin;
     */
    char big[MAXEXPLEN];
    char old[MAXEXPLEN];
    char new[MAXEXPLEN];
    char *my_string;
    int32 is_array = 0;
    if (nnn == 0) {
        form_ode_init_varinfo();
    }
    while (notdone) {
        nstrings = 0;
        if (start || nnn == 1) {
            form_ode_read_a_line(fp, old);

        } else {
            if (loadincludefile) {
                loadincludefile = 0; /*Only do this once*/
                for (int32 j = 0; j < NincludedFiles; j++) {
                    printf("Trying to open %d %s\n", NincludedFiles,
                           includefilename[j]);
                    fnew = fopen(includefilename[j], "r");
                    if (fnew == NULL) {
                        ggets_plintf("Can't open include file <%s>\n",
                                     includefilename[j]);
                        exit(-1);
                        /*continue;*/
                    }
                    ggets_plintf("Including %s \n", includefilename[j]);
                    IN_INCLUDED_FILE++;
                    form_ode_do_new_parser(fnew, includefilename[j], 1);
                    fclose(fnew);
                }
                /*continue;*/
            }

            strcpy(old, first); /* pass the first line ....  */
            start = 1;
        }
        if (IN_INCLUDED_FILE > 0) {
            if (form_ode_if_end_include(old) || feof(fp)) {
                ggets_plintf("Completed include of file %s\n", first);
                IN_INCLUDED_FILE--;
                return 1;
            }
        }
        if (form_ode_if_include_file(old, newfile)) {
            fnew = fopen(newfile, "r");
            if (fnew == NULL) {
                ggets_plintf("Cant open include file <%s>\n", newfile);
                continue;
            }
            ggets_plintf("Including %s...\n", newfile);
            IN_INCLUDED_FILE++;
            form_ode_do_new_parser(fnew, newfile, 1);
            fclose(fnew);
            if (IN_INCLUDED_FILE > 0) {
                if (feof(fp)) {
                    /*ggets_plintf("We are at end of file now
                     * %d\n",IN_INCLUDED_FILE);*/
                }
            } else {
                continue;
            }
        }

        /*    printf("calling search %s \n",old); */
        form_ode_search_array(old, new, &jj1, &jj2, &is_array);
        jj = jj1;
        jjsgn = 1;
        if (jj2 < jj1) {
            jjsgn = -1;
        }

        switch (is_array) {
        case 0: /*  not a for loop so */
        case 1:
            nstrings = 1;
            strings[0] = xmalloc(strlen(new) + 10);
            strcpy(strings[0], new);
            break;
        case 2: /*  a for loop, so we will ignore the first line */
            /* is_array=1; */
            while (true) {
                form_ode_read_a_line(fp, old);
                if (old[0] == '%') {
                    break;
                }
                strings[nstrings] = xmalloc(strlen(old) + 10);
                strcpy(strings[nstrings], old);
                nstrings++;
                if (nstrings > 255) {
                    break;
                }
            }

            break;
        default:
            fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
            exit(EXIT_FAILURE);
        }

        while (true) {
            for (ns = 0; ns < nstrings; ns++) {
                strcpy(new, strings[ns]);
                form_ode_subsk(new, big, jj, is_array);

                done = form_ode_parse_a_string(big, &v);

                if (done == -1) {
                    ggets_plintf(" Error in parsing %s \n", big);
                    return -1;
                }
                if (done == 1) {
                    if (v.type == COMMAND) {
                        strupr(v.lhs);
                    }
                    if (v.type == COMMAND && v.lhs[0] == 'G' &&
                        v.lhs[1] == 'R') {
                        my_string = form_ode_get_first(v.rhs, " ");
                        strcpy(name, my_string);
                        my_string = form_ode_do_fit_get_next(" \n");
                        if (my_string == NULL) {
                            nstates = 0;
                        } else {
                            nstates = atoi(my_string);
                        }
                        if (nstates < 1) {
                            ggets_plintf(
                                "Group %s  must have at least 1 part \n", name);
                            return -1;
                        }
                        ggets_plintf("Group %s has %d parts\n", name, nstates);
                        for (istates = 0; istates < nstates; istates++) {
                            form_ode_read_a_line(fp, old);
                            ggets_plintf("part %d is %s \n", istates, old);
                        }

                        v.type = GROUP;
                    }
                    /* check for Markov to get rid of extra lines */

                    if (v.type == COMMAND && v.lhs[0] == 'M' &&
                        v.lhs[1] == 'A') {
                        my_string = form_ode_get_first(v.rhs, " ");
                        strcpy(name, my_string);
                        my_string = form_ode_do_fit_get_next(" \n");
                        if (my_string == NULL) {
                            nstates = 0;
                        } else {
                            nstates = atoi(my_string);
                        }
                        if (nstates < 2) {
                            ggets_plintf(
                                "Markov variable %s  must have at least 2 "
                                "states \n",
                                name);
                            return -1;
                        }
                        /*nlin=NLINES;
                         */
                        markov_add(nstates, name);
                        if (jj ==
                            jj1) { /* test to see if this is the first one */
                            markovarrays =
                                xmalloc((usize)nstates*sizeof(char *));
                            markovarrays2 =
                                xmalloc((usize)nstates*sizeof(char *));

                            for (istates = 0; istates < nstates; istates++) {
                                markovarrays[istates] = xmalloc(MAXEXPLEN);
                                markovarrays2[istates] = xmalloc(MAXEXPLEN);
                                /* fgets(markovarrays[istates],MAXEXPLEN,fp); */

                                if (is_array == 2) {
                                    strcpy(markovarrays[istates],
                                           strings[ns + 1 + istates]);

                                } else {
                                    form_ode_read_a_line(fp,
                                                         markovarrays[istates]);
                                }
                            }
                        }

                        /*  now we clean up these arrays */
                        for (istates = 0; istates < nstates; istates++) {
                            form_ode_subsk(markovarrays[istates],
                                           markovarrays2[istates], jj,
                                           is_array);
                        }

                        build_markov(markovarrays2, name);
                        v.type = MARKOV_VAR;
                        strcpy(v.lhs, name);
                        /* strcpy(v.rhs,save_eqn[nlin]); */
                        strcpy(v.rhs, "...many states..");
                    }

                    /* take care of special form for SOLVE-VARIABLE */
                    if (v.type == COMMAND && v.lhs[0] == 'S' &&
                        v.lhs[1] == 'O') {
                        if (form_ode_find_char(v.rhs, "=", 0, &i1) < 0) {
                            strcpy(v.lhs, v.rhs);
                            strcpy(v.rhs, "0");
                        } else {
                            form_ode_strpiece(v.lhs, v.rhs, 0, i1 - 1);
                            strcpy(big, v.rhs);
                            form_ode_strpiece(v.rhs, big, i1 + 1,
                                              (int32)strlen(big));
                        }
                        v.type = SOL_VAR;
                    }

                    /* take care of special form for auxiliary */
                    if (v.type == COMMAND && v.lhs[0] == 'A' &&
                        v.lhs[1] == 'U') {
                        form_ode_find_char(v.rhs, "=", 0, &i1);
                        form_ode_strpiece(v.lhs, v.rhs, 0, i1 - 1);
                        strcpy(big, v.rhs);
                        form_ode_strpiece(v.rhs, big, i1 + 1,
                                          (int32)strlen(big));
                        v.type = AUX_VAR;
                    }

                    /* take care of special form for vector */
                    if (v.type == COMMAND && v.lhs[0] == 'V' &&
                        v.lhs[1] == 'E' && v.lhs[5] == 'R') {
                        form_ode_find_char(v.rhs, "=", 0, &i1);
                        form_ode_strpiece(v.lhs, v.rhs, 0, i1 - 1);
                        strcpy(big, v.rhs);
                        form_ode_strpiece(v.rhs, big, i1 + 1,
                                          (int32)strlen(big));
                        v.type = VECTOR;
                    }
                    /* take care of special form for special */
                    if (v.type == COMMAND && v.lhs[0] == 'S' &&
                        v.lhs[1] == 'P' && v.lhs[5] == 'A') {
                        form_ode_find_char(v.rhs, "=", 0, &i1);
                        form_ode_strpiece(v.lhs, v.rhs, 0, i1 - 1);
                        strcpy(big, v.rhs);
                        form_ode_strpiece(v.rhs, big, i1 + 1,
                                          (int32)strlen(big));
                        v.type = SPEC_FUN;
                    }

                    /*   import-export to external C program   */
                    if (v.type == COMMAND && v.lhs[0] == 'E' &&
                        v.lhs[1] == 'X') {
                        v.type = EXPORT;
                        form_ode_find_char(v.rhs, "}", 0, &i1);
                        form_ode_strpiece(v.lhs, v.rhs, 0, i1);
                        strcpy(big, v.rhs);
                        form_ode_strpiece(v.rhs, big, i1 + 1,
                                          (int32)strlen(big));
                    }

                    /*  ONLY save options  */

                    if (v.type == COMMAND && v.lhs[0] == 'O' &&
                        v.lhs[1] == 'N') {
                        form_ode_break_up_list(v.rhs);
                        v.type = ONLY;
                    }

                    /*  forced integral equation form */
                    if (v.type == COMMAND && v.lhs[0] == 'V') {
                        form_ode_find_char(v.rhs, "=", 0, &i1);
                        form_ode_strpiece(v.lhs, v.rhs, 0, i1 - 1);
                        strcpy(big, v.rhs);
                        form_ode_strpiece(v.rhs, big, i1 + 1,
                                          (int32)strlen(big));
                        v.type = VEQ;
                    }
                    /* take care of tables   */

                    if (v.type == COMMAND && v.lhs[0] == 'T' &&
                        v.lhs[1] == 'A') {
                        i0 = 0;
                        form_ode_next_nonspace(v.rhs, i0, &i1);
                        i0 = i1;
                        i2 = form_ode_find_char(v.rhs, " ", i0, &i1);
                        if (i2 != 0) {
                            printf(" Illegal definition of table %s \n", v.rhs);
                            exit(0);
                        }
                        form_ode_strpiece(v.lhs, v.rhs, i0, i1 - 1);
                        strcpy(big, v.rhs);
                        form_ode_strpiece(v.rhs, big, i1 + 1,
                                          (int32)strlen(big));
                        v.type = TABLE;
                    }

                    /* printf("v.lhs=%s v.rhs=%s v.type=%d
                     * v.args=%s\n",v.lhs,v.rhs,v.type,v.args);
                     */
                    form_ode_add_varinfo(v.type, v.lhs, v.rhs, v.nargs, v.args);
                    switch (v.type) {
                    case ODE:
                    case MAP:
                        break;
                    case FIXED:
                        break;
                    case VEQ:
                        break;
                    case MARKOV_VAR:
                        break;
                    case DERIVE_PAR:
                    case PAR_AM:
                        break;
                    case SOL_VAR:
                        break;
                    default:
                        break;
                    }
                }
            } /* end loop for the strings */
            /*     if(nstrings>0){
              for(i=0;i<nstrings;i++)
                 free(strings[i]);
              nstrings=0;

              } */
            if (done == 2) {
                notdone = 0;
            }
            if (feof(fp)) {
                /*if (IN_INCLUDED_FILE>0)
                {
                        ggets_plintf("End of include file reached NOW \n");
                        IN_INCLUDED_FILE--;
                        return 1;
                }*/
                notdone = 0;
            }

            if (jj == jj2) {
                break;
            }

            jj += jjsgn;
        }

        if (v.type == COMMAND && v.lhs[0] == 'M' && v.lhs[1] == 'A') {
            for (istates = 0; istates < nstates; istates++) {
                free(markovarrays[istates]);
                free(markovarrays2[istates]);
            }
            free(markovarrays);
            free(markovarrays2);
        }
    }
    for (ns = 0; ns < nstrings; ns++) {
        free(strings[ns]);
    }
    {
        /* form ode compile em */
        /* Now we try to keep track of markov, fixed, etc as well as their names
         */
        VarInfo *v2;
        char vnames[MAX_ODE1][MAXVNAM];
        char fnames[MAX_ODE1][MAXVNAM];
        char anames[MAX_ODE1][MAXVNAM];
        char mnames[MAX_ODE1][MAXVNAM];
        double z, xlo, xhi;
        char tmp[50];
        char big2[2*MAXEXPLEN + 10];
        char formula[MAXEXPLEN];
        char *junk, *ptr;
        int32 nmark = 0, nfix = 0, naux = 0, nvar = 0;
        int32 nn, alt, in, ntab = 0, nufun = 0;
        int32 in1, in2, iflag;
        int32 fon;
        FILE *fp2 = NULL;

        v2 = my_varinfo;
        /* On this first pass through, all the variable names
           are kept as well as fixed declarations, boundary conds,
           and parameters, functions and tables.  Once this pass is
           completed all the names will be known to the compiler.
        */
        while (true) {
            if (v2->type == COMMAND && v2->lhs[0] == 'P') {
                snprintf(big2, sizeof(big2), "par %s \n", v2->rhs);
                form_ode_compiler(big2, fp2);
            }
            if (v2->type == COMMAND && v2->lhs[0] == 'W') {
                snprintf(big2, sizeof(big2), "wie %s \n", v2->rhs);
                form_ode_compiler(big2, fp2);
            }
            if (v2->type == COMMAND && v2->lhs[0] == 'N') {
                snprintf(big2, sizeof(big2), "num %s \n", v2->rhs);
                form_ode_compiler(big2, fp2);
            }
            if (v2->type == COMMAND && v2->lhs[0] == 'O') {
                snprintf(big2, sizeof(big2), "c %s \n", v2->rhs);
                form_ode_compiler(big2, fp2);
            }
            if (v2->type == COMMAND && v2->lhs[0] == 'S' && v2->lhs[1] == 'E') {
                snprintf(big2, sizeof(big2), "x %s\n", v2->rhs);
                form_ode_compiler(big2, fp2);
            }

            if (v2->type == COMMAND && v2->lhs[0] == 'B') {
                snprintf(big2, sizeof(big2), "b %s \n", v2->rhs);
                form_ode_compiler(big2, fp2);
            }
            if (v2->type == COMMAND && v2->lhs[0] == 'G') {
                snprintf(big2, sizeof(big2), "g %s \n", v2->rhs);
                form_ode_compiler(big2, fp2);
            }
            if (v2->type == MAP || v2->type == ODE || v2->type == VEQ) {
                convert(v2->lhs, tmp);
                if (form_ode_find_the_name(vnames, nvar, tmp) < 0) {
                    strcpy(vnames[nvar], tmp);
                    nvar++;
                } else {
                    ggets_plintf(" %s is a duplicate name \n", tmp);
                    exit(0);
                }
            }

            if (v2->type == MARKOV_VAR) {
                convert(v2->lhs, tmp);
                if (form_ode_find_the_name(mnames, nmark, tmp) < 0) {
                    strcpy(mnames[nmark], tmp);
                    nmark++;
                }
            }
            if (v2->type == EXPORT) {
                extra_add_export_list(v2->lhs, v2->rhs);
            }
            if (v2->type == VECTOR) {
                simplenet_add_vectorizer_name(v2->lhs, v2->rhs);
            }
            if (v2->type == SPEC_FUN) {
                simplenet_add_special_name(v2->lhs, v2->rhs);
            }
            if (v2->type == SOL_VAR) {
                if (dae_fun_add_svar(v2->lhs, v2->rhs) == 1) {
                    exit(0);
                }
            }

            if (v2->type == AUX_VAR) {
                convert(v2->lhs, tmp);
                strcpy(anames[naux], tmp);
                naux++;
                ggets_plintf("%s = %s \n", anames[naux - 1], v2->rhs);
            }
            if (v2->type == DERIVE_PAR) {
                if (derived_add(v2->lhs, v2->rhs) == 1) {
                    exit(0);
                }
            }
            if (v2->type == FIXED) {
                fixinfo[nfix].name = xmalloc(strlen(v2->lhs) + 2);
                fixinfo[nfix].value = xmalloc(strlen(v2->rhs) + 2);
                strcpy(fixinfo[nfix].name, v2->lhs);
                strcpy(fixinfo[nfix].value, v2->rhs);
                convert(v2->lhs, tmp);
                strcpy(fnames[nfix], tmp);
                nfix++;
                ggets_plintf("%s = %s \n", fnames[nfix - 1], v2->rhs);
            }

            if (v2->type == TABLE) {
                convert(v2->lhs, tmp);
                if (parserslow_add_table_name(ntab, tmp) == 1) {
                    printf(" %s is duplicate name \n", tmp);
                    exit(0);
                }
                ggets_plintf("added name %d\n", ntab);
                ntab++;
            }

            if (v2->type == FUNCTION) {
                convert(v2->lhs, tmp);
                if (parserslow_add_ufun_name(tmp, nufun, v2->nargs) == 1) {
                    printf("Duplicate name or too many functions for %s \n",
                           tmp);
                    exit(0);
                }

                nufun++;
            }

            if (v2->next == NULL) {
                break;
            }
            v2 = v2->next;
        }

        /* now we add all the names of the variables and the
           fixed stuff
        */
        for (int32 i = 0; i < nvar; i++) {
            if (parserslow_add_var(vnames[i], 0.0)) {
                printf(" Duplicate name %s \n", vnames[i]);
                exit(0);
            }
            strcpy(uvar_names[i], vnames[i]);
            last_ic[i] = 0.0;
            default_ic[i] = 0.0;
        }
        for (int32 i = 0; i < nfix; i++) {
            if (parserslow_add_var(fnames[i], 0.0)) {
                printf(" Duplicate name %s \n", fnames[i]);
                exit(0);
            }
        }
        for (int32 i = 0; i < nmark; i++) {
            if (parserslow_add_var(mnames[i], 0.0)) {
                printf(" Duplicate name %s \n", mnames[i]);
                exit(0);
            }
            strcpy(uvar_names[i + nvar], mnames[i]);
            last_ic[i + nvar] = 0.0;
            default_ic[i + nvar] = 0.0;
        }
        for (int32 i = 0; i < naux; i++) {
            strcpy(aux_names[i], anames[i]);
        }
        dae_fun_add_svar_names();

        /* NODE = nvars ; Naux = naux ; NEQ = NODE+NMarkov+Naux ; FIX_VAR =
         * nfix; */

        IN_VARS = nvar;
        Naux = naux;
        NEQ = nvar + NMarkov + Naux;
        FIX_VAR = nfix;
        NTable = ntab;
        NFUN = nufun;

        /* Reset all this stuff so we align the indices correctly */

        nvar = 0;
        naux = 0;
        ntab = 0;
        nufun = 0;
        nfix = 0;
        nmark = 0;

        v2 = my_varinfo;
        while (true) {
            if (v2->type == COMMAND && v2->lhs[0] == 'I') {
                snprintf(big2, sizeof(big2), "i %s \n", v2->rhs);
                ptr = big2;
                junk = form_ode_get_first(ptr, " ,");
                if (junk == NULL) {
                    /*No more tokens.  Should this throw an error?*/
                }
                form_ode_advance_past_first_word(&ptr);
                while ((my_string = form_ode_get_next2(&ptr)) != NULL) {
                    form_ode_take_apart(my_string, &z, name);
                    free(my_string);
                    convert(name, tmp);
                    in = form_ode_find_the_name(vnames, IN_VARS, tmp);
                    if (in >= 0) {
                        last_ic[in] = z;
                        default_ic[in] = z;
                        set_val(tmp, z);
                        ggets_plintf(" Initial %s(0)=%g\n", tmp, z);
                    } else {
                        in = form_ode_find_the_name(mnames, NMarkov, tmp);
                        if (in >= 0) {
                            last_ic[in + IN_VARS] = z;
                            default_ic[in + IN_VARS] = z;
                            set_val(tmp, z);
                            ggets_plintf(" Markov %s(0)=%g\n", tmp, z);
                        } else {
                            ggets_plintf(
                                "In initial value statement no variable %s \n",
                                tmp);
                            exit(0);
                        }
                    }
                } /* end take apart */
            } /* end  init  command    */
            if (v2->type == IC) {
                convert(v2->lhs, tmp);
                fon = form_ode_formula_or_number(v2->rhs, &z);

                if (fon == 1) {
                    if (v2->rhs[0] == '-' &&
                        (isdigit(v2->rhs[1]) || (v2->rhs[1] == '.'))) {
                        z = atof(v2->rhs);
                    }
                }

                in = form_ode_find_the_name(vnames, IN_VARS, tmp);
                if (in >= 0) {
                    last_ic[in] = z;
                    default_ic[in] = z;
                    set_val(tmp, z);
                    /* if(fon==1) */
                    strcpy(delay_string[in], v2->rhs);

                    ggets_plintf(" Initial %s(0)=%s\n", tmp, v2->rhs);
                } else {
                    in = form_ode_find_the_name(mnames, NMarkov, tmp);
                    if (in >= 0) {
                        last_ic[in + IN_VARS] = z;
                        default_ic[in + IN_VARS] = z;
                        set_val(tmp, z);
                        ggets_plintf(" Markov %s(0)=%g\n", tmp, z);
                    } else {
                        ggets_plintf(
                            "In initial value statement no variable %s \n",
                            tmp);
                        exit(0);
                    }
                }
            } /* end IC stuff  */

            /*   all that is left is the right-hand sides !!   */
            iflag = 0;
            switch (v2->type) {
            case VEQ:
                iflag = 1;
                __attribute__((fallthrough));
            case ODE:
                __attribute__((fallthrough));
            case MAP:
                EqType[nvar] = iflag;
                nn = (int32)strlen(v2->rhs) + 1;
                if ((ode_names[nvar] = xmalloc((usize)nn + 2)) == NULL ||
                    (my_ode[nvar] = xmalloc(MAXEXPLEN*sizeof(int32))) ==
                        NULL) {
                    ggets_plintf("could not allocate space for %s \n", v2->lhs);
                    exit(0);
                }

                strcpy(ode_names[nvar], v2->rhs);
                form_ode_find_ker(v2->rhs, &alt);
                /*       ode_names[nvar][nn-1]=0; */
                if (parserslow_add_expr(v2->rhs, my_ode[nvar], &leng[nvar])) {
                    printf("A\n");
                    ggets_plintf("ERROR compiling %s' \n", v2->lhs);
                    exit(0);
                }
                /* fpr_command(my_ode[nvar]); */
                if (v2->type == MAP) {
                    ggets_plintf("%s(t+1)=%s\n", v2->lhs, v2->rhs);
                    is_a_map = 1;
                }
                if (v2->type == VEQ) {
                    ggets_plintf("%s(t)=%s\n", v2->lhs, v2->rhs);
                }
                if (v2->type == ODE) {
                    ggets_plintf("%d:d%s/dt=%s\n", nvar, v2->lhs, v2->rhs);
                }
                nvar++;
                break;
            case FIXED:
                form_ode_find_ker(v2->rhs, &alt);
                if ((my_ode[nfix + IN_VARS] =
                         xmalloc(MAXEXPLEN*sizeof(int32))) == NULL ||
                    parserslow_add_expr(v2->rhs, my_ode[nfix + IN_VARS],
                                        &leng[IN_VARS + nfix]) != 0) {
                    ggets_plintf(" Error allocating or compiling %s\n",
                                 v2->lhs);
                    exit(0);
                }
                nfix++;
                ggets_plintf("%s=%s\n", v2->lhs, v2->rhs);
                break;
            case DAE:
                if (dae_fun_add_aeqn(v2->rhs) == 1) {
                    exit(0);
                }
                ggets_plintf(" DAE eqn: %s=0 \n", v2->rhs);
                break;

            case AUX_VAR:
                in1 = IN_VARS + NMarkov + naux;
                in2 = IN_VARS + FIX_VAR + naux;
                nn = (int32)strlen(v2->rhs) + 1;
                if ((ode_names[in1] = xmalloc((usize)nn + 2)) == NULL ||
                    (my_ode[in2] = xmalloc(MAXEXPLEN*sizeof(int32))) ==
                        NULL) {
                    ggets_plintf("could not allocate space for %s \n", v2->lhs);
                    exit(0);
                }

                strcpy(ode_names[in1], v2->rhs);
                /* ode_names[in1][nn]=0; */
                if (parserslow_add_expr(v2->rhs, my_ode[in2], &leng[in2])) {
                    printf("B\n");
                    ggets_plintf("ERROR compiling %s \n", v2->lhs);
                    exit(0);
                }
                naux++;
                ggets_plintf("%s=%s\n", v2->lhs, v2->rhs);
                break;
            case VECTOR:
                if (simplenet_add_vectorizer(v2->lhs, v2->rhs) == 0) {
                    ggets_plintf(" Illegal vector  %s \n", v2->rhs);
                    exit(0);
                }

                break;
            case SPEC_FUN:
                if (simplenet_add_spec_fun(v2->lhs, v2->rhs) == 0) {
                    ggets_plintf(" Illegal special function %s \n", v2->rhs);
                    exit(0);
                }
                break;
            case MARKOV_VAR:
                nn = (int32)strlen(v2->rhs) + 1;

                if ((ode_names[IN_VARS + nmark] = xmalloc((usize)nn + 2)) ==
                    NULL) {
                    ggets_plintf(" Out of memory for  %s \n", v2->lhs);
                    exit(0);
                }
                strncpy(ode_names[IN_VARS + nmark], v2->rhs, (usize)nn);
                ode_names[IN_VARS + nmark][nn] = 0;
                nmark++;
                ggets_plintf("%s: %s", v2->lhs, v2->rhs);
                break;
            case FUNCTION:
                if (parserslow_add_ufun_new(nufun, v2->nargs, v2->rhs,
                                            v2->args) != 0) {
                    ggets_plintf(" Function %s messed up \n", v2->lhs);
                    exit(0);
                }
                nufun++;
                ggets_plintf("%s(%s", v2->lhs, v2->args[0]);
                for (in = 1; in < v2->nargs; in++) {
                    ggets_plintf(",%s", v2->args[in]);
                }
                ggets_plintf(")=%s\n", v2->rhs);
                break;

            case TABLE:
                snprintf(big2, sizeof(big2), "t %s %s ", v2->lhs, v2->rhs);
                ptr = big2;
                junk = form_ode_get_first(ptr, " ,");
                my_string = form_ode_do_fit_get_next(" ");
                my_string = form_ode_do_fit_get_next(" \n");
                if (my_string[0] == '%') {
                    ggets_plintf(" Function form of table....\n");
                    my_string = form_ode_do_fit_get_next(" ");
                    nn = atoi(my_string);
                    my_string = form_ode_do_fit_get_next(" ");
                    xlo = atof(my_string);
                    my_string = form_ode_do_fit_get_next(" ");
                    xhi = atof(my_string);
                    my_string = form_ode_do_fit_get_next("\n");
                    strcpy(formula, my_string);
                    ggets_plintf(" %s has %d pts from %f to %f = %s\n", v2->lhs,
                                 nn, xlo, xhi, formula);
                    if (parserslow_add_form_table(ntab, nn, xlo, xhi,
                                                  formula)) {
                        ggets_plintf("ERROR computing %s\n", v2->lhs);
                        exit(0);
                    }
                    ntab++;
                } else if (my_string[0] == '@') {
                    ggets_plintf(" Two-dimensional array: \n ");
                    my_string = form_ode_do_fit_get_next(" ");
                    strcpy(formula, my_string);
                    ggets_plintf(" %s = %s \n", name, formula);
                } else {
                    strcpy(formula, my_string);
                    ggets_plintf("Lookup table %s = %s \n", v2->lhs, formula);

                    if (parserslow_add_file_table(ntab, formula)) {
                        ggets_plintf("ERROR computing %s", v2->lhs);
                        exit(0);
                    }
                    ntab++;
                }
                break;
            default:
                break;
            }

            if (v2->next == NULL) {
                break;
            }
            v2 = v2->next;
        }
        if (derived_compile() == 1) {
            exit(0);
        }
        if (dae_fun_compile_svars() == 1) {
            exit(0);
        }
        derived_evaluate();
        extra_do_export_list();
        ggets_plintf(" All formulas are valid!!\n");
        NODE = nvar + naux + nfix;
        ggets_plintf(" nvar=%d naux=%d nfix=%d nmark=%d NEQ=%d NODE=%d \n",
                     nvar, naux, nfix, nmark, NEQ, NODE);
    }

    {
        /* form ode free varinfo */
        VarInfo *v3;
        VarInfo *vnew;
        v3 = my_varinfo;
        while (v3->next != NULL) {
            v3 = v3->next;
        }
        while (v3->prev != NULL) {
            vnew = v3->prev;
            v3->next = NULL;
            v3->prev = NULL;
            free(v3);
            v3 = vnew;
        }
        form_ode_init_varinfo();
    }
    return 1;
}

void
form_ode_create_plot_list(void) {
    int32 j = 0;
    int32 k;
    if (N_only == 0) {
        return;
    }
    plotlist = xmalloc(sizeof(*plotlist)*(usize)(N_only + 1));
    for (int32 i = 0; i < N_only; i++) {
        browse_find_variable(onlylist[i], &k);
        if (k >= 0) {
            plotlist[j] = k;
            j++;
        }
        N_plist = j;
    }
    return;
}

void
form_ode_add_only(char *s) {
    if (strlen(s) < 1) {
        return;
    }
    if (N_only >= MAXONLY) {
        return;
    }
    onlylist[N_only] = xmalloc(11);
    strcpy(onlylist[N_only], s);

    N_only++;
    return;
}

void
form_ode_break_up_list(char *rhs) {
    int32 i = 0;
    int32 j = 0;
    int32 l = (int32)strlen(rhs);
    char s[20];
    char c;
    while (i < l) {
        c = rhs[i];
        if (c == ' ' || c == ',') {
            s[j] = 0;
            form_ode_add_only(s);
            j = 0;
        } else {
            s[j] = c;
            j++;
        }
        i++;
    }
    s[j] = 0;
    form_ode_add_only(s);
    return;
}

int32
form_ode_find_the_name(char list[MAX_ODE1][MAXVNAM], int32 n, char *name) {
    for (int32 i = 0; i < n; i++) {
        if (strcmp(list[i], name) == 0) {
            return i;
        }
    }
    return -1;
}

/* this code checks if the right-hand side for an initial
 * condition is a formula (for delays) or a number
 */
int32
form_ode_formula_or_number(char *expr, double *z) {
    char num[80];
    char form[80];
    int32 flag;
    int32 i = 0;
    int32 olderr = ERROUT;
    ERROUT = 0;
    *z = 0.0; /* initial it to 0 */
    convert(expr, form);
    flag = do_num(form, num, z, &i);
    if (i < (int32)strlen(form)) {
        flag = 1;
    }
    ERROUT = olderr;
    if (flag == 0) {
        return 0; /* 0 is a number */
    }
    return 1; /* 1 is a formula */
}

void
form_ode_strpiece(char *dest, char *src, int32 i0, int32 ie) {
    for (int32 i = i0; i <= ie; i++) {
        dest[i - i0] = src[i];
    }
    dest[ie - i0 + 1] = 0;
    return;
}

int32
form_ode_parse_a_string(char *s1, VarInfo *v) {
    int32 i0 = 0;
    int32 i1;
    int32 i2;
    int32 i3;
    char lhs[MAXEXPLEN], rhs[MAXEXPLEN], args[MAXARG][NAMLEN + 1];
    int32 type;
    int32 type2;
    int32 narg = 0;
    int32 n1 = (int32)strlen(s1) - 1;
    char s1old[MAXEXPLEN];
    char ch;
    if (s1[0] == '"') {
        form_ode_add_comment(s1);
        return 0;
    }
    if (s1[0] == '@') {
        /*    printf("internopts from parse string\n"); */
        load_eqn_stor_internopts(s1);
        return 0;
    }
    form_ode_remove_blanks(s1);

    strcpy(s1old, s1);
    strupr(s1);
    if (strlen(s1) < 1) {
        return 0;
    }
    if (s1[0] == '0' && s1[1] == '=') { /* ||(s1[1]==' '&&s1[2]=='='))) */
        type2 = DAE;
        snprintf(lhs, sizeof(lhs), "0=");
        form_ode_strpiece(rhs, s1, 2, n1);
        v->type = type2;
        strcpy(v->lhs, lhs);
        strcpy(v->rhs, rhs);
        goto good_type;
    }
    if (s1[0] == '#') {
        return 0;
    }

    type = form_ode_find_char(s1, " =/'(", i0, &i1);
    switch (type) {
    case 0:
        i0 = i1;
        ch = (char)form_ode_next_nonspace(s1, i0, &i2);
        switch (ch) {
        case '=':
            if (s1[0] == '!') {
                form_ode_strpiece(lhs, s1, 1, i1 - 1);
                form_ode_strpiece(rhs, s1, i2 + 1, n1);
                type2 = DERIVE_PAR;
                break;
            }
            form_ode_strpiece(lhs, s1, 0, i1 - 1);
            form_ode_strpiece(rhs, s1, i2 + 1, n1);
            type2 = FIXED;
            break;
        default:
            type2 = COMMAND;
            form_ode_strpiece(lhs, s1, 0, i1 - 1);
            form_ode_strpiece(rhs, s1old, i2, n1);
            break;
        }
        break;
    case 1:
        if (s1[0] == '!') {
            form_ode_strpiece(lhs, s1, 1, i1 - 1);
            form_ode_strpiece(rhs, s1, i1 + 1, n1);
            type2 = DERIVE_PAR;
            break;
        }

        type2 = FIXED;
        form_ode_strpiece(lhs, s1, 0, i1 - 1);
        form_ode_strpiece(rhs, s1, i1 + 1, n1);
        break;
    case 2:
        if (s1[0] != 'D') {
            return -1;
        }
        if (form_ode_extract(s1, &i2, i1)) {
            form_ode_strpiece(lhs, s1, 1, i1 - 1);
            form_ode_strpiece(rhs, s1, i2, n1);
            type2 = ODE;
        } else {
            return -1;
        }
        break;
    case 3:
        if (form_ode_extract(s1, &i2, i1)) {
            form_ode_strpiece(lhs, s1, 0, i1 - 1);
            form_ode_strpiece(rhs, s1, i2, n1);
            type2 = ODE;
        } else {
            return -1;
        }
        break;

    case 4:
        i0 = i1;
        if (form_ode_strparse(s1, "T+1)=", i0, &i2)) {
            type2 = MAP;
            is_a_map = 1;
            form_ode_strpiece(lhs, s1, 0, i1 - 1);
            form_ode_strpiece(rhs, s1, i2, n1);
            break;
        }
        if (form_ode_strparse(s1, "(0)=", i0 - 1, &i2)) {
            type2 = IC;
            form_ode_strpiece(lhs, s1, 0, i1 - 1);
            form_ode_strpiece(rhs, s1, i2, n1);
            break;
        }
        if (form_ode_strparse(s1, "T)=", i0, &i2)) {
            if (form_ode_strparse(s1, "INT{", 0, &i3) == 1 ||
                form_ode_strparse(s1, "INT[", 0, &i3) == 1) {
                type2 = VEQ;
                form_ode_strpiece(lhs, s1, 0, i1 - 1);
                form_ode_strpiece(rhs, s1, i2, n1);
                break;
            } else {
                type2 = FUNCTION;
                form_ode_extract_args(s1, i0 + 1, &i2, &narg, args);
                form_ode_strpiece(lhs, s1, 0, i0 - 1);
                form_ode_strpiece(rhs, s1, i2, n1);
                break;
            }
        }
        i0++;
        form_ode_extract_args(s1, i0, &i2, &narg, args);
        type2 = FUNCTION;
        form_ode_strpiece(lhs, s1, 0, i0 - 2);
        form_ode_strpiece(rhs, s1, i2, n1);
        break;
    default:
        return -1;
    }

good_type:
    v->type = type2;
    strcpy(v->lhs, lhs);
    strcpy(v->rhs, rhs);
    v->nargs = narg;
    for (int32 i = 0; i < narg; i++) {
        strcpy(v->args[i], args[i]);
    }

    if (lhs[0] == 'D' && type2 == COMMAND) {
        return 2;
    }
    return 1;
}

void
form_ode_init_varinfo(void) {
    my_varinfo = xmalloc(sizeof(VarInfo));
    my_varinfo->next = NULL;
    my_varinfo->prev = NULL;
    start_var_info = 0;
    return;
}

void
form_ode_add_varinfo(int32 type, char *lhs, char *rhs, int32 nargs,
                     char args[MAXARG][NAMLEN + 1]) {
    VarInfo *v;
    VarInfo *vnew;
    v = my_varinfo;
    if (start_var_info == 0) {
        v->type = type;
        v->nargs = nargs;
        strcpy(v->lhs, lhs);
        strcpy(v->rhs, rhs);
        for (int32 i = 0; i < nargs; i++) {
            strcpy(v->args[i], args[i]);
        }
        start_var_info = 1;
    } else {
        while (v->next != NULL) {
            v = (v->next);
        }
        v->next = xmalloc(sizeof(VarInfo));
        vnew = v->next;
        vnew->type = type;
        vnew->nargs = nargs;
        strcpy(vnew->lhs, lhs);
        strcpy(vnew->rhs, rhs);
        for (int32 i = 0; i < nargs; i++) {
            strcpy(vnew->args[i], args[i]);
        }
        vnew->next = NULL;
        vnew->prev = v;
    }
    return;
}

int32
form_ode_extract(/* name is char 1-i1  ie is start of rhs */
                 char *s1, int32 *ie, int32 i1) {
    int32 i = 0;
    int32 n = (int32)strlen(s1);

    i = i1;
    while (i < n) {
        if (s1[i] == '=') {
            *ie = i + 1;
            return 1;
        }
        i++;
    }
    return 0;
}

int32
form_ode_strparse(char *s1, char *s2, int32 i0, int32 *i1) {
    int32 i = i0;
    int32 n = (int32)strlen(s1);
    int32 m = (int32)strlen(s2);
    int32 j = 0;
    char ch;
    int32 start = 0;

    while (i < n) {
        ch = s1[i];
        if (start == 1) {
            if (ch == s2[j] || ch == ' ') {
                if (ch == s2[j]) {
                    j++;
                }
                i++;
                if (j == m) {
                    *i1 = i;
                    return 1;
                }
            } else {
                start = 0;
                j = 0;
            }
        } else /* just starting */
        {

            if (ch == s2[0]) {
                j++;
                i++;
                start = 1;
                if (j == m) { /* only one char */
                    *i1 = i;
                    return 1;
                }
            } else {
                i++;
            }
        }
    }
    return 0;
}

int32
form_ode_extract_args(char *s1, int32 i0, int32 *ie, int32 *narg,
                      char args[MAXARG][NAMLEN + 1]) {
    int32 i = i0;
    int32 n = (int32)strlen(s1);
    int32 type;
    int32 na = 0;
    int32 i1;
    while (i < n) {
        type = form_ode_find_char(s1, ",)", i, &i1);
        if (type == 0) {
            for (int32 k = i; k < i1; k++) {
                args[na][k - i] = s1[k];
            }
            args[na][i1 - i] = 0;
            na++;
            i = i1 + 1;
        }
        if (type == 1) {
            for (int32 k = i; k < i1; k++) {
                args[na][k - i] = s1[k];
            }
            args[na][i1 - i] = 0;
            na++;
            i = i1 + 1;
            form_ode_find_char(s1, "=", i, &i1);
            *ie = i1 + 1;
            *narg = na;
            return 1;
        }
    }
    return 0;
}

int32
form_ode_find_char(char *s1, char *s2, int32 i0, int32 *i1) {
    int32 m = (int32)strlen(s2);
    int32 n = (int32)strlen(s1);
    int32 i = i0;
    char ch;
    while (i < n) {
        ch = s1[i];
        for (int32 j = 0; j < m; j++) {
            if (ch == s2[j]) {
                *i1 = i;
                return j;
            }
        }
        i++;
    }
    return -1;
}

int32
form_ode_next_nonspace(char *s1, int32 i0, int32 *i1) {
    int32 i = i0;
    int32 n = (int32)strlen(s1);
    char ch;
    *i1 = n - 1;
    while (i < n) {
        ch = s1[i];
        if (ch != ' ') {
            *i1 = i;
            return (int32)ch;
        }
        i++;
    }
    return -1;
}

/* removes starting blanks from s  */
void
form_ode_remove_blanks(char *s1) {
    int32 i = 0;
    int32 n = (int32)strlen(s1);
    int32 l;
    char ch;
    while (i < n) {
        ch = s1[i];
        if (isspace(ch)) {
            i++;
        } else {
            break;
        }
    }
    if (i == n) {
        s1[0] = 0;
    } else {
        l = n - i;
        for (int32 j = 0; j < l; j++) {
            s1[j] = s1[j + i];
        }
        s1[l] = 0;
    }
    return;
}

void
form_ode_read_a_line(FILE *fp, char *s) {
    char temp[MAXEXPLEN];
    int32 n;
    int32 nn;
    int32 ok;
    int32 ihat = 0;
    s[0] = 0;
    ok = 1;

    while (ok) {
        ok = 0;
        fgets(temp, MAXEXPLEN, fp);

        nn = (int32)strlen(temp) + 1;
        if ((save_eqn[NLINES] = xmalloc((usize)nn)) == NULL) {
            exit(0);
        }
        strncpy(save_eqn[NLINES++], temp, (usize)nn);
        n = (int32)strlen(temp);
        for (int32 i = n - 1; i >= 0; i--) {
            if (temp[i] == '\\') {
                ok = 1;
                ihat = i;
            }
        }
        if (ok == 1) {
            temp[ihat] = 0;
        }
        strcat(s, temp);
    }
    n = (int32)strlen(s);
    /*  if((s[n-1]=='\n')||(s[n-1]=='\r'))
      {
        s[n-1]=0;
        n=strlen(s);
      }
    if((s[n-1]=='\n')||(s[n-1]=='\r'))
      {
        s[n-1]=0;
        n=strlen(s);
      }
    */

    if (s[n - 1] == '\n' || s[n - 1] == '\r') {
        s[n - 1] = ' ';
    }
    s[n] = ' ';
    s[n + 1] = 0;
    return;
}

int32
form_ode_search_array(char *old, char *new, int32 *i1, int32 *i2, int32 *flag) {
    int32 j;
    int32 l;
    int32 ileft;
    int32 iright;
    int32 n = (int32)strlen(old);
    char num1[20];
    char num2[20];
    char ch;
    char chp;
    ileft = n - 1;
    iright = -1;
    *i1 = 0;
    *i2 = 0;
    *flag = 0;
    strcpy(num1, "0");
    strcpy(num2, "0");
    if (old[0] == '#' || old[1] == '#') { /* check for comments */

        strcpy(new, old);

        return 1;
    }
    if (form_ode_check_if_ic(old) == 1) {
        integrate_extract_ic_data(old);
        strcpy(new, old);
        return 1;
    }
    for (int32 i = 0; i < n; i++) {
        ch = old[i];
        chp = old[i + 1];
        if (ch == '.' && chp == '.') {
            j = 0;
            *flag = 1;
            if (old[0] == '%') {
                *flag = 2; /*   FOR LOOP CONSTRUCTION  */
            }
            while (true) {
                ch = old[i + j];
                if (ch == '[') {
                    ileft = i + j;
                    l = 0;
                    for (int32 k = i + j + 1; k < i; k++) {
                        num1[l] = old[k];
                        l++;
                    }
                    num1[l] = 0;
                    break;
                }
                j--;
                if ((i + j) <= 0) {
                    *i1 = 0;
                    *i2 = 0;
                    strcpy(new, old);
                    ggets_plintf(
                        " Possible error in array %s -- ignoring it \n", old);
                    return 0; /* error in array  */
                }
            }
            j = 2;
            while (true) {
                ch = old[i + j];
                if (ch == ']') {
                    iright = i + j;
                    l = 0;
                    for (int32 k = 2; k < j; k++) {
                        num2[l] = old[i + k];
                        l++;
                    }
                    num2[l] = 0;
                    break;
                }
                j++;
                if ((i + j) >= n) {
                    *i1 = 0;
                    *i2 = 0;
                    strcpy(new, old);
                    ggets_plintf(
                        " Possible error in array  %s -- ignoring it \n", old);
                    return 0; /* error again   */
                }
            }
        }
    }
    /*  printf(" I have extracted [%s] and [%s] \n",num1,num2); */
    *i1 = atoi(num1);
    *i2 = atoi(num2);
    /* now we have the numbers and will get rid of the junk inbetween */
    l = 0;
    for (int32 i = 0; i <= ileft; i++) {
        new[l] = old[i];
        l++;
    }
    if (iright > 0) {
        new[l] = 'j';
        l++;
        for (int32 i = iright; i < n; i++) {
            new[l] = old[i];
            l++;
        }
    }
    new[l] = 0;
    return 1;
}

int32
form_ode_check_if_ic(char *big) {
    char c;
    int32 n = (int32)strlen(big);
    int32 j;
    j = 0;
    while (true) {
        c = big[j];
        if (c == ']') {
            if ((big[j + 1] == '(') && (big[j + 2] == '0') &&
                (big[j + 3] == ')')) {
                return 1;
            }
        }
        j++;
        if (j >= n) {
            break;
        }
    }
    return 0;
}

int32
form_ode_not_ker(/* returns 1 if string is not 'int32[' */
                 char *s, int32 i) {
    if (i < 3) {
        return 1;
    }
    if (s[i - 3] == 'i' && s[i - 2] == 'n' && s[i - 1] == 't') {
        return 0;
    }
    return 1;
}

int32
form_ode_is_comment(char *s) {
    int32 n = (int32)strlen(s);
    int32 i = 0;
    char c;
    while (true) {
        c = s[i];
        if (c == '#') {
            return 1;
        }
        if (isspace(c)) {
            i++;

            if (i >= n) {
                return 0;
            }
        } else {
            return 0;
        }
    }
}

void
form_ode_subsk(char *big, char *new, int32 k, int32 flag) {
    int32 i, n = (int32)strlen(big), inew, add, inum, m, isign, ok,
             multflag = 0;
    char ch;
    char chp;
    char num[20];
    inew = 0;
    i = 0;
    /*  if(big[0]=='#'){   */
    if (form_ode_is_comment(big)) {
        strcpy(new, big);
        return;
    }

    while (true) {
        ch = big[i];
        chp = big[i + 1];
        if (ch == '[' && chp != 'j' && form_ode_not_ker(big, i)) {
            ok = 1;
            add = 0;
            inum = 0;
            isign = 1;
            i++;
            while (ok) {
                ch = big[i];
                if (ch == ']') {
                    i++;
                    num[inum] = 0;
                    add = atoi(num);
                    snprintf(num, sizeof(num), "%d", add);
                    m = (int32)strlen(num);
                    for (int32 j = 0; j < m; j++) {
                        new[inew] = num[j];
                        inew++;
                    }
                    ok = 0;
                } else {
                    i++;
                    num[inum] = ch;
                    inum++;
                }
            }
        } else

            if (ch == '[' && chp == 'j') {
            if (flag == 0) {
                printf(" Illegal use of [j] at %s \n", big);
                exit(0);
            }
            add = 0;
            inum = 0;
            isign = 1;
            i += 2;
            ok = 1;
            while (ok) {
                if (i >= n) {
                    new[inew] = 0;
                    ggets_plintf(
                        "Error in %s The expression does not terminate. "
                        "Perhaps a ] "
                        "is missing.\n",
                        big);
                    exit(0);
                }
                ch = big[i];
                switch (ch) {
                case '+':
                    isign = 1;
                    i++;
                    break;
                case '-':
                    isign = -1;
                    i++;
                    break;
                case '*':
                    i++;
                    isign = 1;
                    multflag = 1;
                    break;
                case ']':
                    i++;
                    num[inum] = 0;
                    if (multflag == 0) {
                        add = atoi(num)*isign + k;
                    } else {
                        add = atoi(num)*k;
                        multflag = 0;
                    }
                    snprintf(num, sizeof(num), "%d", add);
                    m = (int32)strlen(num);
                    for (int32 j = 0; j < m; j++) {
                        new[inew] = num[j];
                        inew++;
                    }
                    ok = 0;
                    break;
                default:
                    i++;
                    num[inum] = ch;
                    inum++;
                    break;
                }
            }
        } else {
            new[inew] = ch;

            i++;
            inew++;
        }

        if (i >= n) {
            break;
        }
    }
    new[inew] = 0;
    return;
}

void
form_ode_add_comment(char *s) {
    char text[256];
    char action[256];
    char ch;
    int32 n = (int32)strlen(s);
    int32 j1 = 0;
    int32 ja = 0;
    int32 noact = 1;
    if (n_comments >= MAXCOMMENTS) {
        return;
    }
    for (int32 i = 0; i < n; i++) {
        if (s[i] == '{') {
            j1 = i + 1;
            noact = 0;
            break;
        }
    }
    if (noact) {
        comments[n_comments].text = xmalloc(strlen(s) + 1);
        strcpy(comments[n_comments].text, s + 1);
        comments[n_comments].aflag = 0;
    } else {
        text[0] = '*';
        text[1] = ' ';
        action[0] = '$';
        action[1] = ' ';
        ja = 2;
        for (int32 i = j1; i < n; i++) {
            ch = s[i];
            if (ch == ',') {
                action[ja] = ' ';
                ja++;
            }
            if (ch == '}') {
                action[ja] = ' ';
                action[ja + 1] = 0;
                j1 = i + 1;
                break;
            }
            if (ch != ',') {
                action[ja] = ch;
                ja++;
            }
        }
        ja = 2;
        for (int32 i = j1; i < n; i++) {
            text[ja] = s[i];
            ja++;
        }
        text[ja] = 0;
        comments[n_comments].text = xmalloc(strlen(text) + 1);
        strcpy(comments[n_comments].text, text);
        comments[n_comments].action = xmalloc(strlen(action) + 1);
        strcpy(comments[n_comments].action, action);
        comments[n_comments].aflag = 1;
    }
    ggets_plintf("text=%s \n", comments[n_comments].text);
    if (comments[n_comments].aflag == 1) {
        ggets_plintf("action=%s \n", comments[n_comments].action);
    }
    n_comments++;
    return;
}

void
form_ode_advance_past_first_word(char **sptr) {
    /* changes the string pointed to by sptr to start after the end of the
       string... this may seem odd, but it has to do with avoiding \0's added by
       strtok */
    int32 len = (int32)strlen(*sptr);
    (*sptr) += len + 1;
    return;
}

char *
form_ode_new_string2(char *old, int32 length) {
    /*cout << "form_ode_new_string2(\"" << old << "\", " << length << ")\n"; */
    char *s = xmalloc((usize)(length + 1)*sizeof(char));
    strncpy(s, old, (usize)length);
    s[length] = '\0';
    if (length > 0 && s[length - 1] == ',') {
        s[length - 1] = '\0';
    }
    /* printf("s = %s; length = %d\n", s, length); */
    return s;
}

char *
form_ode_get_next2(char **tokens_ptr) {
    /* grabs (a copy of) the next block of the form var = val, ending with a \n,
     * space, or comma */
    /* importantly, this supports white space around the equal sign */
    /* returns NULL if no more text */
    /* advances tokens_ptr */
    /* modified 2012-10-12 to also work if no = */
    int32 success = 0;
    int32 i = 0;
    char *tokens = *tokens_ptr;
    int32 len;
    for (; *tokens; tokens++) {
        if (!isspace(tokens[i])) {
            break;
        }
    }
    if (!(*tokens)) {
        (*tokens_ptr) = tokens;
        return NULL;
    }
    len = (int32)strlen(tokens);
    /* advance past space/the equal sign/comma */
    success = 0;
    for (i = 1; i < len; i++) {
        if (tokens[i] == '=' || isspace(tokens[i]) || tokens[i] == ',') {
            success = 1;
            break;
        }
    }

    if (!success) {
        /* this is either a variable alone or a syntax error */
        *tokens_ptr = &tokens[len];
        return form_ode_new_string2(tokens, len);
    }

    /* advance past any spaces */
    success = 0;
    for (; i < len; i++) {
        if (!isspace(tokens[i])) {
            success = 1;
            break;
        }
    }

    if (!success) {
        /* this is either a variable alone or a syntax error */
        *tokens_ptr = &tokens[len];
        return form_ode_new_string2(tokens, len);
    }

    if (tokens[i] != '=') {
        if (tokens[i] == ',') {
            *tokens_ptr = &tokens[i + 1];
        } else {
            *tokens_ptr = &tokens[i];
        }
        return form_ode_new_string2(tokens, i);
    }

    /* advance until the first non-space */
    success = 0;
    for (i = i + 1; i < len; i++) {
        if (!isspace(tokens[i])) {
            success = 1;
            break;
        }
    }
    if (!success) {
        /* also a syntax error */
        *tokens_ptr = &tokens[len];
        return form_ode_new_string2(tokens, len);
    }

    /* advance past the nonspaces and non-commas */
    for (; i < len; i++) {
        if (isspace(tokens[i]) || tokens[i] == ',') {
            break;
        }
    }

    /* advance past any spaces */
    for (; i < len; i++) {
        if (!isspace(tokens[i])) {
            break;
        }
    }

    /* advance past a comma, if any */
    if (i < len) {
        if (tokens[i] == ',') {
            i++;
        }
    }

    /* advance the pointer to point to the next character, or the null character
     * if no more */
    *tokens_ptr = &tokens[i];
    return form_ode_new_string2(tokens, i);
}

void
form_ode_strcpy_trim(char *dest, char *source) {
    /* like strcpy, except removes leading and trailing whitespace */
    int32 len;
    int32 i;
    while (*source && isspace(*source)) {
        source++;
    }
    len = (int32)strlen(source);
    for (i = len - 1; i >= 0; i--) {
        if (!isspace(source[i])) {
            break;
        }
    }
    strncpy(dest, source, (usize)i + 1);
    dest[i + 1] = '\0';
    return;
}

void
form_ode_strncpy_trim(char *dest, char *source, int32 n) {
    int32 i;
    while (*source && isspace(*source)) {
        source++;
        n--;
    }
    for (i = n - 1; i >= 0; i--) {
        if (!isspace(source[i])) {
            break;
        }
    }
    if (i + 1 > n) {
        i = n - 1;
    }
    strncpy(dest, source, (usize)i + 1);
    dest[i + 1] = '\0';
    return;
}
