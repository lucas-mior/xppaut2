
#ifndef _form_ode_h
#define _form_ode_h
#include "integers.h"

#include "xpplim.h"
#include "newpars.h"
#include <stdio.h>

#define MAXVNAM 33
#define MAXLINES 5000

/*void break_up_list(char *rhs);
void compile_em();
void free_varinfo();
void remove_blanks(char *s1);
void read_a_line(FILE *fp,char *s);

void subsk(char *big,char *new,int32 k,int32 flag);
void free_comments();

void add_comment(char *s);
void init_varinfo();
void add_varinfo(int32 type,char *lhs,char *rhs,int32 nargs,char
args[MAXARG][NAMLEN+1]); void stor_internopts(char *s1);
*/

typedef struct {
    char *name, *value;
} FIXINFO;

int32 make_eqn(void);
void strip_saveqn(void);
int32 disc(char *string);
void dump_src(void);
void dump_comments(void);
void format_list(char **s, int32 n);
int32 get_a_filename(char *filename, char *wild);
void list_em(char *wild);
int32 read_eqn(void);
int32 get_eqn(FILE *fptr);
int32 compiler(char *bob, FILE *fptr);
void list_upar(void);
void welcome(void);
void show_syms(void);
void take_apart(char *bob, double *value, char *name);
char *get_first(char *string, char *src);
char *get_next(char *src);
void find_ker(char *string, int32 *alt);
void pos_prn(char *s, int32 x, int32 y);
void clrscr(void);
int32 getuch(void);
int32 getchi(void);
int32 if_include_file(char *old, char *nf);
int32 if_end_include(char *old);
int32 do_new_parser(FILE *fp, char *first, int32 nnn);
void create_plot_list(void);
void add_only(char *s);
void break_up_list(char *rhs);
int32 find_the_name(char list[1949][33], int32 n, char *name);
void compile_em(void);
int32 formula_or_number(char *expr, double *z);
void strpiece(char *dest, char *src, int32 i0, int32 ie);
int32 parse_a_string(char *s1, VAR_INFO *v);
void init_varinfo(void);
void add_varinfo(int32 type, char *lhs, char *rhs, int32 nargs,
                 char args[20][13 + 1]);
void free_varinfo(void);
int32 extract_ode(char *s1, int32 *ie, int32 i1);
int32 strparse(char *s1, char *s2, int32 i0, int32 *i1);
int32 extract_args(char *s1, int32 i0, int32 *ie, int32 *narg,
                   char args[20][13 + 1]);
int32 find_char(char *s1, char *s2, int32 i0, int32 *i1);
int32 next_nonspace(char *s1, int32 i0, int32 *i1);
void remove_blanks(char *s1);
void read_a_line(FILE *fp, char *s);
int32 search_array(char *old, char *new, int32 *i1, int32 *i2, int32 *flag);
int32 check_if_ic(char *big);
int32 not_ker(char *s, int32 i);
int32 is_comment(char *s);
void subsk(char *big, char *new, int32 k, int32 flag);
void keep_orig_comments(void);
void default_comments(void);
void free_comments(void);
void new_comment(FILE *f);
void add_comment(char *s);

/* for parsing par, init with whitespace correctly */
char *new_string2(char *old, int32 length);
void advance_past_first_word(char **sptr);
char *get_next2(char **tokens_ptr);
void strcpy_trim(char *dest, char *source);
void strncpy_trim(char *dest, char *source, int32 n);
#endif
