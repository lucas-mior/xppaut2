#ifndef _comline_h_
#define _comline_h_
#include "integers.h"

typedef struct {
    char *name;
    char *does;
    uint32 use;
} INTERN_SET;

typedef struct {
    char *name;
    struct SET_NAME *next;
} SET_NAME;

int32 is_set_name(SET_NAME *set, char *nam);
SET_NAME *add_set(SET_NAME *set, char *nam);
SET_NAME *rm_set(SET_NAME *set, char *nam);
void do_comline(int32 argc, char **argv);
int32 if_needed_select_sets(void);
int32 if_needed_load_set(void);
int32 if_needed_load_par(void);
int32 if_needed_load_ic(void);
int32 if_needed_load_ext_options(void);
int32 parse_it(char *com);

#endif
