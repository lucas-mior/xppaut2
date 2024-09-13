#include "functions.h"

#include <stdlib.h>
#include <string.h>
#include "parserslow.h"
#include "integers.h"

/* Derived parameter stuff !!  */
#define MAXDERIVED 200
extern double constants[];
extern int32 NCON;
typedef struct {
    int32 index, *form;
    char *rhs;
    double value;
} DERIVED;

DERIVED derived[MAXDERIVED];
int32 nderived = 0;

/* clean up derived stuff */
void
free_derived(void) {
    int32 i;
    for (i = 0; i < nderived; i++) {
        free(derived[i].form);
        free(derived[i].rhs);
    }
    nderived = 0;
    return;
}

/* This compiles all of the formulae
It is called only once during the session
*/
int32
compile_derived(void) {
    int32 i, k;
    int32 f[256], n;
    for (i = 0; i < nderived; i++) {
        if (add_expr(derived[i].rhs, f, &n) == 1) {
            plintf(" Bad right-hand side for derived parameters \n");
            return 1;
        }
        derived[i].form = malloc(sizeof(int32) * (n + 2));
        for (k = 0; k < n; k++)
            derived[i].form[k] = f[k];
    }
    evaluate_derived();
    return 0;
}

/* This evaluates all derived quantities in order of definition
called before any integration or numerical computation
and after changing parameters and constants
*/
void
evaluate_derived(void) {
    int32 i;
    for (i = 0; i < nderived; i++) {
        derived[i].value = evaluate(derived[i].form);
        constants[derived[i].index] = derived[i].value;
    }
    return;
}

/* this adds a derived quantity  */
int32
add_derived(char *name, char *rhs) {
    int32 n = strlen(rhs) + 2;
    int32 i0;
    if (nderived >= MAXDERIVED) {
        plintf(" Too many derived constants! \n");
        return 1;
    }
    i0 = nderived;
    derived[i0].rhs = malloc(n);
    /* save the right hand side */
    strcpy(derived[i0].rhs, rhs);
    /* this is the constant to which it addresses */
    derived[i0].index = NCON;
    /* add the name to the recognized symbols */
    plintf(" derived constant[%d] is %s = %s\n", NCON, name, rhs);
    nderived++;
    return add_con(name, 0.0);
}
