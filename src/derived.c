#include "functions.h"

#include <string.h>
#include "parserslow.h"
#include "integers.h"

/* Derived parameter stuff !!  */
#define MAXDERIVED 200
typedef struct {
    int32 index;
    int32 *form;
    char *rhs;
    double value;
} Derived;

static Derived derived[MAXDERIVED];
static int32 nderived = 0;

/* This compiles all of the formulae
It is called only once during the session
*/
int32
derived_compile(void) {
    int32 f[256];
    int32 n;
    for (int64 i = 0; i < nderived; i++) {
        if (parserslow_add_expr(derived[i].rhs, f, &n) == 1) {
            ggets_plintf(" Bad right-hand side for derived parameters \n");
            return 1;
        }
        derived[i].form = xmalloc(sizeof(*(derived[i].form))*(usize)(n + 2));
        for (int32 k = 0; k < n; k++) {
            derived[i].form[k] = f[k];
        }
    }
    derived_evaluate();
    return 0;
}

/* This evaluates all derived quantities in order of definition
called before any integration or numerical computation
and after changing parameters and constants
*/
void
derived_evaluate(void) {
    for (int32 i = 0; i < nderived; i++) {
        derived[i].value = evaluate(derived[i].form);
        constants[derived[i].index] = derived[i].value;
    }
    return;
}

/* this adds a derived quantity  */
int32
derived_add(char *name, char *rhs) {
    int32 n = (int32)strlen(rhs) + 2;
    int32 i0;

    if (nderived >= MAXDERIVED) {
        ggets_plintf(" Too many derived constants! \n");
        return 1;
    }
    i0 = nderived;
    derived[i0].rhs = xmalloc((usize)n);
    /* save the right hand side */
    strcpy(derived[i0].rhs, rhs);
    /* this is the constant to which it addresses */
    derived[i0].index = NCON;
    /* add the name to the recognized symbols */
    ggets_plintf(" derived constant[%d] is %s = %s\n", NCON, name, rhs);
    nderived++;
    return parserslow_add_con(name, 0.0);
}
