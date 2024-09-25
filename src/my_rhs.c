#include "functions.h"
#include "integers.h"
#include <stdlib.h>
#include "xpplim.h"
#include "parserslow.h"

int32
main(int32 argc, char **argv) {
    do_main(argc, argv);
}

void
my_rhs_extra(double *y__y, double t, int32 nod, int32 neq) {
    int32 i;
    if (nod >= neq)
        return;
    SETVAR(0, t);
    for (i = 0; i < nod; i++)
        SETVAR(i + 1, y__y[i]);
    for (i = nod + FIX_VAR; i < nod + FIX_VAR + NMarkov; i++)
        SETVAR(i + 1, y__y[i - FIX_VAR]);
    for (i = nod; i < nod + FIX_VAR; i++)
        SETVAR(i + 1, evaluate(my_ode[i]));

    for (i = nod + NMarkov; i < neq; i++)
        y__y[i] = evaluate(my_ode[i + FIX_VAR - NMarkov]);
    return;
}

void
set_fix_rhs(double t, double *y) {
    int32 i;
    SETVAR(0, t);
    for (i = 0; i < NODE; i++)
        SETVAR(i + 1, y[i]);
    for (i = 0; i < NMarkov; i++)
        SETVAR(i + 1 + NODE + FIX_VAR, y[i + NODE]);
    for (i = NODE; i < NODE + FIX_VAR; i++)
        SETVAR(i + 1, evaluate(my_ode[i]));
    simplenet_eval_all_nets();

    extra_do_in_out();
    return;
}

int32
my_rhs(double t, double *y, double *ydot, int32 neq) {
    (void)neq;
    SETVAR(0, t);
    for (int32 i = 0; i < NODE; i++)
        SETVAR(i + 1, y[i]);

    for (int32 i = NODE; i < NODE + FIX_VAR; i++) {
        SETVAR(i + 1, evaluate(my_ode[i]));
    }
    simplenet_eval_all_nets();

    dae_fun_do_daes();

    extra_do_in_out();
    for (int32 i = 0; i < NODE; i++) {
        ydot[i] = evaluate(my_ode[i]);
    }

    return 1;
}

void
update_based_on_current(void) {
    int32 i;
    for (i = NODE; i < NODE + FIX_VAR; i++)
        SETVAR(i + 1, evaluate(my_ode[i]));

    simplenet_eval_all_nets();
    extra_do_in_out();
    return;
}

void
fix_only(void) {
    int32 i;
    for (i = NODE; i < NODE + FIX_VAR; i++)
        SETVAR(i + 1, evaluate(my_ode[i]));
    return;
}

void
rhs_only(double *ydot) {
    int32 i;
    for (i = 0; i < NODE; i++) {
        ydot[i] = evaluate(my_ode[i]);
    }
    return;
}

/*
 * This is the order in which quantities are evaluated
 *
 * 1.  Fixed variables
 * 2.  network stuff
 * 3.  DAEs
 * 4.  External C code
 * 5.  RH sides of the ODEs

 * For Auxilliary stuff

 * external C code is not evaluated but fixed are */
