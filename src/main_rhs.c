#include "functions.h"
#include "integers.h"
#include "parserslow.h"

int32
main(int32 argc, char **argv) {
    do_main(argc, argv);
}

void
main_rhs_extra(double *y__y, double t, int32 nod, int32 neq) {
    if (nod >= neq)
        return;
    SETVAR(0, t);
    for (int32 i = 0; i < nod; i++)
        SETVAR(i + 1, y__y[i]);
    for (int32 i = nod + FIX_VAR; i < nod + FIX_VAR + NMarkov; i++)
        SETVAR(i + 1, y__y[i - FIX_VAR]);
    for (int32 i = nod; i < nod + FIX_VAR; i++)
        SETVAR(i + 1, evaluate(my_ode[i]));

    for (int32 i = nod + NMarkov; i < neq; i++)
        y__y[i] = evaluate(my_ode[i + FIX_VAR - NMarkov]);
    return;
}

void
main_rhs_set_fix(double t, double *y) {
    SETVAR(0, t);
    for (int32 i = 0; i < NODE; i++)
        SETVAR(i + 1, y[i]);
    for (int32 i = 0; i < NMarkov; i++)
        SETVAR(i + 1 + NODE + FIX_VAR, y[i + NODE]);
    for (int32 i = NODE; i < NODE + FIX_VAR; i++)
        SETVAR(i + 1, evaluate(my_ode[i]));
    simplenet_eval_all_nets();

    extra_do_in_out();
    return;
}

int32
main_rhs(double t, double *y, double *ydot, int32 neq) {
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
main_rhs_update_based_on_current(void) {
    for (int32 i = NODE; i < NODE + FIX_VAR; i++)
        SETVAR(i + 1, evaluate(my_ode[i]));

    simplenet_eval_all_nets();
    extra_do_in_out();
    return;
}

void
main_rhs_fix_only(void) {
    for (int32 i = NODE; i < NODE + FIX_VAR; i++)
        SETVAR(i + 1, evaluate(my_ode[i]));
    return;
}

void
main_rhs_only(double *ydot) {
    for (int32 i = 0; i < NODE; i++) {
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
