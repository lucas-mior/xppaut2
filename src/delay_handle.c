#include "functions.h"
#include "parserslow.h"

#include <stdlib.h>
#include <string.h>
/*   This handles the delay stuff    */

#include <stdio.h>
#include <math.h>
#include "xpplim.h"
#include "integers.h"

double alpha_max = 2;
double omega_max = 2;
int32 delay_flag = 0;
int32 delay_grid;

static double *delay_work;
static int32 LatestDelay;
static int32 MaxDelay;

int32 NDelay;
int32 del_stab_flag;
int32 which_delay;
int32 delay_grid = 1000;

double variable_shift[2][MAX_ODE];
double delay_list[MAX_DELAY];

static void delay_handle_polint(double *xa, double *ya, int32 n, double x, double *y, double *dy);

double
delay_handle_stab_eval(
    // this returns appropriate values for delay jacobian
    double delay2, int32 var) {

    if (del_stab_flag == 0) {
        // search for all delays
        for (int32 i = 0; i < NDelay; i++) {
            if (delay2 == delay_list[i]) {
                return GETVAR(var);
            }
        }
        delay_list[NDelay] = delay2;
        NDelay++;
        return GETVAR(var);
    }
    //  now we must determine the value to return
    //  del_stab_flag =-1
    for (int32 i = 0; i < NDelay; i++) {
        if (delay2 == delay_list[i]) {
            if (i == which_delay) {
                return variable_shift[1][var - 1];
            }
        }
    }
    return variable_shift[0][var - 1];
}

int32
delay_handle_alloc_delay(double big) {
    int32 n;

    n = (int32)(big / fabs(delta_t)) + 1;

    MaxDelay = n;
    LatestDelay = 1;
    delay_flag = 0;
    delay_work = malloc((usize)(n*NODE)*sizeof(*delay_work));

    if (delay_work == NULL) {
        ggets_err_msg("Could not allocate memory for Delay");
        return 0;
    }

    memset(delay_work, 0, (usize)(n*NODE));
    delay_flag = 1;
    NDelay = 0;
    which_delay = -1;
    del_stab_flag = 1;
    for (int32 i = 0; i < n*(NODE); i++) {
        delay_work[i] = 0.0;
    }
    return 1;
}

void
delay_handle_free_delay(void) {
    if (delay_flag) {
        free(delay_work);
    }
    delay_flag = 0;
    return;
}

void
delay_handle_stor_delay(double *y) {
    int32 in;
    int32 nodes = NODE;

    if (delay_flag == 0) {
        return;
    }
    --LatestDelay;
    if (LatestDelay < 0) {
        LatestDelay += MaxDelay;
    }
    in = LatestDelay*(nodes);
    for (int32 i = 0; i < (nodes); i++) {
        delay_work[i + in] = y[i];
    }
    return;
}

void
delay_handle_polint(double *xa, double *ya, int32 n, double x, double *y, double *dy) {
    int32 m;
    int32 ns = 1;
    double den;
    double dif;
    double dift;
    double h0;
    double hp;
    double w;
    double c[10];
    double d[10];
    dif = fabs(x - xa[0]);
    for (int32 i = 1; i <= n; i++) {
        if ((dift = fabs(x - xa[i - 1])) < dif) {
            ns = i;
            dif = dift;
        }
        c[i - 1] = ya[i - 1];
        d[i - 1] = ya[i - 1];
    }
    *y = ya[(ns--) - 1];
    for (m = 1; m < n; m++) {
        for (int32 i = 1; i <= n - m; i++) {
            h0 = xa[i - 1] - x;
            hp = xa[i + m - 1] - x;
            w = c[i] - d[i - 1];
            if ((den = h0 - hp) == 0.0) {
                return;
            }
            den = w / den;
            d[i - 1] = hp*den;
            c[i - 1] = h0*den;
        }
        *y += (*dy = (2*ns < (n - m) ? c[ns] : d[ns-- - 1]));
    }
    return;
}

/* this is like get_delay but uses cubic interpolation */
double
delay_handle_get_delay(int32 in, double tau) {
    double x = tau / fabs(delta_t);
    double dd = fabs(delta_t);
    double y;
    double ya[4];
    double xa[4];
    double dy;
    int32 n1 = (int32)x;
    int32 n2 = n1 + 1;
    int32 nodes = NODE;
    int32 n0 = n1;
    int32 n3 = n2 + 1;
    int32 i0;
    int32 i1;
    int32 i2;
    int32 i3;

    if (tau < 0.0 || tau > delay) {
        ggets_err_msg("Delay negative or too large");
        integrate_stop_integration();
        return 0.0;
    }
    if (tau == 0.0) {  // check fro zero delay and ignore the rest
        return delay_work[in + nodes*(LatestDelay % MaxDelay)];
    }
    xa[1] = n1*dd;
    xa[0] = xa[1] - dd;
    xa[2] = xa[1] + dd;
    xa[3] = xa[2] + dd;
    i1 = (n1 + LatestDelay) % MaxDelay;
    i2 = (n2 + LatestDelay) % MaxDelay;
    i0 = (n0 + LatestDelay) % MaxDelay;
    i3 = (n3 + LatestDelay) % MaxDelay;
    if (i1 < 0) {
        i1 += MaxDelay;
    }
    if (i2 < 0) {
        i2 += MaxDelay;
    }
    if (i3 < 0) {
        i3 += MaxDelay;
    }
    if (i0 < 0) {
        i0 += MaxDelay;
    }

    ya[1] = delay_work[in + (nodes)*i1];
    ya[2] = delay_work[in + (nodes)*i2];
    ya[0] = delay_work[in + (nodes)*i0];
    ya[3] = delay_work[in + (nodes)*i3];
    delay_handle_polint(xa, ya, 4, tau, &y, &dy);

    return y;
}

/*  Handling of the initial data  */
int32
delay_handle_do_init_delay(double big) {
    double t = T0;
    double old_t;
    double y[MAX_ODE];
    int32 nt;
    int32 len;

    int32 *del_form[MAX_ODE];
    nt = (int32)(big / fabs(delta_t));
    NCON = NCON_START;
    NSYM = NSYM_START;
    for (int32 i = 0; i < (NODE); i++) {
        del_form[i] = (int32 *)calloc(200, sizeof(int32));
        if (del_form[i] == NULL) {
            ggets_err_msg("Failed to allocate delay formula ...");
            for (int32 j = 0; j < i; j++) {
                free(del_form[j]);
            }
            NCON = NCON_START;
            NSYM = NSYM_START;
            return 0;
        }

        if (parserslow_add_expr(delay_string[i], del_form[i], &len)) {
            ggets_err_msg("Illegal delay expression");
            for (int32 j = 0; j <= i; j++) {
                free(del_form[j]);
            }
            NCON = NCON_START;
            NSYM = NSYM_START;
            return 0;
        }
    }  //  Okay all formulas are cool...
    LatestDelay = 1;

    get_val("t", &old_t);

    for (int32 i = nt; i >= 0; i--) {
        t = T0 - fabs(delta_t)*i;
        set_val("t", t);
        for (int32 j = 0; j < (NODE); j++) {
            y[j] = evaluate(del_form[j]);
        }
        delay_handle_stor_delay(y);
    }
    for (int32 j = 0; j < (NODE); j++) {
        free(del_form[j]);
    }
    NCON = NCON_START;
    NSYM = NSYM_START;
    set_val("t", old_t);
    return 1;
}
