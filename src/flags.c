#include "functions.h"
#include "integers.h"
#include <stdbool.h>

#include "cv2.h"
#include "parserslow.h"

#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "xpplim.h"
#define MY_DBL_EPS 5e-16

/*  this is a new (Summer 1995) addition to XPP that allows one to
    do things like delta functions and other discontinuous
    stuff.

    The conditions are set up as part of the "ODE" file:

global sign condition {event1;....;eventn}
global sign {condition} {event1;....;eventn}

the {} and ;  are required for the events

condition is anything that when it evaluates to 0 means the flag should be
set.  The sign is like in Poincare maps, thus let C(t1) and C(t2) be
the value of the condition at t1  and t2.
sign = 0 ==>  just find when C(t)=0
sign = 1 ==>  C(t1)<0 C(t2)>0
sign = -1==>  C(t1)>0 C(t2)<0

To get the time of the event, we use linear interpolation:

 t* = t1 + (t2-t1)
           -------   (0-C(t1))
          C(t2)-C(t1)
This yields  the variables, etc at that time

Now what are the events:

They are of the form:
   name = expression  the variable  <name> is replaced by the value of
  <expression>

Note that there may be several "conditions" defined and that
these must also be checked to see if they have been switched
and in what order.  This is particularly true for "delta" function
type things.

Here is a simple example -- the kicked cycle:
dx/dt = y
dy/dy = -x -c y

if y=0 and y goes from pos to neg then x=x+b
here is how it would work:

global -1 y {x=x+b}

Here is Tysons model:

global -1 u-.2 {m=.5*m}

*/

/*

type =0 variable
type =1 parameter
type =2 output
type =3 halt
*/

#define MAX_EVENTS 20 /*  this is the maximum number of events per flag */

typedef struct {
    double f0;
    double f1;
    double tstar;
    int32 lhs[MAX_EVENTS];
    double vrhs[MAX_EVENTS];
    char lhsname[MAX_EVENTS][11];
    char *rhs[MAX_EVENTS];
    int32 *comrhs[MAX_EVENTS];
    char *cond;
    int32 *comcond;
    int32 sign;
    int32 nevents;
    int32 hit, type[MAX_EVENTS];
    int32 anypars;
    int32 nointerp;
} FLAG;

#define IC 2
#define PARAM 1
/* #define Set_ivar(a,b) variables[(a)]=(b) */
static FLAG flag[MAX_FLAG];
int32 NFlags = 0;

double STOL = 1.e-10;
extern double variables[];
extern int32 NVAR;

int32
flags_add_global(char *cond, int32 sign, char *rest) {
    char temp[256];
    int32 nevents, ii, k, l, lt, j = NFlags;
    char ch;
    if (NFlags >= MAX_FLAG) {
        ggets_plintf("Too many global conditions\n");
        return 1;
    }
    l = (int32)strlen(cond);
    flag[j].cond = xmalloc((usize)l + 1);
    strcpy(flag[j].cond, cond);
    nevents = 0;
    flag[j].lhsname[0][0] = 0;
    k = 0;
    l = (int32)strlen(rest);
    for (ii = 0; ii < l; ii++) {
        ch = rest[ii];
        if (ch == '{' || ch == ' ')
            continue;
        if (ch == '}' || ch == ';') {
            if (nevents == MAX_EVENTS) {
                printf(" Too many events per flag \n");
                return 1;
            }
            temp[k] = 0;
            lt = (int32)strlen(temp);
            if (flag[j].lhsname[nevents][0] == 0) {
                printf(" No event variable named for %s \n", temp);
                return 1;
            }
            flag[j].rhs[nevents] = xmalloc((usize)lt + 1);
            strcpy(flag[j].rhs[nevents], temp);
            nevents++;
            k = 0;
            if (ch == '}')
                break;
            continue;
        }
        if (ch == '=') {
            temp[k] = 0;
            strcpy(flag[j].lhsname[nevents], temp);

            k = 0;
            if (nevents < MAX_EVENTS - 1)
                flag[j].lhsname[nevents + 1][0] = 0;
            continue;
        }

        temp[k] = ch;
        k++;
    }
    if (nevents == 0) {
        ggets_plintf(" No events for condition %s \n", cond);
        return 1;
    }
    /*  we now have the condition, the names, and the formulae */
    flag[j].sign = sign;
    flag[j].nevents = nevents;
    NFlags++;
    return 0;
}

void
flags_show(void) {
    /* uncomment for debugging */
    /*
    for(i=0;i<NFlags;i++){
      n=flag[i].nevents;
      ggets_plintf(" Flag %d has sign %d and %d events and condition %s \n",
             i+1,flag[i].sign,n,flag[i].cond);
      for(j=0;j<n;j++)
        ggets_plintf("%d:  %s [%d] = %s \n",j+1,flag[i].lhsname[j],flag[i].lhs[j],
               flag[i].rhs[j]);
    }
    */
    return;
}

int32
flags_compile(void) {
    int32 j;
    int32 i, k, index, nc;
    int32 command[256];
    if (NFlags == 0)
        return 0;
    for (j = 0; j < NFlags; j++) {
        if (add_expr(flag[j].cond, command, &nc)) {
            ggets_plintf("Illegal global condition:  %s\n", flag[j].cond);
            return 1;
        }
        flag[j].anypars = 0;
        flag[j].nointerp = 0;
        flag[j].comcond = xmalloc(sizeof(*(flag[j].comcond))*(usize)(nc + 1));
        for (k = 0; k <= nc; k++)
            flag[j].comcond[k] = command[k];
        for (i = 0; i < flag[j].nevents; i++) {
            index = init_conds_find_user_name(IC, flag[j].lhsname[i]);
            if (index < 0) {
                index = init_conds_find_user_name(PARAM, flag[j].lhsname[i]);
                if (index < 0) {
                    if (strcasecmp(flag[j].lhsname[i], "out_put") == 0) {
                        flag[j].type[i] = 2;
                        flag[j].lhs[i] = 0;
                    } else {
                        if (strcasecmp(flag[j].lhsname[i], "arret") == 0) {
                            flag[j].type[i] = 3;
                            flag[j].lhs[i] = 0;

                        } else {
                            if (strcasecmp(flag[j].lhsname[i], "no_interp") ==
                                0) {
                                flag[j].nointerp = 1;
                                flag[j].type[i] = 0;
                                flag[j].lhs[i] = 0;
                            }

                            else {
                                ggets_plintf(" <%s> is not a valid "
                                       "variable/parameter name \n",
                                       flag[j].lhsname[i]);
                                return 1;
                            }
                        }
                    }
                } else {
                    flag[j].lhs[i] = index;
                    flag[j].type[i] = 1;
                    flag[j].anypars = 1;
                }
            } else {
                flag[j].lhs[i] = index;
                flag[j].type[i] = 0;
            }
            if (add_expr(flag[j].rhs[i], command, &nc)) {
                printf("Illegal event %s for global %s\n", flag[j].rhs[i],
                       flag[j].cond);
                return 1;
            }
            flag[j].comrhs[i] =
                xmalloc(sizeof(*(flag[j].comrhs[i]))*(usize)(nc + 1));
            for (k = 0; k <= nc; k++)
                flag[j].comrhs[i][k] = command[k];
        }
    }
    return 0;
}

/*  here is the shell code for a loop around  integration step  */

int32
one_flag_step(double *yold, double *ynew, int32 *istart, double told,
              double *tnew, int32 neq, double *s) {
    double dt = *tnew - told;
    double f0, f1, tol, tolmin = 1e-10;
    double smin = 2;
    int32 sign, i, j, in, ncycle = 0, newhit, nevents;

    if (NFlags == 0)
        return 0;
    /* printf("dt=%g yold= %g ynew = %g \n",dt,yold[0],ynew[0]); */
    /*  if(abs(dt)<MY_DBL_EPS) return 0;  */
    for (i = 0; i < NFlags; i++) {
        flag[i].tstar = 2.0;
        flag[i].hit = 0;
    }
    /* If this is the first call, then need f1  */
    if (*istart == 1) {
        for (i = 0; i < neq; i++)
            SETVAR(i + 1, yold[i]);
        SETVAR(0, told);
        for (i = 0; i < NFlags; i++)
            *istart = 0;
    }
    for (i = 0; i < NFlags; i++) {
        sign = flag[i].sign;
        flag[i].f0 = flag[i].f1;
        f0 = flag[i].f0;
        for (j = 0; j < neq; j++)
            SETVAR(j + 1, ynew[j]);
        SETVAR(0, *tnew);
        f1 = evaluate(flag[i].comcond);
        flag[i].f1 = f1;
        tol = fabs(f1 - f0);
        switch (sign) {
        case 1:
            if ((((f0 < 0.0) && (f1 > 0.0)) || ((f0 < 0.0) && (f1 > 0.0))) &&
                tol > tolmin) {
                flag[i].hit = ncycle + 1;
                flag[i].tstar = f0 / (f0 - f1);
            }
            break;
        case -1:
            if (f0 > 0.0 && f1 <= 0.0 && tol > tolmin) {
                flag[i].hit = ncycle + 1;
                flag[i].tstar = f0 / (f0 - f1);
            }
            break;
        case 0:
            /* if(f1==0.0){ */
            if (fabs(f1) < MY_DBL_EPS) {
                flag[i].hit = ncycle + 1;
                flag[i].tstar = told;
            }
            break;
        default:
            break;
        }
        if (flag[i].nointerp == 1) {
            flag[i].tstar = 1.0;
        }

        if (smin > flag[i].tstar)
            smin = flag[i].tstar;
    } /* run through flags */

    if (smin < STOL)
        smin = STOL;
    else
        smin = (1 + STOL)*smin;
    if (smin > 1.0)
        return 0;

    *tnew = told + dt*smin;
    SETVAR(0, *tnew);
    for (i = 0; i < neq; i++) {
        ynew[i] = yold[i] + smin*(ynew[i] - yold[i]);
        SETVAR(i + 1, ynew[i]);
    }
    for (i = 0; i < NFlags; i++)
        flag[i].f0 = evaluate(flag[i].comcond);
    while (true) { /* run through all possible events  */
        ncycle++;
        newhit = 0;
        for (i = 0; i < NFlags; i++) {
            nevents = flag[i].nevents;
            if (flag[i].hit == ncycle && flag[i].tstar <= smin) {
                for (j = 0; j < nevents; j++) {
                    flag[i].vrhs[j] = evaluate(flag[i].comrhs[j]);
                    in = flag[i].lhs[j];
                    if (flag[i].type[j] == 0)
                        SETVAR(in + 1, flag[i].vrhs[j]);
                }
            }
        }
        /* printf("step 7 \n");
        for(i=0;i<neq;i++)
        printf("%d %g %g\n",i,ynew[i],GETVAR(i+1)); */
        for (i = 0; i < NFlags; i++) {
            nevents = flag[i].nevents;
            if (flag[i].hit == ncycle && flag[i].tstar <= smin) {
                for (j = 0; j < nevents; j++) {
                    in = flag[i].lhs[j];
                    if (flag[i].type[j] == 0) {
                        ynew[in] = flag[i].vrhs[j];
                        /* SETVAR(in+1,ynew[in]); if this screws up */
                    } else {
                        if (flag[i].type[j] == 1)
                            set_val(upar_names[in], flag[i].vrhs[j]);
                        else {
                            if ((flag[i].type[j] == 2) && (flag[i].vrhs[j] > 0))
                                integrate_send_output(ynew, *tnew);
                            if ((flag[i].type[j] == 3) && (flag[i].vrhs[j] > 0))
                                integrate_send_halt();
                        }
                    }
                }
                if (flag[i].anypars) {
                    derived_evaluate();
                    init_conds_redraw_params();
                }
            }
        }

        for (i = 0; i < neq; i++) {
            /* printf("step 8 %d %g %g\n",i,ynew[i],GETVAR(i+1)); */
            /*  SETVAR(i+1,ynew[i]); */
            ynew[i] = GETVAR(i + 1); /* if this screws up */
            /*      printf("step 9 %d %g %g\n",i,ynew[i],GETVAR(i+1)); */
        }
        for (i = 0; i < NFlags; i++) {
            flag[i].f1 = evaluate(flag[i].comcond);
            if (flag[i].hit > 0)
                continue; /* already hit so dont do anything */
            f1 = flag[i].f1;
            sign = flag[i].sign;
            f0 = flag[i].f0;
            tol = fabs(f1 - f0);
            switch (sign) {
            case 1:
                if (f0 <= 0.0 && f1 >= 0.0 && tol > tolmin) {
                    flag[i].tstar = smin;
                    flag[i].hit = ncycle + 1;
                    newhit = 1;
                }
                break;
            case -1:
                if (f0 >= 0.0 && f1 <= 0.0 && tol > tolmin) {
                    flag[i].tstar = smin;
                    flag[i].hit = ncycle + 1;
                    newhit = 1;
                }
                break;
            case 0:
                if (f0*f1 <= 0 && (f1 != 0 || f0 != 0) && tol > tolmin) {
                    flag[i].tstar = smin;
                    flag[i].hit = ncycle + 1;
                    newhit = 1;
                }
                break;
            default:
                break;
            }
        }
        if (newhit == 0)
            break;
    }

    *s = smin;
    return 1;
}

/*  here are the ODE drivers */

int32
one_flag_step_symp(double *y, double dt, double *work, int32 neq, double *tim,
                   int32 *istart) {
    double yold[MAX_ODE], told;
    int32 i;
    int32 hit;
    double s, dtt = dt;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < neq; i++)
            yold[i] = y[i];
        told = *tim;
        odesol2_one_step_symp(y, dtt, work, neq, tim);
        if ((hit = one_flag_step(yold, y, istart, told, tim, neq, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;
        dtt = (1 - s)*dt;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard?? ");
            ggets_plintf("smin=%g\n", s);
            break;
        }
    }

    return 1;
}

int32
one_flag_step_euler(double *y, double dt, double *work, int32 neq, double *tim,
                    int32 *istart) {
    double yold[MAX_ODE], told;
    int32 i;
    int32 hit;
    double s, dtt = dt;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < neq; i++)
            yold[i] = y[i];
        told = *tim;
        odesol2_one_step_euler(y, dtt, work, neq, tim);
        if ((hit = one_flag_step(yold, y, istart, told, tim, neq, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;
        dtt = (1 - s)*dt;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard?? ");
            ggets_plintf("smin=%g\n", s);
            break;
        }
    }

    return 1;
}

int32
one_flag_step_discrete(double *y, double dt, double *work, int32 neq,
                       double *tim, int32 *istart) {
    double yold[MAX_ODE], told;
    int32 i;
    int32 hit;
    double s, dtt = dt;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < neq; i++)
            yold[i] = y[i];
        told = *tim;
        odesol2_one_step_discrete(y, dtt, work, neq, tim);
        if ((hit = one_flag_step(yold, y, istart, told, tim, neq, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;
        dtt = (1 - s)*dt;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard?? ");
            ggets_plintf("smin=%g\n", s);
            break;
        }
    }
    return 1;
}

int32
one_flag_step_heun(double *y, double dt, double *yval[2], int32 neq,
                   double *tim, int32 *istart) {
    double yold[MAX_ODE], told;
    int32 i;
    int32 hit;
    double s, dtt = dt;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < neq; i++)
            yold[i] = y[i];
        told = *tim;
        one_step_heun(y, dtt, yval, neq, tim);
        if ((hit = one_flag_step(yold, y, istart, told, tim, neq, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;
        dtt = (1 - s)*dt;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard? ");
            ggets_plintf(" smin=%g\n", s);
            break;
        }
    }
    return 1;
}

int32
one_flag_step_rk4(double *y, double dt, double *yval[3], int32 neq, double *tim,
                  int32 *istart) {
    double yold[MAX_ODE], told;
    int32 i;
    int32 hit;
    double s, dtt = dt;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < neq; i++)
            yold[i] = y[i];
        told = *tim;
        one_step_rk4(y, dtt, yval, neq, tim);
        if ((hit = one_flag_step(yold, y, istart, told, tim, neq, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;
        dtt = (1 - s)*dt;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard?");
            ggets_plintf("smin=%g\n", s);
            break;
        }
    }
    return 1;
}

int32
one_flag_step_gear(int32 neq, double *t, double tout, double *y, double hmin,
                   double hmax, double eps, int32 mf, double *error,
                   int32 *kflag, int32 *jstart, double *work, int32 *iwork) {
    double yold[MAX_ODE], told;
    int32 i;
    int32 hit;
    double s;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < neq; i++)
            yold[i] = y[i];
        told = *t;
        ggear(neq, t, tout, y, hmin, hmax, eps, mf, error, kflag, jstart, work,
              iwork);
        if (*kflag < 0)
            break;
        if ((hit = one_flag_step(yold, y, jstart, told, t, neq, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;
        *jstart = 0; /* for gear always reset  */
        if (*t == tout)
            break;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard? ");
            ggets_plintf("smin=%g\n", s);
            break;
        }
    }
    return 0;
}

int32
one_flag_step_rosen(double *y, double *tstart, double tfinal, int32 *istart,
                    int32 n, double *work, int32 *ierr) {
    double yold[MAX_ODE], told;
    int32 i, ok, hit;
    double s;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < n; i++)
            yold[i] = y[i];
        told = *tstart;
        ok = rosen(y, tstart, tfinal, istart, n, work, ierr);
        if (ok == -1)
            break;
        if ((hit = one_flag_step(yold, y, istart, told, tstart, n, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;

        if (*tstart == tfinal)
            break;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard? ");
            ggets_plintf("smin=%g\n", s);
            *ierr = -2;
            return 1;
        }
    }
    return 0;
}

int32
one_flag_step_dp(int32 *istart, double *y, double *t, int32 n, double tout,
                 double *tol, double *atol, int32 flag2, int32 *kflag) {
    double yold[MAX_ODE], told;
    int32 i;
    int32 hit;
    double s;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < n; i++)
            yold[i] = y[i];
        told = *t;
        dormprin(istart, y, t, n, tout, tol, atol, flag2, kflag);
        if (*kflag != 1)
            break;
        if ((hit = one_flag_step(yold, y, istart, told, t, n, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;

        if (*t == tout)
            break;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard? ");
            ggets_plintf("smin=%g\n", s);
            return 1;
        }
    }
    return 0;
}

#ifdef CVODE_YES
int32
one_flag_step_cvode(
    /* command =0 continue, 1 is start 2 finish */
    int32 *command, double *y, double *t, int32 n, double tout, int32 *kflag,
    double *atol, double *rtol) {
    double yold[MAX_ODE], told;
    int32 i, hit, neq = n;
    double s;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < neq; i++)
            yold[i] = y[i];
        told = *t;
        ccvode(command, y, t, n, tout, kflag, atol, rtol);
        if (*kflag < 0)
            break;
        if ((hit = one_flag_step(yold, y, command, told, t, neq, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;
        end_cv();
        *command = 1; /* for cvode always reset  */
        if (*t == tout)
            break;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard? ");
            ggets_plintf("smin=%g\n", s);
            return 1;
        }
    }
    return 0;
}

#endif
int32
one_flag_step_adap(double *y, int32 neq, double *t, double tout, double eps,
                   double *hguess, double hmin, double *work, int32 *ier,
                   double epjac, int32 iflag, int32 *jstart) {
    double yold[MAX_ODE], told;
    int32 i;
    int32 hit;
    double s;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < neq; i++)
            yold[i] = y[i];
        told = *t;
        gadaptive(y, neq, t, tout, eps, hguess, hmin, work, ier, epjac, iflag);
        if (*ier)
            break;
        if ((hit = one_flag_step(yold, y, jstart, told, t, neq, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;

        if (*t == tout)
            break;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard? ");
            ggets_plintf("smin=%g\n", s);
            break;
        }
    }
    return 0;
}

int32
one_flag_step_backeul(double *y, double *t, double dt, int32 neq, double *yg,
                      double *yp, double *yp2, double *ytemp, double *errvec,
                      double *jac, int32 *istart) {
    double yold[MAX_ODE], told;
    int32 i, hit, j;
    double s;
    double dtt = dt;
    int32 nstep = 0;
    while (true) {
        for (i = 0; i < neq; i++)
            yold[i] = y[i];
        told = *t;
        if ((j = one_bak_step(y, t, dtt, neq, yg, yp, yp2, ytemp, errvec,
                              jac)) != 0)
            return j;
        if ((hit = one_flag_step(yold, y, istart, told, t, neq, &s)) == 0)
            break;
        /* Its a hit !! */
        nstep++;
        dtt = (1 - s)*dt;
        if (nstep > (NFlags + 2)) {
            ggets_plintf(" Working too hard?");
            ggets_plintf("smin=%g\n", s);
            break;
        }
    }
    return 0;
}
