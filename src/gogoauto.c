#include "auto_nox.h"
#include "auto_c.h"
#include "x_auto.h"
#include "integers.h"
#include "functions.h"

FILE *fp3;
FILE *fp7;
FILE *fp9;
FILE *fp12;
int32 global_conpar_type = CONPAR_DEFAULT;
int32 global_setubv_type = SETUBV_DEFAULT;
int32 global_num_procs = 1;
int32 global_verbose_flag = 0;

int32
go_go_auto(void) {
    // this is the entry at this point, x_auto has been set
    int64 icp[NPARX2];
    double par[NPARX2];
    double thl[NPARX];
    double *thu;
    int64 iuz[100];
    double vuz[100];
    iap_type iap;
    rap_type rap;
    function_list list;
    int32 irs = x_auto.irs;

    if (irs > 0) {
        fp3 = fopen(fort3, "r");
    } else {
        fp3 = fopen(fort3, "w+");
    }

    fp7 = fopen(fort7, "w");
    fp9 = fopen(fort9, "w");

    // Initialization :

    iap.mynode = mynode();
    iap.numnodes = numnodes();
    if (iap.numnodes > 1) {
        iap.parallel_flag = 1;
    } else {
        iap.parallel_flag = 0;
    }

    // here is the feeder code from x_auto structure

    init(&iap, &rap, par, icp, thl, &thu, iuz, vuz);

    // Find restart label and determine type of restart point.
    if (iap.irs > 0) {
        int64 found = false;

        findlb(&iap, &rap, iap.irs, &(iap.nfpr), &found);
        if (!found) {
            if (iap.mynode == 0) {
                fprintf(stderr, "\nRestart label %4ld not found\n", iap.irs);
            }
            return 0;  // bad retrun
        }
    }

    // this is good for debugging and writes all the auto parameters
    gogoauto_set_function_pointers(iap, &list);
    autlib1_init(&iap, &rap, icp, par);
    autlib_check_dimensions(&iap);

    /* Create the allocations for the global structures used in
       autlib3.c and autlib5.c.  There are purely an efficiency thing.
       The allocation and deallocation of these scratch areas takes
       up a nontrivial amount of time if done directly in the
       wrapper functions in autlib3.c*/
    allocate_global_memory(iap);

    // ----------------------------------------------------------
    //  One-parameter continuations
    // ----------------------------------------------------------

    if (list.type == AUTOAE) {
        autoae(&iap, &rap, par, icp, list.aelist.funi, list.aelist.stpnt, list.aelist.pvli, thl,
               thu, iuz, vuz);
    }
    if (list.type == AUTOBV) {
        autobv(&iap, &rap, par, icp, list.bvlist.funi, list.bvlist.bcni, list.bvlist.icni,
               list.bvlist.stpnt, list.bvlist.pvli, thl, thu, iuz, vuz);
    }

    free(thu);
    fclose(fp3);
    fclose(fp7);
    fclose(fp9);
    return 1;
}

int32
gogoauto_set_function_pointers(iap_type iap, function_list *data) {
    if ((iap.ips == 0 || iap.ips == 1) && ABS(iap.isw) != 2) {
        //	** Algebraic systems.
        if (iap.irs == 0) {
            data->type = AUTOAE;
            data->aelist.funi = funi;
            data->aelist.stpnt = stpnus;
            data->aelist.pvli = pvlsae;
        } else {
            data->type = AUTOAE;
            data->aelist.funi = funi;
            data->aelist.stpnt = stpnae;
            data->aelist.pvli = pvlsae;
        }
    } else if (iap.ips == 11 && ABS(iap.isw) != 2) {
        //	** Waves : Spatially homogeneous solutions,
        if (iap.irs == 0) {
            data->type = AUTOAE;
            data->aelist.funi = fnws;
            data->aelist.stpnt = stpnus;
            data->aelist.pvli = pvlsae;
        } else {
            data->type = AUTOAE;
            data->aelist.funi = fnws;
            data->aelist.stpnt = stpnae;
            data->aelist.pvli = pvlsae;
        }
    } else if (iap.ips == -1 && ABS(iap.isw) != 2) {
        //	** Discrete dynamical systems : fixed points.
        if (iap.irs == 0) {
            data->type = AUTOAE;
            data->aelist.funi = fnds;
            data->aelist.stpnt = stpnus;
            data->aelist.pvli = pvlsae;
        } else {
            data->type = AUTOAE;
            data->aelist.funi = fnds;
            data->aelist.stpnt = stpnae;
            data->aelist.pvli = pvlsae;
        }
    } else if (iap.ips == -2) {
        //	** Time integration.
        if (iap.irs == 0) {
            data->type = AUTOAE;
            data->aelist.funi = fnti;
            data->aelist.stpnt = stpnus;
            data->aelist.pvli = pvlsae;
        } else {
            data->type = AUTOAE;
            data->aelist.funi = fnti;
            data->aelist.stpnt = stpnae;
            data->aelist.pvli = pvlsae;
        }
    } else if (iap.ips == 2 && ABS(iap.isw) != 2) {
        //	** Periodic solutions
        if (iap.itp != 3 && ABS(iap.itp / 10) != 3) {
            if (iap.irs > 0) {
                data->type = AUTOBV;
                data->bvlist.funi = fnps;
                data->bvlist.bcni = bcps;
                data->bvlist.icni = icps;
                data->bvlist.stpnt = stpnbv;
                data->bvlist.pvli = pvlsbv;
            } else {
                data->type = AUTOBV;
                data->bvlist.funi = fnps;
                data->bvlist.bcni = bcps;
                data->bvlist.icni = icps;
                data->bvlist.stpnt = stpnub;
                data->bvlist.pvli = pvlsbv;
            }
        } else {
            data->type = AUTOBV;
            data->bvlist.funi = fnps;
            data->bvlist.bcni = bcps;
            data->bvlist.icni = icps;
            data->bvlist.stpnt = stpnps;
            data->bvlist.pvli = pvlsbv;
        }
    } else if (iap.ips == 12 && ABS(iap.isw) != 2) {
        //	** Wave train solutions to parabolic systems.
        if (iap.itp != 3) {
            if (iap.irs > 0) {
                data->type = AUTOBV;
                data->bvlist.funi = fnwp;
                data->bvlist.bcni = bcps;
                data->bvlist.icni = icps;
                data->bvlist.stpnt = stpnbv;
                data->bvlist.pvli = pvlsbv;
            } else {
                data->type = AUTOBV;
                data->bvlist.funi = fnwp;
                data->bvlist.bcni = bcps;
                data->bvlist.icni = icps;
                data->bvlist.stpnt = stpnub;
                data->bvlist.pvli = pvlsbv;
            }
        } else {
            data->type = AUTOBV;
            data->bvlist.funi = fnwp;
            data->bvlist.bcni = bcps;
            data->bvlist.icni = icps;
            data->bvlist.stpnt = stpnwp;
            data->bvlist.pvli = pvlsbv;
        }
    } else if (iap.ips == 4 && ABS(iap.isw) != 2) {
        //	** Boundary value problems.
        if (iap.irs > 0) {
            data->type = AUTOBV;
            data->bvlist.funi = funi;
            data->bvlist.bcni = bcni;
            data->bvlist.icni = icni;
            data->bvlist.stpnt = stpnbv;
            data->bvlist.pvli = pvlsbv;
        } else {
            data->type = AUTOBV;
            data->bvlist.funi = funi;
            data->bvlist.bcni = bcni;
            data->bvlist.icni = icni;
            data->bvlist.stpnt = stpnub;
            data->bvlist.pvli = pvlsbv;
        }
    } else if (iap.ips == 7 && ABS(iap.isw) != 2) {
        //	** Boundary value problems with Floquet multipliers.
        if (iap.irs > 0) {
            data->type = AUTOBV;
            data->bvlist.funi = funi;
            data->bvlist.bcni = bcni;
            data->bvlist.icni = icni;
            data->bvlist.stpnt = stpnbv;
            data->bvlist.pvli = pvlsbv;
        } else {
            data->type = AUTOBV;
            data->bvlist.funi = funi;
            data->bvlist.bcni = bcni;
            data->bvlist.icni = icni;
            data->bvlist.stpnt = stpnub;
            data->bvlist.pvli = pvlsbv;
        }
    } else if (iap.ips == 9 && ABS(iap.isw) != 2) {
        //	** Homoclinic bifurcation analysis.
        if (iap.irs > 0) {
            data->type = AUTOBV;
            data->bvlist.funi = fnho;
            data->bvlist.bcni = bcho;
            data->bvlist.icni = icho;
            data->bvlist.stpnt = stpnbv;
            data->bvlist.pvli = pvlsho;
        } else {
            data->type = AUTOBV;
            data->bvlist.funi = fnho;
            data->bvlist.bcni = bcho;
            data->bvlist.icni = icho;
            data->bvlist.stpnt = stpnho;
            data->bvlist.pvli = pvlsho;
        }
    } else if (iap.ips == 14) {
        //	** Evolution calculations for parabolic systems.
        //	   (Periodic boundary conditions.)
        if (iap.irs > 0) {
            data->type = AUTOBV;
            data->bvlist.funi = fnpe;
            data->bvlist.bcni = bcps;
            data->bvlist.icni = icpe;
            data->bvlist.stpnt = stpnbv;
            data->bvlist.pvli = pvlsbv;
        } else {
            data->type = AUTOBV;
            data->bvlist.funi = fnpe;
            data->bvlist.bcni = bcps;
            data->bvlist.icni = icpe;
            data->bvlist.stpnt = stpnub;
            data->bvlist.pvli = pvlsbv;
        }
    } else if (iap.ips == 15 && ABS(iap.isw) == 1) {
        //	** Optimization of periodic solutions.
        if (iap.nfpr < 6) {
            data->type = AUTOBV;
            data->bvlist.funi = fnpo;
            data->bvlist.bcni = bcpo;
            data->bvlist.icni = icpo;
            data->bvlist.stpnt = stpnpo;
            data->bvlist.pvli = pvlsbv;
        } else {
            data->type = AUTOBV;
            data->bvlist.funi = fnpo;
            data->bvlist.bcni = bcpo;
            data->bvlist.icni = icpo;
            data->bvlist.stpnt = stpnbv;
            data->bvlist.pvli = pvlsbv;
        }
    } else if (iap.ips == 16) {
        //	** Evolution calculations for parabolic systems.
        //	   (User supplied boundary conditions.)
        if (iap.irs > 0) {
            data->type = AUTOBV;
            data->bvlist.funi = fnpe;
            data->bvlist.bcni = bcni;
            data->bvlist.icni = icpe;
            data->bvlist.stpnt = stpnbv;
            data->bvlist.pvli = pvlsbv;
        } else {
            data->type = AUTOBV;
            data->bvlist.funi = fnpe;
            data->bvlist.bcni = bcni;
            data->bvlist.icni = icpe;
            data->bvlist.stpnt = stpnub;
            data->bvlist.pvli = pvlsbv;
        }
    } else if (iap.ips == 17) {
        //	** Continuation of stationary states of parabolic systems.
        //	   (User supplied boundary conditions.)
        if (iap.irs > 0) {
            data->type = AUTOBV;
            data->bvlist.funi = fnsp;
            data->bvlist.bcni = bcni;
            data->bvlist.icni = icpe;
            data->bvlist.stpnt = stpnbv;
            data->bvlist.pvli = pvlsbv;
        } else {
            data->type = AUTOBV;
            data->bvlist.funi = fnsp;
            data->bvlist.bcni = bcni;
            data->bvlist.icni = icpe;
            data->bvlist.stpnt = stpnub;
            data->bvlist.pvli = pvlsbv;
        }
    } else if (iap.ips == 5) {
        //	** Algebraic optimization problems.
        int32 nfpr = (int32)iap.nfpr;
        if (iap.itp % 10 == 2 || iap.irs == 0) {
            nfpr++;
        }
        if (nfpr == 2) {
            if (iap.irs > 0) {
                data->type = AUTOAE;
                data->aelist.funi = fnc1;
                data->aelist.stpnt = stpnae;
                data->aelist.pvli = pvlsae;
            } else {
                data->type = AUTOAE;
                data->aelist.funi = fnc1;
                data->aelist.stpnt = stpnc1;
                data->aelist.pvli = pvlsae;
            }
        } else {
            if (iap.itp % 10 != 2) {
                data->type = AUTOAE;
                data->aelist.funi = fnc2;
                data->aelist.stpnt = stpnae;
                data->aelist.pvli = pvlsae;
            } else {
                data->type = AUTOAE;
                data->aelist.funi = fnc2;
                data->aelist.stpnt = stpnc2;
                data->aelist.pvli = pvlsae;
            }
        }
    }
    // ---------------------------------------------------
    // ---------------------------------------------------
    //	Two-Parameter Continuation.
    // ---------------------------------------------------
    // ---------------------------------------------------

    else if (iap.ips <= 1 && ABS(iap.isw) == 2 && (iap.itp == 1 || iap.itp == 2)) {
        //	** Fold continuation (algebraic problems).
        data->type = AUTOAE;
        data->aelist.funi = fnlp;
        data->aelist.stpnt = stpnlp;
        data->aelist.pvli = pvlsae;
    } else if (iap.ips <= 1 && ABS(iap.isw) == 2 &&
               (ABS(iap.itp) / 10 == 1 || ABS(iap.itp) / 10 == 2)) {
        //	** Fold continuation (algebraic problems, restart).
        data->type = AUTOAE;
        data->aelist.funi = fnlp;
        data->aelist.stpnt = stpnae;
        data->aelist.pvli = pvlsae;
    } else if ((iap.ips == 0 || iap.ips == 1) && ABS(iap.isw) == 2 && iap.itp == 3) {
        //	** Hopf bifurcation continuation (ODE).
        data->type = AUTOAE;
        data->aelist.funi = fnhb;
        data->aelist.stpnt = stpnhb;
        data->aelist.pvli = pvlsae;
    } else if (ABS(iap.ips) == 1 && ABS(iap.isw) == 2 && ABS(iap.itp) / 10 == 3) {
        //	** Hopf bifurcation continuation (ODE, restart).
        data->type = AUTOAE;
        data->aelist.funi = fnhb;
        data->aelist.stpnt = stpnae;
        data->aelist.pvli = pvlsae;
    } else if (iap.ips == 11 && ABS(iap.isw) == 2 && iap.itp == 3) {
        //	** Hopf bifurcation continuation (Waves).
        data->type = AUTOAE;
        data->aelist.funi = fnhw;
        data->aelist.stpnt = stpnhw;
        data->aelist.pvli = pvlsae;
    } else if (iap.ips == 11 && ABS(iap.isw) == 2 && ABS(iap.itp) / 10 == 3) {
        //	** Hopf bifurcation continuation (Waves, restart).
        data->type = AUTOAE;
        data->aelist.funi = fnhw;
        data->aelist.stpnt = stpnae;
        data->aelist.pvli = pvlsae;
    } else if (iap.ips == -1 && ABS(iap.isw) == 2 && iap.itp == 3) {
        //	** Hopf bifurcation continuation (Maps).
        data->type = AUTOAE;
        data->aelist.funi = fnhd;
        data->aelist.stpnt = stpnhd;
        data->aelist.pvli = pvlsae;
    } else if (iap.ips == -1 && ABS(iap.isw) == 2 && ABS(iap.itp) / 10 == 3) {
        //	** Hopf bifurcation continuation (Maps).
        data->type = AUTOAE;
        data->aelist.funi = fnhd;
        data->aelist.stpnt = stpnae;
        data->aelist.pvli = pvlsae;
    } else if (iap.ips == 2 && ABS(iap.isw) == 2 && (iap.itp == 5 || iap.itp == 6)) {
        //	** Fold continuation (Periodic solutions, start).
        data->type = AUTOBV;
        data->bvlist.funi = fnpl;
        data->bvlist.bcni = bcpl;
        data->bvlist.icni = icpl;
        data->bvlist.stpnt = stpnpl;
        data->bvlist.pvli = pvlsbv;
    } else if (iap.ips == 2 && ABS(iap.isw) == 2 &&
               (ABS(iap.itp) / 10 == 5 || ABS(iap.itp) / 10 == 6)) {
        //        ** Fold continuation (Periodic solutions, restart).
        data->type = AUTOBV;
        data->bvlist.funi = fnpl;
        data->bvlist.bcni = bcpl;
        data->bvlist.icni = icpl;
        data->bvlist.stpnt = stpnbv;
        data->bvlist.pvli = pvlsbv;
    } else if (iap.ips == 2 && ABS(iap.isw) == 2 && iap.itp == 7) {
        //        ** Continuation of period doubling bifurcations (start).
        data->type = AUTOBV;
        data->bvlist.funi = fnpd;
        data->bvlist.bcni = bcpd;
        data->bvlist.icni = icpd;
        data->bvlist.stpnt = stpnpd;
        data->bvlist.pvli = pvlsbv;
    } else if (iap.ips == 2 && ABS(iap.isw) == 2 && ABS(iap.itp) / 10 == 7) {
        //        ** Continuation of period doubling bifurcations (restart).
        data->type = AUTOBV;
        data->bvlist.funi = fnpd;
        data->bvlist.bcni = bcpd;
        data->bvlist.icni = icpd;
        data->bvlist.stpnt = stpnbv;
        data->bvlist.pvli = pvlsbv;
    } else if (iap.ips == 2 && ABS(iap.isw) == 2 && iap.itp == 8) {
        //        ** Continuation of torus bifurcations (start).
        data->type = AUTOBV;
        data->bvlist.funi = fntr;
        data->bvlist.bcni = bctr;
        data->bvlist.icni = ictr;
        data->bvlist.stpnt = stpntr;
        data->bvlist.pvli = pvlsbv;
    } else if (iap.ips == 2 && ABS(iap.isw) == 2 && ABS(iap.itp) / 10 == 8) {
        //        ** Continuation of torus bifurcations (restart).
        data->type = AUTOBV;
        data->bvlist.funi = fntr;
        data->bvlist.bcni = bctr;
        data->bvlist.icni = ictr;
        data->bvlist.stpnt = stpnbv;
        data->bvlist.pvli = pvlsbv;
    } else if (iap.ips == 4 && ABS(iap.isw) == 2 && (iap.itp == 5 || iap.itp == 6)) {
        //        ** Continuation of folds (BVP, start).
        data->type = AUTOBV;
        data->bvlist.funi = fnbl;
        data->bvlist.bcni = bcbl;
        data->bvlist.icni = icbl;
        data->bvlist.stpnt = stpnbl;
        data->bvlist.pvli = pvlsbv;
    } else if (iap.ips == 4 && ABS(iap.isw) == 2 &&
               (ABS(iap.itp) / 10 == 5 || ABS(iap.itp) / 10 == 6)) {
        //        ** Continuation of folds (BVP, restart).
        data->type = AUTOBV;
        data->bvlist.funi = fnbl;
        data->bvlist.bcni = bcbl;
        data->bvlist.icni = icbl;
        data->bvlist.stpnt = stpnbv;
        data->bvlist.pvli = pvlsbv;
    } else {
        // Error in INIT.
        printf("\nInitialization Error CRASH!!\n");
        printf("itp=%ld ips=%ld isw=%ld\n", iap.itp, iap.ips, iap.isw);
    }
    // -----------------------------------------------------------------------

    return 0;
}
