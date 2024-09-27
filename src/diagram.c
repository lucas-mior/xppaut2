#include "functions.h"
#include "autevd.h"

#include <stdbool.h>
#include "auto_nox.h"
#include <stdlib.h>
#include <stdio.h>
#include "autlim.h"
#include "integers.h"

int32 NBifs = 0;
Diagram *bifd;

void
diagram_start(int32 n) {
    NBifs = 1;
    bifd = xmalloc(sizeof(*bifd));
    bifd->prev = NULL;
    bifd->next = NULL;
    bifd->index = 0;
    bifd->uhi = xmalloc((usize)n*sizeof(double));
    bifd->ulo = xmalloc((usize)n*sizeof(double));
    bifd->u0 = xmalloc((usize)n*sizeof(double));
    bifd->ubar = xmalloc((usize)n*sizeof(double));
    bifd->evr = xmalloc((usize)n*sizeof(double));
    bifd->evi = xmalloc((usize)n*sizeof(double));
    bifd->norm = 0;
    bifd->lab = 0;

    DiagFlag = 0;
    return;
}

void
diagram_edit_start(int32 ibr, int32 ntot, int32 itp, int32 lab, int32 nfpar,
                   double a, double *uhi, double *ulo, double *u0, double *ubar,
                   double *par, double per, int32 n, int32 icp1, int32 icp2,
                   int32 icp3, int32 icp4, double *evr, double *evi) {
    edit_diagram(bifd, ibr, ntot, itp, lab, nfpar, a, uhi, ulo, u0, ubar, par,
                 per, n, icp1, icp2, icp3, icp4, AutoTwoParam, evr, evi,
                 blrtn.torper);
    return;
}

void
edit_diagram(Diagram *d, int32 ibr, int32 ntot, int32 itp, int32 lab,
             int32 nfpar, double a, double *uhi, double *ulo, double *u0,
             double *ubar, double *par, double per, int32 n, int32 icp1,
             int32 icp2, int32 icp3, int32 icp4, int32 flag2, double *evr,
             double *evi, double tp) {
    d->calc = TypeOfCalc;
    d->ibr = ibr;
    d->ntot = ntot;
    d->itp = itp;
    d->lab = lab;
    d->nfpar = nfpar;
    d->norm = a;
    for (int32 i = 0; i < 8; i++) {
        d->par[i] = par[i];
    }

    d->per = per;

    d->icp1 = icp1;
    d->icp2 = icp2;
    d->icp3 = icp3;
    d->icp4 = icp4;
    d->flag2 = flag2;
    for (int32 i = 0; i < n; i++) {
        d->ulo[i] = ulo[i];
        d->uhi[i] = uhi[i];
        d->ubar[i] = ubar[i];
        d->u0[i] = u0[i];
        d->evr[i] = evr[i];
        d->evi[i] = evi[i];
    }
    d->torper = tp;
    return;
}

void
add_diagram(int32 ibr, int32 ntot, int32 itp, int32 lab, int32 nfpar, double a,
            double *uhi, double *ulo, double *u0, double *ubar, double *par,
            double per, int32 n, int32 icp1, int32 icp2, int32 icp3, int32 icp4,
            int32 flag2, double *evr, double *evi) {
    Diagram *d;
    Diagram *dnew;

    d = bifd;
    while (d->next != NULL) {
        d = (d->next);
    }
    d->next = xmalloc(sizeof(*(d->next)));
    dnew = d->next;
    dnew->next = NULL;
    dnew->prev = d;
    dnew->uhi = xmalloc((usize)n*sizeof(double));
    dnew->ulo = xmalloc((usize)n*sizeof(double));
    dnew->u0 = xmalloc((usize)n*sizeof(double));
    dnew->ubar = xmalloc((usize)n*sizeof(double));
    dnew->evr = xmalloc((usize)n*sizeof(double));
    dnew->evi = xmalloc((usize)n*sizeof(double));
    dnew->index = NBifs;
    NBifs++;
    edit_diagram(dnew, ibr, ntot, itp, lab, nfpar, a, uhi, ulo, u0, ubar, par,
                 per, n, icp1, icp2, icp3, icp4, flag2, evr, evi, blrtn.torper);
    return;
}

void
kill_diagrams(void) {
    Diagram *d;
    Diagram *dnew;
    d = bifd;
    while (d->next != NULL) { /*  Move to the end of the tree  */
        d = d->next;
    }
    while (d->prev != NULL) {
        dnew = d->prev;
        d->next = NULL;
        d->prev = NULL;
        free(d->uhi);
        free(d->ulo);
        free(d->u0);
        free(d->ubar);
        free(d->evr);
        free(d->evi);
        free(d);
        d = dnew;
    }
    free(bifd->uhi);
    free(bifd->ulo);
    free(bifd->u0);
    free(bifd->ubar);
    free(bifd->evr);
    free(bifd->evi);
    free(bifd);
    diagram_start(NODE);
    return;
}

void
diagram_redraw(void) {
    Diagram *d;
    int32 type;
    int32 flag = 0;
    auto_nox_draw_bix_axes();
    d = bifd;
    if (d->next == NULL) {
        return;
    }
    while (true) {
        type = autevd_get_bif_type(d->ibr, d->ntot);

        if (d->ntot == 1) {
            flag = 0;
        } else {
            flag = 1;
        }
        auto_nox_add_point(d->par, d->per, d->uhi, d->ulo, d->ubar, d->norm,
                           type, flag, d->lab, d->icp1, d->icp2, d->flag2,
                           d->evr, d->evi);
        d = d->next;
        if (d == NULL) {
            break;
        }
    }
    return;
}

void
diagram_write_info_out(void) {
    char filename[XPP_MAX_NAME];
    Diagram *d;
    int32 type;
    /*int32 flag=0
     */
    int32 status;
    int32 icp1;
    int32 icp2;
    double *par;
    double par1, par2 = 0, *uhigh, *ulow, per;
    FILE *fp;
    sprintf(filename, "allinfo.dat");
    /* status=dialog_box_get("Write all
     * info","Filename",filename,"Ok","Cancel",60);
     */
    status = init_conds_file_selector("Write all info", filename, "*.dat");

    if (status == 0) {
        return;
    }
    fp = fopen(filename, "w");
    if (fp == NULL) {
        ggets_err_msg("Can't open file");
        return;
    }

    d = bifd;
    if (d->next == NULL) {
        return;
    }
    while (true) {
        type = autevd_get_bif_type(d->ibr, d->ntot);

        icp1 = d->icp1;
        icp2 = d->icp2;
        par = d->par;
        per = d->per;
        uhigh = d->uhi;
        ulow = d->ulo;
        par1 = par[icp1];
        if (icp2 < NAutoPar) {
            par2 = par[icp2];
        } else {
            par2 = par1;
        }

        fprintf(fp, "%d %d %d %g %g %g ", type, d->ibr, d->flag2, par1, par2,
                per);
        for (int32 i = 0; i < NODE; i++) {
            fprintf(fp, "%g ", uhigh[i]);
        }
        for (int32 i = 0; i < NODE; i++) {
            fprintf(fp, "%g ", ulow[i]);
        }
        for (int32 i = 0; i < NODE; i++) {
            fprintf(fp, "%g %g ", d->evr[i], d->evi[i]);
        }
        fprintf(fp, "\n");
        d = d->next;
        if (d == NULL) {
            break;
        }
    }
    fclose(fp);
    return;
}

void
diagram_load_browser_with_branch(int32 ibr, int32 pts, int32 pte) {
    Diagram *d;
    int32 i;
    int32 j;
    int32 pt;
    int32 icp1;
    double *par;
    double par1;
    double *u0;
    int32 first;
    int32 last;
    int32 nrows;
    first = ABS(pts);
    last = ABS(pte);
    if (last <
        first) { /* reorder the points so that we will store in right range*/
        i = first;
        first = last;
        last = i;
    }
    nrows = last - first + 1;
    d = bifd;
    if (d->next == NULL) {
        return;
    }
    j = 0;
    while (true) {
        pt = ABS(d->ntot);
        if ((d->ibr == ibr) && (pt >= first) && (pt <= last)) {
            icp1 = d->icp1;
            par = d->par;
            u0 = d->u0;

            par1 = par[icp1];
            storage[0][j] = par1;
            for (i = 0; i < NODE; i++) {
                storage[i + 1][j] = u0[i];
            }
            j++;
        }
        d = d->next;
        if (d == NULL) {
            break;
        }
    }
    storind = nrows;
    refresh_browser(nrows);
}

void
diagram_write_init_data_file(void) {
    char filename[XPP_MAX_NAME];
    Diagram *d;
    int32 status;
    int32 icp1;
    double *par;
    double par1;
    double *u0;
    FILE *fp;
    sprintf(filename, "initdata.dat");
    /* status=dialog_box_get("Write all
     * info","Filename",filename,"Ok","Cancel",60);
     */
    status =
        init_conds_file_selector("Write init data file", filename, "*.dat");

    if (status == 0) {
        return;
    }
    fp = fopen(filename, "w");
    if (fp == NULL) {
        ggets_err_msg("Can't open file");
        return;
    }

    d = bifd;
    if (d->next == NULL) {
        return;
    }
    while (true) {
        icp1 = d->icp1;
        par = d->par;
        /*
        uhigh=d->uhi;
        ulow=d->ulo;
        ubar=d->ubar;
        Unused here??
        */
        u0 = d->u0;

        /*
        a=d->norm;

        Unused here??
        */
        par1 = par[icp1];

        /* fprintf(fp,"%d %d %g %g %g ",
           type,d->ibr,par1,par2,per); */
        fprintf(fp, "%g ", par1);
        for (int32 i = 0; i < NODE; i++) {
            fprintf(fp, "%g ", u0[i]);
        }
        fprintf(fp, "\n");
        d = d->next;
        if (d == NULL) {
            break;
        }
    }
    fclose(fp);
    return;
}

void
diagram_write_pts(void) {
    char filename[XPP_MAX_NAME];
    Diagram *d;
    int32 type;
    int32 status;
    int32 icp1;
    int32 icp2;
    double *par;
    double x, y1, y2, par1, par2 = 0, a, *uhigh, *ulow, *ubar, per;
    FILE *fp;
    sprintf(filename, "diagram.dat");
    status = init_conds_file_selector("Write points", filename, "*.dat");
    if (status == 0) {
        return;
    }
    fp = fopen(filename, "w");
    if (fp == NULL) {
        ggets_err_msg("Can't open file");
        return;
    }

    d = bifd;
    if (d->next == NULL) {
        return;
    }
    while (true) {
        type = autevd_get_bif_type(d->ibr, d->ntot);

        icp1 = d->icp1;
        icp2 = d->icp2;
        par = d->par;
        per = d->per;
        uhigh = d->uhi;
        ulow = d->ulo;
        ubar = d->ubar;
        a = d->norm;
        par1 = par[icp1];
        if (icp2 < NAutoPar) {
            par2 = par[icp2];
        }

        /* now we have to check is the diagram parameters correspond to the
           current view
        */
        if (auto_nox_check_plot_type(d->flag2, icp1, icp2) == 1) {
            auto_nox_xy_plot(&x, &y1, &y2, par1, par2, per, uhigh, ulow, ubar,
                             a);
            fprintf(fp, "%g %g %g %d %d %d\n", x, y1, y2, type, ABS(d->ibr),
                    d->flag2);
        }
        d = d->next;
        if (d == NULL) {
            break;
        }
    }
    fclose(fp);
    return;
}

void
diagram_post_auto(void) {
    char filename[XPP_MAX_NAME];
    Diagram *d;
    int32 type;
    int32 flag = 0;
    int32 status;
    sprintf(filename, "auto.ps");
    status = init_conds_file_selector("Postscript", filename, "*.ps");
    if (status == 0) {
        return;
    }
    if (!ps_init(filename, PS_Color)) {
        return;
    }
    auto_nox_draw_ps_axes();
    d = bifd;
    if (d->next == NULL) {
        return;
    }
    while (true) {
        type = autevd_get_bif_type(d->ibr, d->ntot);
        if (type < 0) {
            ggets_plintf("Unable to get bifurcation type.\n");
        }
        if (d->ntot == 1) {
            flag = 0;
        } else {
            flag = 1;
        }
        auto_nox_add_ps_point(d->par, d->per, d->uhi, d->ulo, d->ubar, d->norm,
                              type, flag, d->icp1, d->icp2, d->flag2);
        d = d->next;
        if (d == NULL) {
            break;
        }
    }
    ps_end();
    graphics_set_normal_scale();
    return;
}

void
diagram_svg_auto(void) {
    char filename[XPP_MAX_NAME];
    Diagram *d;
    int32 type;
    int32 flag = 0;
    int32 status;
    sprintf(filename, "auto.svg");
    status = init_conds_file_selector("SVG", filename, "*.svg");
    if (status == 0) {
        return;
    }
    if (!svg_init(filename)) {
        return;
    }
    auto_nox_draw_svg_axes();
    d = bifd;
    if (d->next == NULL) {
        return;
    }
    while (true) {
        type = autevd_get_bif_type(d->ibr, d->ntot);
        if (type < 0) {
            ggets_plintf("Unable to get bifurcation type.\n");
        }
        if (d->ntot == 1) {
            flag = 0;
        } else {
            flag = 1;
        }
        auto_nox_add_ps_point(d->par, d->per, d->uhi, d->ulo, d->ubar, d->norm,
                              type, flag, d->icp1, d->icp2, d->flag2);
        d = d->next;
        if (d == NULL) {
            break;
        }
    }
    svg_end();

    graphics_set_normal_scale();
    return;
}

void
diagram_bound(double *xlo, double *xhi, double *ylo, double *yhi) {
    Diagram *d;
    int32 type;

    double x;
    double y1;
    double y2;
    double par1;
    double par2 = 0.0;
    d = bifd;
    if (d->next == NULL) {
        return;
    }
    *xlo = 1.e16;
    *ylo = *xlo;
    *xhi = -*xlo;
    *yhi = -*ylo;
    while (true) {
        type = autevd_get_bif_type(d->ibr, d->ntot);
        if (type < 1) {
            ggets_plintf("Unable to get bifurcation type.\n");
        }
        par1 = d->par[d->icp1];
        if (d->icp2 < NAutoPar) {
            par2 = d->par[d->icp2];
        }
        auto_nox_xy_plot(&x, &y1, &y2, par1, par2, d->per, d->uhi, d->ulo,
                         d->ubar, d->norm);
        if (x < *xlo) {
            *xlo = x;
        }
        if (x > *xhi) {
            *xhi = x;
        }
        if (y2 < *ylo) {
            *ylo = y2;
        }
        if (y1 > *yhi) {
            *yhi = y1;
        }
        d = d->next;
        if (d == NULL) {
            break;
        }
    }
    return;
}

int32
diagram_save(FILE *fp, int32 n) {
    Diagram *d;
    fprintf(fp, "%d\n", NBifs - 1);
    if (NBifs == 1) {
        return -1;
    }
    d = bifd;
    while (true) {
        fprintf(fp, "%d %d %d %d %d %d %d %d %d %d %d %d\n", d->calc, d->ibr,
                d->ntot, d->itp, d->lab, d->index, d->nfpar, d->icp1, d->icp2,
                d->icp3, d->icp4, d->flag2);
        for (int32 i = 0; i < 8; i++) {
            fprintf(fp, "%g ", d->par[i]);
        }
        fprintf(fp, "%g %g \n", d->norm, d->per);

        for (int32 i = 0; i < n; i++) {
            fprintf(fp, "%f %f %f %f %f %f\n", d->u0[i], d->uhi[i], d->ulo[i],
                    d->ubar[i], d->evr[i], d->evi[i]);
        }
        d = d->next;
        if (d == NULL) {
            break;
        }
    }
    return 1;
}

int32
diagram_load(FILE *fp, int32 node) {
    double u0[NAUTO];
    double uhi[NAUTO];
    double ulo[NAUTO];
    double ubar[NAUTO];
    double evr[NAUTO];
    double evi[NAUTO];
    double norm;
    double par[8];
    double per;

    int32 flag = 0;
    int32 n;
    int32 calc, ibr, ntot, itp, lab, index, nfpar, icp1, icp2, icp3, icp4,
        flag2;
    fscanf(fp, "%d", &n);
    if (n == 0) {
        return -1;
    }

    while (true) {
        fscanf(fp, "%d %d %d %d %d %d %d %d %d %d %d %d", &calc, &ibr, &ntot,
               &itp, &lab, &index, &nfpar, &icp1, &icp2, &icp3, &icp4, &flag2);
        for (int32 i = 0; i < 8; i++) {
            fscanf(fp, "%lg ", &par[i]);
        }
        fscanf(fp, "%lg %lg ", &norm, &per);
        for (int32 i = 0; i < node; i++) {
            fscanf(fp, "%lg %lg %lg %lg %lg %lg", &u0[i], &uhi[i], &ulo[i],
                   &ubar[i], &evr[i], &evi[i]);
        }
        if (flag == 0) {
            diagram_edit_start(ibr, ntot, itp, lab, nfpar, norm, uhi, ulo, u0,
                               ubar, par, per, node, icp1, icp2, icp3, icp4,
                               evr, evi);
            flag = 1;
            DiagFlag = 1;
        } else {
            add_diagram(ibr, ntot, itp, lab, nfpar, norm, uhi, ulo, u0, ubar,
                        par, per, node, icp1, icp2, icp3, icp4, flag2, evr,
                        evi);
        }
        if (index >= n) {
            break;
        }
    }
    return 1;
}
