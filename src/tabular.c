#include "integers.h"
#include "functions.h"
#include "xmalloc.h"
#include "read_dir.h"
#include <stdbool.h>

#include "parserslow.h"

#include <stdlib.h>
#include <string.h>

/*********************************************************
     This is code for read-in tables in XPP
     This should probably be accessible from within the program
     as well.  It will probably be added to the Numerics Menu

     The files consist of y-values of a function evaluated at
     equally spaced points as well as some header information.
     They are ascii files of the form:

     npts <-- Integer
     xlo <--- fp
     xhi <--  fp
     y1  <-- fp
     y2
     ...
     yn

     Thus   dx = (xhi-xlo)/(npts-1)

    If the first line of the file says "xyvals" then the table is of the
    form: x1 < x2 < .... < xn
    npts
    x1 y1
    x2 y2
     ...
    xn yn

  In the creation of the file, one can instead use the following:

 table <name> % numpts xlo xhi formula

 to create a "formula" table which is linearly interpolated

 table <name> @ filename creates an array for two-valued
                functions

 filename has the following info:
 nxpts
 nypts
 xlo
 xhi
 ylo
 yhi

 nx*ny points as follows

 f(x1,y1), f(x2,y1),....,f(xn,y1),
 ...
 f(x1,ym), ..., f(xn,ym)

to be added later
**************************************************************/

#include <stdio.h>

Tabular my_table[MAX_TAB];

static int32 tabular_eval_fun(int32 n, double xlo, double xhi, char *formula,
                              double *y);
static double tabular_interp(double xlo, double h, double x, double *y,
                             int32 i);
static double tabular_lookup_xy(double x, int32 n, double *xv, double *yv);

void
tabular_set_auto_eval_flags(int32 f) {
    for (int32 i = 0; i < MAX_TAB; i++) {
        my_table[i].autoeval = f;
    }
    return;
}

void
tabular_set_table_name(char *name, int32 index) {
    strcpy(my_table[index].name, name);
    return;
}

void
tabular_new_lookup_com(int32 i) {
    char file[128];
    int32 index;
    int32 ok;
    int32 status;
    double xlo;
    double xhi;
    int32 npts;
    char newform[80];

    index = many_pops_select_table();
    if (index == -1) {
        return;
    }
    if (i == 1) {
        // tabular view
        int32  n = my_table[index].n;
        int32  len;
        double *y = my_table[index].y;
        double dx = my_table[index].dx;
        xlo = my_table[index].xlo;
        len = n;
        if (len >= max_stor) {
            len = max_stor - 1;
        }
        for (int32 i2 = 0; i2 < len; i2++) {
            storage[0][i2] = xlo + i2*dx;
            storage[1][i2] = y[i2];
        }
        browser_refresh(len);
        return;
    }
    if (my_table[index].flag == 1) {
        strcpy(file, my_table[index].filename);
        status = init_conds_file_selector("Load table", file, "*.tab");
        if (status == 0) {
            return;
        }
        ok = tabular_load_table(file, index);
        if (ok == 1) {
            strcpy(my_table[index].filename, file);
        }
    }
    if (my_table[index].flag == 2) {
        npts = my_table[index].n;

        xlo = my_table[index].xlo;
        xhi = my_table[index].xhi;
        strcpy(newform, my_table[index].filename);
        ggets_new_int("Auto-evaluate? (1/0)", &my_table[index].autoeval);
        ggets_new_int("NPts: ", &npts);
        ggets_new_float("Xlo: ", &xlo);
        ggets_new_float("Xhi: ", &xhi);
        ggets_new_string("Formula :", newform);
        tabular_create_fun(npts, xlo, xhi, newform, index);
    }
    return;
}

double
tabular_lookup_xy(double x, int32 n, double *xv, double *yv) {
    double dx;
    double dy;
    double x1;
    double y1;
    double x2;
    double y2;

    if (x <= xv[0]) {
        return yv[0] + (yv[1] - yv[0])*(x - xv[0]) / (xv[1] - xv[0]);
    }
    if (x >= xv[n - 1]) {
        return (yv[n - 1] + (yv[n - 2] - yv[n - 1])*(x - xv[n - 1]) /
                                (xv[n - 1] - xv[n - 2]));
    }
    x1 = xv[0];
    y1 = yv[0];
    for (int32 i = 1; i < n; i++) {
        if (x <= xv[i]) {
            x2 = xv[i];
            y2 = yv[i];
            dx = x2 - x1;
            dy = y2 - y1;
            return y1 + dy*(x - x1) / dx;
        }
        x1 = xv[i];
        y1 = yv[i];
    }
    return yv[n - 1];
}

double
tabular_interp(double xlo, double h, double x, double *y, int32 i) {
    double a;
    double b;
    double c;
    double d;
    double ym;
    double y0;
    double y1;
    double y2;
    double tt;
    ym = y[i - 1];
    y0 = y[i];
    y1 = y[i + 1];
    y2 = y[i + 2];
    d = y0;
    b = .5*(y1 + ym - 2*y0);
    a = (3*(y0 - y1) + y2 - ym) / 6;
    c = (6*y1 - y2 - 3*y0 - 2*ym) / 6;
    tt = (x - xlo) / h - i;
    return d + tt*(c + tt*(b + tt*a));
}

double
tabular_lookup(double x, int32 index) {
    double xlo = my_table[index].xlo, xhi = my_table[index].xhi,
           dx = my_table[index].dx;
    double *y;
    double x1;
    double y1;
    double y2;
    int32 i1;
    int32 i2;
    int32 n = my_table[index].n;
    y = my_table[index].y;

    if (my_table[index].flag == 0) {
        return 0.0;  // Not defined
    }
    if (my_table[index].xyvals == 1) {
        return tabular_lookup_xy(x, n, my_table[index].x, y);
    }

    i1 = (int32)((x - xlo) / dx);  // (int32)floor(x) instead of (int32)x ???
    if (my_table[index].interp == 2 && i1 > 0 && i1 < (n - 2)) {
        // if it is on the edge - use linear
        return tabular_interp(xlo, dx, x, y, i1);
    }
    i2 = i1 + 1;
    if (i1 > -1 && i2 < n) {
        x1 = dx*i1 + xlo;
        y1 = y[i1];
        y2 = y[i2];
        if (my_table[index].interp == 0 || my_table[index].interp == 2) {
            return y1 + (y2 - y1)*(x - x1) / dx;
        } else {
#ifdef DEBUG
            ggets_plintf(
                "index=%d; x=%lg; i1=%d; i2=%d; x1=%lg; y1=%lg; y2=%lg\n",
                index, x, i1, i2, x1, y1, y2);
#endif
            return y1;
        }
    }
    if (i1 < 0) {
        return y[0] + (y[1] - y[0])*(x - xlo) / dx;
    }
    if (i2 >= n) {
        return y[n - 1] + (y[n - 1] - y[n - 2])*(x - xhi) / dx;
    }

    return 0.0;
}

void
tabular_init_table(void) {
    for (int32 i = 0; i < MAX_TAB; i++) {
        my_table[i].flag = 0;
        my_table[i].autoeval = 1;
        my_table[i].interp = 0;
    }
    return;
}

void
tabular_redo_all_fun_tables(void) {
    for (int32 i = 0; i < NTable; i++) {
        if (my_table[i].flag == 2 && my_table[i].autoeval == 1) {
            tabular_eval_fun(my_table[i].n, my_table[i].xlo, my_table[i].xhi,
                             my_table[i].filename, my_table[i].y);
        }
    }
    simplenet_update_all_ffts();
    return;
}

int32
tabular_eval_fun(int32 n, double xlo, double xhi, char *formula, double *y) {
    int32 i;

    double dx;
    double oldt;
    int32 command[200];
    int32 ncold = NCON;
    int32 nsym = NSYM;
    if (parserslow_add_expr(formula, command, &i)) {
        ggets_err_msg("Illegal formula...");
        NCON = ncold;
        NSYM = nsym;
        return 0;
    }
    oldt = get_ivar(0);
    dx = (xhi - xlo) / ((double)(n - 1));
    for (i = 0; i < n; i++) {
        set_ivar(0, dx*i + xlo);
        y[i] = evaluate(command);
    }
    set_ivar(0, oldt);
    NCON = ncold;
    NSYM = nsym;
    return 1;
}

int32
tabular_create_fun(int32 npts, double xlo, double xhi, char *formula,
                   int32 index) {
    int32 length = npts;

    if (my_table[index].flag == 1) {
        ggets_err_msg("Not a function table...");
        return 0;
    }
    if (xlo > xhi) {
        ggets_err_msg("Xlo > Xhi ???");
        return 0;
    }
    if (npts < 2) {
        ggets_err_msg("Too few points...");
        return 0;
    }
    if (my_table[index].flag == 0) {
        my_table[index].y =
            xmalloc((usize)length*sizeof(*(my_table[index].y)));
    } else {
        my_table[index].y = (double *)realloc(
            my_table[index].y, (usize)length*sizeof(*(my_table[index].y)));
    }
    if (my_table[index].y == NULL) {
        ggets_err_msg("Unable to allocate table");
        return 0;
    }
    my_table[index].flag = 2;
    if (tabular_eval_fun(npts, xlo, xhi, formula, my_table[index].y)) {
        my_table[index].xlo = xlo;
        my_table[index].xhi = xhi;
        my_table[index].n = npts;
        my_table[index].dx = (xhi - xlo) / ((double)(npts - 1));
        strcpy(my_table[index].filename, formula);
        return 1;
    }
    return 0;
}

int32
tabular_load_table(char *filename, int32 index) {
    char bobtab[100];
    char *bob;
    int32 length;
    double xlo;
    double xhi;
    FILE *fp;
    char filename2[512];
    char ch;
    char error[sizeof(filename2) + sizeof(cur_dir) + 20];
    int32 n = (int32)strlen(filename);
    int32 j = 0;
    int32 flag = 0;
    for (int32 i = 0; i < n; i++) {
        ch = filename[i];
        if ((ch == '"') && flag == 1) {
            break;
        }
        if ((ch == '"') && (flag == 0)) {
            flag = 1;
        }
        if (ch != '"') {
            filename2[j] = ch;
            j++;
        }
    }
    filename2[j] = 0;

    bob = bobtab;

    if (my_table[index].flag == 2) {
        ggets_err_msg("Not a file table...");
        return 0;
    }
    fp = fopen(filename2, "r");
    if (fp == NULL) {
        read_dir_get_directory(cur_dir);
        snprintf(error, sizeof(error), "File<%s> not found in %s", filename2,
                 cur_dir);
        ggets_err_msg(error);
        return 0;
    }
    my_table[index].interp = 0;
    fgets(bob, 100, fp);
    if (bob[0] == 'i')  // closest step value
    {
        my_table[index].interp = 1;
        bob++;  // skip past initial "i" to length
    }
    if (bob[0] == 's')  // cubic spline
    {
        my_table[index].interp = 2;
        bob++;  // skip past initial "i" to length
    }
    length = atoi(bob);
    if (length < 2) {
        ggets_err_msg("Length too small");
        fclose(fp);
        return 0;
    }
    fgets(bob, 100, fp);
    xlo = atof(bob);
    fgets(bob, 100, fp);
    xhi = atof(bob);
    if (xlo >= xhi) {
        ggets_err_msg("xlo >= xhi ??? ");
        fclose(fp);
        return 0;
    }
    if (my_table[index].flag == 0) {
        my_table[index].y =
            xmalloc((usize)length*sizeof(*(my_table[index].y)));
        if (my_table[index].y == NULL) {
            ggets_err_msg("Unable to allocate table");
            fclose(fp);
            return 0;
        }
        for (int32 i = 0; i < length; i++) {
            fgets(bob, 100, fp);
            my_table[index].y[i] = atof(bob);
        }
        my_table[index].xlo = xlo;
        my_table[index].xhi = xhi;
        my_table[index].n = length;
        my_table[index].dx = (xhi - xlo) / (length - 1);
        my_table[index].flag = 1;
        strcpy(my_table[index].filename, filename2);
        fclose(fp);
        return 1;
    }
    my_table[index].y =
        (double *)realloc((void *)my_table[index].y,
                          (usize)length*sizeof(*(my_table[index].y)));
    if (my_table[index].y == NULL) {
        ggets_err_msg("Unable to reallocate table");
        fclose(fp);
        return 0;
    }
    for (int32 i = 0; i < length; i++) {
        fgets(bob, 100, fp);
        my_table[index].y[i] = atof(bob);
    }
    my_table[index].xlo = xlo;
    my_table[index].xhi = xhi;
    my_table[index].n = length;
    my_table[index].dx = (xhi - xlo) / (length - 1);
    my_table[index].flag = 1;
    fclose(fp);
    return 1;
}

int32
tabular_get_lookup_len(int32 i) {
    return my_table[i].n;
}

/*   network stuff

table name <type> ... arguments ...
           conv   npts  weight variable_name klo khi end_cond
           sparse npts  variable filename

      name(0 ... npts-1)

conv:
        name(i) = simplenet_sum(k=klo,khi) weight(k-klo)*variable(i+k)
        with end_cond = zero means skip if off end
                      = periodic means wrap around

sparse:
       need a file with the structure:
       ncon i_1 w_1 ... i_ncon w_ncon
for npts lines
 name(i) = simplenet_sum(j=1,ncon_i) w_j name(i_j)

*/
