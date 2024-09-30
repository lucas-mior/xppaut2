#include "functions.h"
#include "xmalloc.h"
#include "integers.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "xpplim.h"

#include "parserslow.h"

int32 spec_col = 1;
int32 spec_wid = 512;
int32 spec_win = 2;
int32 spec_col2 = 1;
static int32 spec_type = 0;
/* type =0 for PSD
 * type =1 for crossspectrum
 * type =2 for coherence

*/

int32 post_process = 0;

HistInfo hist_inf = {100, 100, 0, 1, 1, 0, 0, 1, 0, 1, ""};

static int32 hist_len;
static int32 four_len;
static double *my_hist[MAX_ODE + 1];
static double *my_four[MAX_ODE + 1];
static int32 hist_here;
int32 four_here;

static int32 histogram_two_d2(void);

static void histogram_just_fourier(int32 flag);
static void histogram_mycor2(double *x, double *y, int32 n, int32 nbins, double *z, int32 flag);
static int32 histogram_spectrum(double *data, int32 nr, int32 win, int32 w_type, double *pow);
static int32 histogram_get_col_info(int32 *col, char *prompt);
static void histogram_new_four(int32 nmodes, int32 col);
static void histogram_fft(double *data, double *ct, double *st, int32 nmodes, int32 length);

int32
histogram_two_d(int32 col1, int32 col2, int32 ndat, int32 n1, int32 n2, double xlo, double xhi,
                double ylo, double yhi) {
    /*
      col1,2 are the data you want to histogram
      ndat - number of points in the data
      n1,2 number of bins for two data streams
      xlo,xhi - range of first column
      ylo,yhi - range of second column
      val[0] = value of first data
      val[1] = value of second data
      val[3] = number of points - which will be normalized by ndat
    EXAMPLE of binning
      if xl0 = 0 and xhi=1 and nbin=10
      dx=1/10
      then bins are [0,1/10), [1/10,2/10), ....,[9/10,1)
      thus  bin j = int32 ((x-xlo)/dx)
            bin k = int32 ((y-ylo)/dy)
            if j<0 or j>=nxbin then skip etc
    */
    int32 i;
    int32 j;
    double dx;
    double dy;
    double norm;
    double x;
    double y;
    dx = (xhi - xlo) / (double)n1;
    dy = (yhi - ylo) / (double)n2;
    norm = 1. / (double)ndat;
    /* now fill the data with the bin values - take the midpoints of
       each bin
    */
    for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
            my_hist[0][i + j*n1] = xlo + (i + .5)*dx;
            my_hist[1][i + j*n1] = ylo + (j + .5)*dy;
            my_hist[2][i + j*n1] = 0.0;
        }
    }
    for (int32 k = 0; k < ndat; k++) {
        x = (storage[col1][k] - xlo) / dx;
        y = (storage[col2][k] - ylo) / dy;
        i = (int32)x;
        j = (int32)y;
        if ((i >= 0) && (i < n1) && (j >= 0) && (j < n2)) {
            my_hist[2][i + j*n1] += norm;
        }
    }
    return 0;
}

void
histogram_back(void) {
    if (hist_here) {
        browser_set_data(my_hist, 1);
        /*
        browser_my.data=my_hist;
        browser_my.col0=1; */
        browser_refresh(hist_len);
    }
    return;
}

void
histogram_new_four(int32 nmodes, int32 col) {
    int32 length = nmodes + 1;
    double total = storage[0][storind - 1] - storage[0][0];
    double *bob;
    if (four_here) {
        adjoints_data_back();
        free(my_four[0]);
        free(my_four[1]);
        free(my_four[2]);
        four_here = 0;
    }
    four_len = nmodes;
    my_four[0] = xmalloc(sizeof(*(my_four[0]))*(usize)length);
    my_four[1] = xmalloc(sizeof(*(my_four[1]))*(usize)length);
    my_four[2] = xmalloc(sizeof(*(my_four[2]))*(usize)length);
    if (my_four[2] == NULL) {
        free(my_four[1]);
        free(my_four[2]);
        ggets_err_msg("Cant allocate enough memory...");
        return;
    }
    four_here = 1;
    for (int32 i = 3; i <= NEQ; i++) {
        my_four[i] = storage[i];
    }
    for (int32 i = 0; i < length; i++) {
        my_four[0][i] = (double)i / total;
    }
    bob = browser_get_data_col(col);
    histogram_fft(bob, my_four[1], my_four[2], nmodes, storind);
    if (four_here) {
        // histogram four back
        browser_set_data(my_four, 1);
        browser_refresh(four_len);
    }
    ggets_ping();
    return;
}

void
histogram_post_process_stuff(void) {
    if (post_process == 0) {
        return;
    }
    if (N_plist < 1) {
        plotlist = xmalloc(sizeof(*plotlist)*10);
    }
    N_plist = 2;
    plotlist[0] = 0;
    plotlist[1] = 1;
    if (post_process == 7) {  // two-d histogram stuff
        histogram_two_d2();
        return;
    }
    if (post_process == 1) {
        histogram_new(hist_inf.nbins, hist_inf.xlo, hist_inf.xhi, hist_inf.col, 0, "", 0);
        return;
    }
    if (post_process == 2) {
        histogram_just_fourier(0);
        return;
    }
    if (post_process == 3) {
        histogram_just_fourier(1);
        return;
    }
    if (post_process > 3 && post_process < 7) {
        // histogram just sd
        int32 flag = post_process - 4;
        int32 length;
        double total = storage[0][storind - 1] - storage[0][0];
        spec_type = flag;
        if (hist_here) {
            adjoints_data_back();
            free(my_hist[0]);
            free(my_hist[1]);
            if (hist_here == 2) {
                free(my_hist[2]);
            }
            hist_here = 0;
        }
        hist_len = spec_wid / 2;
        length = hist_len + 2;
        my_hist[0] = xmalloc(sizeof(*(my_hist[0]))*(usize)length);
        my_hist[1] = xmalloc(sizeof(*(my_hist[1]))*(usize)length);
        if (my_hist[1] == NULL) {
            free(my_hist[0]);
            ggets_err_msg("Cannot allocate enough...");
            return;
        }
        hist_here = 1;
        for (int32 i = 2; i <= NEQ; i++) {
            my_hist[i] = storage[i];
        }
        for (int32 j = 0; j < hist_len; j++) {
            my_hist[0][j] = ((double)j*storind / spec_wid) / total;
        }
        if (spec_type == 0) {
            histogram_spectrum(storage[spec_col], storind, spec_wid, spec_win, my_hist[1]);
        } else {
            histogram_cross_spectrum(storage[spec_col], storage[spec_col2], storind, spec_wid,
                                     spec_win, my_hist[1], spec_type);
        }
        histogram_back();
        ggets_ping();
        return;
    }
    return;
}

int32
histogram_two_d2(void) {
    int32 length;
    length = hist_inf.nbins*hist_inf.nbins2;
    if (length >= max_stor) {
        length = max_stor - 1;
    }

    if (hist_here) {
        adjoints_data_back();
        free(my_hist[0]);
        free(my_hist[1]);
        if (hist_here == 2) {
            free(my_hist[2]);
        }
        hist_here = 0;
    }

    hist_len = length;
    my_hist[0] = xmalloc(sizeof(*(my_hist[0]))*(usize)length);
    my_hist[1] = xmalloc(sizeof(*(my_hist[1]))*(usize)length);
    my_hist[2] = xmalloc(sizeof(*(my_hist[2]))*(usize)length);
    if (my_hist[2] == NULL) {
        free(my_hist[0]);
        free(my_hist[1]);
        ggets_err_msg("Cannot allocate enough...");
        return -1;
    }
    hist_here = 2;
    for (int32 i = 3; i <= NEQ; i++) {
        my_hist[i] = storage[i];
    }
    hist_len = length;
    histogram_two_d(hist_inf.col, hist_inf.col2, storind, hist_inf.nbins, hist_inf.nbins2,
                    hist_inf.xlo, hist_inf.xhi, hist_inf.ylo, hist_inf.yhi);

    histogram_back();

    ggets_ping();

    return 1;
}

int32
histogram_new_2d(void) {
    if ((NEQ < 2) || (storind < 3)) {
        ggets_err_msg("Need more data and at least 3 columns");
        return 0;
    }
    if (histogram_get_col_info(&hist_inf.col, "Variable 1 ") == 0) {
        return -1;
    }
    ggets_new_int("Number of bins ", &hist_inf.nbins);
    ggets_new_float("Low ", &hist_inf.xlo);
    ggets_new_float("Hi ", &hist_inf.xhi);
    if (hist_inf.nbins < 2) {
        ggets_err_msg("At least 2 bins\n");
        return 0;
    }
    if (hist_inf.xlo >= hist_inf.xhi) {
        ggets_err_msg("Low must be less than hi");
        return 0;
    }

    if (histogram_get_col_info(&hist_inf.col2, "Variable 2 ") == 0) {
        return -1;
    }
    ggets_new_int("Number of bins ", &hist_inf.nbins2);
    ggets_new_float("Low ", &hist_inf.ylo);
    ggets_new_float("Hi ", &hist_inf.yhi);

    if (hist_inf.nbins2 < 2) {
        ggets_err_msg("At least 2 bins\n");
        return 0;
    }
    if (hist_inf.ylo >= hist_inf.yhi) {
        ggets_err_msg("Low must be less than hi");
        return 0;
    }

    return histogram_two_d2();
}

void
histogram_new(int32 nbins, double zlo, double zhi, int32 col, int32 col2, char *condition,
              int32 which) {
    int32 i;
    int32 index;
    int32 command[256];
    int32 cond = 0;
    int32 flag = 1;
    double z;
    double y;
    double dz;
    int32 length = nbins + 1;

    if (length >= max_stor) {
        length = max_stor - 1;
    }
    dz = (zhi - zlo) / (double)(length - 1);
    if (hist_here) {
        adjoints_data_back();
        free(my_hist[0]);
        free(my_hist[1]);
        if (hist_here == 2) {
            free(my_hist[2]);
        }
        hist_here = 0;
    }
    hist_len = length;
    my_hist[0] = xmalloc(sizeof(*(my_hist[0]))*(usize)length);
    my_hist[1] = xmalloc(sizeof(*(my_hist[1]))*(usize)length);
    if (my_hist[1] == NULL) {
        free(my_hist[0]);
        ggets_err_msg("Cannot allocate enough...");
        return;
    }
    hist_here = 1;
    for (i = 2; i <= NEQ; i++) {
        my_hist[i] = storage[i];
    }
    for (i = 0; i < length; i++) {
        my_hist[0][i] = (double)(zlo + dz*i);
        my_hist[1][i] = 0.0;
    }
    if (which == 0) {
        if (strlen(condition) == 0) {
            cond = 0;
        } else {
            if (parserslow_add_expr(condition, command, &i)) {
                ggets_err_msg("Bad condition. Ignoring...");

            } else {
                cond = 1;
            }
        }
        for (i = 0; i < storind; i++) {
            flag = 1;
            if (cond) {
                for (int32 j = 0; j < NODE + 1; j++) {
                    set_ivar(j, (double)storage[j][i]);
                }
                for (int32 j = 0; j < NMarkov; j++) {
                    set_ivar(j + NODE + 1 + fix_var, (double)storage[j + NODE + 1][i]);
                }
                z = evaluate(command);
                if (fabs(z) > 0.0) {
                    flag = 1;
                } else {
                    flag = 0;
                }
            }
            z = (storage[col][i] - zlo) / dz;
            index = (int32)z;
            if (index >= 0 && index < length && flag == 1) {
                my_hist[1][index] += 1.0;
            }
        }
        NCON = NCON_START;
        NSYM = NSYM_START;
        histogram_back();
        ggets_ping();
        return;
    }
    if (which == 1) {
        for (i = 0; i < storind; i++) {
            for (int32 j = 0; j < storind; j++) {
                y = storage[col][i] - storage[col][j];
                z = (y - zlo) / dz;
                index = (int32)z;
                if (index >= 0 && index < length) {
                    my_hist[1][index] += 1.0;
                }
            }
        }
        histogram_back();
        ggets_ping();
        return;
    }
    if (which == 2) {
        histogram_mycor2(storage[col], storage[col2], storind, nbins, my_hist[1], 1);
        histogram_back();
        ggets_ping();
        return;
    }
    if (which == 3) {
        histogram_fft_xcorr(storage[col], storage[col2], storind, (nbins - 1) / 2, my_hist[1], 1);
        histogram_back();
        ggets_ping();
        return;
    }
    return;
}

void
histogram_column_mean(void) {
    char bob[100];
    double sum;
    double sum2;
    double ss;
    double mean;
    double sdev;
    if (storind <= 1) {
        ggets_err_msg("Need at least 2 data points!");
        return;
    }
    if (histogram_get_col_info(&hist_inf.col, "Variable ") == 0) {
        return;
    }
    sum = 0.0;
    sum2 = 0.0;
    for (int32 i = 0; i < storind; i++) {
        ss = storage[hist_inf.col][i];
        sum += ss;
        sum2 += (ss*ss);
    }
    mean = sum / (double)storind;
    sdev = sqrt(sum2 / (double)storind - mean*mean);
    sprintf(bob, "Mean=%g Std. Dev. = %g ", mean, sdev);
    ggets_err_msg(bob);
    return;
}

int32
histogram_get_col_info(int32 *col, char *prompt) {
    char variable[20];
    if (*col == 0) {
        strcpy(variable, "t");
    } else {
        strcpy(variable, uvar_names[*col - 1]);
    }
    ggets_new_string(prompt, variable);
    browser_find_variable(variable, col);
    if (*col < 0) {
        ggets_err_msg("No such variable...");
        return 0;
    }
    return 1;
}

void
histogram_compute_power(void) {
    double s;
    double c;
    double *datx, *daty, ptot = 0;
    histogram_compute_fourier();
    if ((NEQ < 2) || (storind <= 1)) {
        return;
    }
    datx = browser_get_data_col(1);
    daty = browser_get_data_col(2);

    for (int32 i = 0; i < four_len; i++) {
        c = datx[i];
        s = daty[i];
        datx[i] = sqrt(s*s + c*c);
        daty[i] = atan2(s, c);
        ptot += (datx[i]*datx[i]);
    }
    ggets_plintf("a0=%g L2norm= %g  \n", datx[0], sqrt(ptot));
}
/* short-term fft
 * first apply a window
 * give data, window size, increment size,
 * window type - 0-square, 1=par 2=hamming,4-hanning,3- bartlet
 * returns a two vectors, real and imaginary
 * which have each of the data appended to them
 * data is data (not destroyed)
 * nr=number of points in data
 * win=window size
 * w_type=windowing
 * pow returns the power
 *     size = win/2 */

int32
histogram_spectrum(double *data, int32 nr, int32 win, int32 w_type, double *pow) {
    // assumes 50% overlap
    int32 shift = win / 2;
    int32 kwin = (nr - win + 1) / shift;
    int32 kk;
    double *ct, *st, *f, *d, x, nrmf;
    if (nr < 2) {
        return 0;
    }
    if (kwin < 1) {
        return 0;
    }
    ct = xmalloc(sizeof(*ct)*(usize)win);
    d = xmalloc(sizeof(*d)*(usize)win);
    st = xmalloc(sizeof(*st)*(usize)win);
    f = xmalloc(sizeof(*f)*(usize)win);
    nrmf = 0.0;
    for (int32 i = 0; i < win; i++) {
        x = (double)i / ((double)win);
        switch (w_type) {
        case 0:
            f[i] = 1;
            break;
        case 1:
            f[i] = x*(1 - x)*4.0;
            break;
        case 2:
            f[i] = .54 - .46*cos(2*M_PI*x);
            break;
        case 4:
            f[i] = .5*(1 - cos(2*M_PI*x));
            break;
        case 3:
            f[i] = 1 - 2*fabs(x - .5);
            break;
        default:
            break;
        }
        nrmf += (f[i]*f[i] / win);
    }
    for (int32 i = 0; i < shift; i++) {
        pow[i] = 0.0;
    }

    for (int32 j = 0; j < kwin; j++) {
        for (int32 i = 0; i < win; i++) {
            kk = (j*shift + i + nr) % nr;
            d[i] = f[i]*data[kk];
        }
        histogram_fft(d, ct, st, shift, win);
        for (int32 i = 0; i < shift; i++) {
            x = ct[i]*ct[i] + st[i]*st[i];
            pow[i] = pow[i] + sqrt(x);
        }
    }
    for (int32 i = 0; i < shift; i++) {
        pow[i] = pow[i] / ((kwin)*sqrt(nrmf));
    }
    free(f);
    free(ct);
    free(st);
    free(d);

    return 1;
}

/*  here is what we do - I think it is what MatLab does as well

    psd(x) breaks data into chunks, takes FFT of each chunk,
    power of each chunk, and averages this.

 * csd(x,y)
 * break into chunks
 * compute for each frequency histogram_fft(y)*histogram_fft(x)^*
 * now average these - note that this will be floatcomplex
 * what I call the cross spectrum is |Pxy|
 *the coherence is
 * |Pxy|^2/|Pxx||Pyy|

*/

int32
histogram_cross_spectrum(double *data, double *data2, int32 nr, int32 win, int32 w_type,
                         double *pow, int32 type) {
    int32 shift = win / 2;
    int32 kwin = (nr - win + 1) / shift;
    int32 kk;
    double *ct, *st, *f, *d, x, nrmwin;
    double *ct2, *st2, *d2;
    double *pxx;
    double *pyy;
    double *pxyr, *pxym;
    if (nr < 2) {
        return 0;
    }
    if (kwin < 1) {
        return 0;
    }
    ct = xmalloc(sizeof(*ct)*(usize)win);
    d = xmalloc(sizeof(*d)*(usize)win);
    st = xmalloc(sizeof(*st)*(usize)win);
    f = xmalloc(sizeof(*f)*(usize)win);
    ct2 = xmalloc(sizeof(*(ct2))*(usize)win);
    d2 = xmalloc(sizeof(*(d2))*(usize)win);
    st2 = xmalloc(sizeof(*(st2))*(usize)win);
    pxx = xmalloc(sizeof(*pxx)*(usize)win);
    pyy = xmalloc(sizeof(*pyy)*(usize)win);
    pxyr = xmalloc(sizeof(*pxyr)*(usize)win);
    pxym = xmalloc(sizeof(*pxym)*(usize)win);
    nrmwin = 0.0;
    for (int32 i = 0; i < win; i++) {
        x = (double)i / ((double)win);
        switch (w_type) {
        case 0:
            f[i] = 1;
            break;
        case 1:
            f[i] = x*(1 - x)*4.0;
            break;
        case 4:
            f[i] = .5*(1 - cos(2*M_PI*x));
            break;
        case 2:
            f[i] = .54 - .46*cos(2*M_PI*x);
            break;
        case 3:
            f[i] = 1 - 2*fabs(x - .5);
            break;
        default:
            break;
        }
        nrmwin += f[i]*f[i];
    }
    for (int32 i = 0; i < shift; i++) {
        pxx[i] = 0.0;
        pyy[i] = 0.0;
        pxyr[i] = 0.0;
        pxym[i] = 0.0;
    }

    for (int32 j = 0; j <= kwin; j++) {
        for (int32 i = 0; i < win; i++) {
            kk = (i + j*shift) % nr;
            d[i] = f[i]*data[kk];
            d2[i] = f[i]*data2[kk];
        }
        histogram_fft(d, ct, st, shift, win);
        histogram_fft(d2, ct2, st2, shift, win);
        for (int32 i = 0; i < shift; i++) {
            pxyr[i] += (ct[i]*ct2[i] + st[i]*st2[i]);
            pxym[i] += (ct[i]*st2[i] - ct2[i]*st[i]);
            pxx[i] += (ct[i]*ct[i] + st[i]*st[i]);
            pyy[i] += (ct2[i]*ct2[i] + st2[i]*st2[i]);
        }
    }
    for (int32 i = 0; i < shift; i++) {
        pxx[i] = pxx[i] / ((kwin)*nrmwin);
        pyy[i] = pyy[i] / ((kwin)*nrmwin);
        pxyr[i] = pxyr[i] / ((kwin)*nrmwin);
        pxym[i] = pxym[i] / ((kwin)*nrmwin);
        pxyr[i] = pxyr[i]*pxyr[i] + pxym[i]*pxym[i];
        if (type == 1) {
            pow[i] = log(pxyr[i]);
        } else {
            pow[i] = pxyr[i] / (pxx[i]*pyy[i]);
        }
    }
    free(f);
    free(ct);
    free(st);
    free(d);
    free(ct2);
    free(st2);
    free(d2);
    free(pxx);
    free(pyy);
    free(pxyr);
    free(pxym);

    return 1;
}

void
histogram_compute_sd(void) {
    int32 length;
    double total = storage[0][storind - 1] - storage[0][0];
    ggets_new_int("(0) PSDx, (1) PSDxy, (2) COHxy:", &spec_type);

    if (histogram_get_col_info(&spec_col, "Variable ") == 0) {
        return;
    }
    if (spec_type > 0) {
        if (histogram_get_col_info(&spec_col2, "Variable 2 ") == 0) {
            return;
        }
    }
    ggets_new_int("Window length ", &spec_wid);
    ggets_new_int("0:sqr 1:par 2:ham 3:bart 4:han ", &spec_win);
    if (hist_here) {
        adjoints_data_back();
        free(my_hist[0]);
        free(my_hist[1]);
        if (hist_here == 2) {
            free(my_hist[2]);
        }
        hist_here = 0;
    }
    hist_len = spec_wid / 2;
    length = hist_len + 2;
    my_hist[0] = xmalloc(sizeof(*(my_hist[0]))*(usize)length);
    my_hist[1] = xmalloc(sizeof(*(my_hist[1]))*(usize)length);
    if (my_hist[1] == NULL) {
        free(my_hist[0]);
        ggets_err_msg("Cannot allocate enough...");
        return;
    }
    hist_here = 1;
    for (int32 i = 2; i <= NEQ; i++) {
        my_hist[i] = storage[i];
    }
    for (int32 j = 0; j < hist_len; j++) {
        my_hist[0][j] = ((double)j*storind / spec_wid) / total;
    }
    if (spec_type == 0) {
        histogram_spectrum(storage[spec_col], storind, spec_wid, spec_win, my_hist[1]);
    } else {
        histogram_cross_spectrum(storage[spec_col], storage[spec_col2], storind, spec_wid, spec_win,
                                 my_hist[1], spec_type);
    }
    histogram_back();
    ggets_ping();
    return;
}

void
histogram_just_fourier(int32 flag) {
    double s;
    double c;
    double *datx, *daty;
    int32 nmodes = storind / 2 - 1;
    if (NEQ < 2 || storind <= 1) {
        return;
    }
    histogram_new_four(nmodes, spec_col);
    if (flag) {
        datx = browser_get_data_col(1);
        daty = browser_get_data_col(2);

        for (int32 i = 0; i < four_len; i++) {
            c = datx[i];
            s = daty[i];
            datx[i] = sqrt(s*s + c*c);
            daty[i] = atan2(s, c);
        }
    }
    return;
}

void
histogram_compute_fourier(void) {
    int32 nmodes = 10;
    if (NEQ < 2) {
        ggets_err_msg("Need at least three data columns");
        return;
    }
    if (storind <= 1) {
        ggets_err_msg("No data!");
        return;
    }
    if (histogram_get_col_info(&spec_col, "Variable ") == 1) {
        nmodes = storind / 2 - 1;
    }
    histogram_new_four(nmodes, spec_col);
    return;
}

void
histogram_compute_correl(void) {
    int32 lag;
    double total = storage[0][storind - 1] - storage[0][0], dta;
    dta = total / (double)(storind - 1);

    ggets_new_int("Number of bins ", &hist_inf.nbins);
    ggets_new_int("(0)Direct or (1) FFT ", &hist_inf.fftc);
    if (hist_inf.nbins > (storind / 2 - 1)) {
        hist_inf.nbins = storind / 2 - 2;
    }

    hist_inf.nbins = 2*(hist_inf.nbins / 2) + 1;
    lag = hist_inf.nbins / 2;

    // lets try to get the lags correct for plotting
    hist_inf.xlo = -lag*dta;
    hist_inf.xhi = lag*dta;

    if (histogram_get_col_info(&hist_inf.col, "Variable 1 ") == 0) {
        return;
    }
    if (histogram_get_col_info(&hist_inf.col2, "Variable 2 ") == 0) {
        return;
    }
    histogram_new(hist_inf.nbins, hist_inf.xlo, hist_inf.xhi, hist_inf.col, hist_inf.col2,
                  hist_inf.cond, 2 + hist_inf.fftc);
    return;
}

void
histogram_compute_stacor(void) {
    ggets_new_int("Number of bins ", &hist_inf.nbins);
    ggets_new_float("Low ", &hist_inf.xlo);
    ggets_new_float("Hi ", &hist_inf.xhi);
    if (histogram_get_col_info(&hist_inf.col, "Variable ") == 0) {
        return;
    }
    histogram_new(hist_inf.nbins, hist_inf.xlo, hist_inf.xhi, hist_inf.col, 0, hist_inf.cond, 1);
    return;
}

void
histogram_mycor(double *x, double *y, int32 n, double zlo, double zhi, int32 nbins, double *z,
                int32 flag) {
    int32 k;
    int32 count = 0;
    double sum;
    double avx = 0.0;
    double avy = 0.0;
    double dz = (zhi - zlo) / (double)nbins, jz;
    if (flag) {
        for (int32 i = 0; i < n; i++) {
            avx += x[i];
            avy += y[i];
        }
        avx = avx / (double)n;
        avy = avy / (double)n;
    }
    for (int32 j = 0; j <= nbins; j++) {
        sum = 0.0;
        count = 0;
        jz = dz*j + zlo;
        for (int32 i = 0; i < n; i++) {
            k = i + (int32)jz;
            if ((k >= 0) && (k < n)) {
                count++;
                sum += (x[i] - avx)*(y[k] - avy);
            }
        }
        if (count > 0) {
            sum = sum / count;
        }
        z[j] = sum;
    }
    return;
}

void
histogram_mycor2(double *x, double *y, int32 n, int32 nbins, double *z, int32 flag) {
    int32 k, count = 0, lag = nbins / 2;
    double sum;
    double avx = 0.0;
    double avy = 0.0;
    if (flag) {
        for (int32 i = 0; i < n; i++) {
            avx += x[i];
            avy += y[i];
        }
        avx = avx / (double)n;
        avy = avy / (double)n;
    }
    for (int32 j = 0; j <= nbins; j++) {
        sum = 0.0;
        count = 0;
        for (int32 i = 0; i < n; i++) {
            k = i + j - lag;
            k = (k + n) % n;
            if ((k >= 0) && (k < n)) {
                count++;
                sum += (x[i] - avx)*(y[k] - avy);
            }
        }
        if (count > 0) {
            sum = sum / count;
        }
        z[j] = sum;
    }
    return;
}

void
histogram_compute(void) {
    ggets_new_int("Number of bins ", &hist_inf.nbins);
    ggets_new_float("Low ", &hist_inf.xlo);
    ggets_new_float("Hi ", &hist_inf.xhi);
    if (histogram_get_col_info(&hist_inf.col, "Variable ") == 0) {
        return;
    }
    ggets_new_string("Condition ", hist_inf.cond);
    histogram_new(hist_inf.nbins, hist_inf.xlo, hist_inf.xhi, hist_inf.col, 0, hist_inf.cond, 0);
    return;
}

/* experimental -- does it work */
/* nlag should be less than length/2 */
void
histogram_fft_xcorr(double *data1, double *data2, int32 length, int32 nlag, double *cr,
                    int32 flag) {
    double *re1, *re2, *im1, *im2, x, y, sum;
    double av1 = 0.0;
    double av2 = 0.0;
    int32 dim[2];
    if (flag) {
        for (int32 i = 0; i < length; i++) {
            av1 += data1[i];
            av2 += data2[i];
        }
        av1 = av1 / (double)length;
        av2 = av2 / (double)length;
    }

    dim[0] = length;
    re1 = xmalloc((usize)length*sizeof(*(re1)));
    im1 = xmalloc((usize)length*sizeof(*(im1)));
    re2 = xmalloc((usize)length*sizeof(*(re2)));
    im2 = xmalloc((usize)length*sizeof(*(im2)));

    for (int32 i = 0; i < length; i++) {
        im1[i] = 0.0;
        re1[i] = (data1[i] - av1);
        im2[i] = 0.0;
        re2[i] = (data2[i] - av2);
    }

    fftn(1, dim, re1, im1, 1, -1);
    fftn(1, dim, re2, im2, 1, -1);
    for (int32 i = 0; i < length; i++) {
        x = re1[i]*re2[i] + im1[i]*im2[i];
        y = im1[i]*re2[i] - im2[i]*re1[i];
        re1[i] = x;
        im1[i] = -y;
    }
    fftn(1, dim, re1, im1, -1, -1);
    /* now lets order these
       I think!  */
    sum = 0.0;
    for (int32 i = 0; i < nlag; i++) {
        sum += fabs(im1[i]);
        cr[nlag + i] = (double)re1[i]*length;  // positive part of the correlation
    }
    for (int32 i = 0; i < nlag; i++) {
        sum += fabs(im1[length - nlag + i]);
        cr[i] = (double)re1[length - nlag + i]*length;
    }
    free(re1);
    free(re2);
    free(im1);
    free(im2);
    ggets_plintf("residual = %g\n", sum);
    return;
}

void
histogram_fft(double *data, double *ct, double *st, int32 nmodes, int32 length) {
    double *im;
    double *re;
    int32 dim[2];
    dim[0] = length;
    re = xmalloc((usize)length*sizeof(*re));
    im = xmalloc((usize)length*sizeof(*im));
    for (int32 i = 0; i < length; i++) {
        im[i] = 0.0;
        re[i] = data[i];
    }

    fftn(1, dim, re, im, 1, -1);
    ct[0] = re[0];
    st[0] = 0.0;
    for (int32 i = 1; i < nmodes; i++) {
        ct[i] = re[i]*2.0;
        st[i] = im[i]*2.0;
    }
    free(im);
    free(re);
}
