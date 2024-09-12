#ifndef _histogram_h_
#define _histogram_h_
#include "integers.h"

int32 two_d_hist(int32 col1, int32 col2, int32 ndat, int32 n1, int32 n2,
                 double xlo, double xhi, double ylo, double yhi);
void four_back(void);
void hist_back(void);
void new_four(int32 nmodes, int32 col);
int32 new_2d_hist(void);
void new_hist(int32 nbins, double zlo, double zhi, int32 col, int32 col2,
              char *condition, int32 which);
void column_mean(void);
int32 get_col_info(int32 *col, char *prompt);
void compute_power(void);
int32 spectrum(float *data, int32 nr, int32 win, int32 w_type, float *pow);
int32 cross_spectrum(float *data, float *data2, int32 nr, int32 win,
                     int32 w_type, float *pow, int32 type);
void compute_sd(void);
void compute_fourier(void);
void compute_correl(void);
void compute_stacor(void);
void mycor(float *x, float *y, int32 n, double zlo, double zhi, int32 nbins,
           float *z, int32 flag);
void mycor2(float *x, float *y, int32 n, int32 nbins, float *z, int32 flag);
void compute_hist(void);
void sft(float *data, float *ct, float *st, int32 nmodes, int32 grid);
void fftxcorr(float *data1, float *data2, int32 length, int32 nlag, float *cr,
              int32 flag);
void fft(float *data, float *ct, float *st, int32 nmodes, int32 length);
void post_process_stuff();
void just_fourier(int32 flag);
void just_sd(int32 flag);

#endif
