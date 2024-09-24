#ifndef _cv2_h
#define _cv2_h
#include "integers.h"

void start_cv(double *y, double t, int32 n, double *atol, double *rtol);
void cv_end(void);
void cvode_err_msg(int32 kflag);
int32 cvode(int32 *command, double *y, double *t, int32 n, double tout,
            int32 *kflag, double *atol, double *rtol);
int32 ccvode(int32 *command, double *y, double *t, int32 n, double tout,
             int32 *kflag, double *atol, double *rtol);

#endif
