#include "integers.h"

typedef struct {
    int32 ndim;
    int32 ips;
    int32 irs;
    int32 ilp;
    int32 nicp;
    int32 icp[100];
    int32 ntst;
    int32 ncol;
    int32 iad;
    int32 isp;
    int32 isw;
    int32 iplt;
    int32 nbc;
    int32 nint;
    int32 nmx;
    double rl0;
    double rl1;
    double a0;
    double a1;
    int32 npr;
    int32 mxbf;
    int32 iid;
    int32 itmx;
    int32 itnw;
    int32 nwtn;
    int32 jac;
    double epsl;
    double epss;
    double epsu;
    double ds;
    double dsmax;
    double dsmin;
    int32 iads;
    int32 nthl;
    int32 ithl[100];
    double thl[100];
    int32 nuzr;
    int32 iuz[20];
    double vuz[20];
    //  these are for the homoclinic stuff
    int32 nunstab;
    int32 nstab;
    int32 iequib;  // +1 homoclinic   -2 heteroclinic

} XAuto;

extern XAuto x_auto;

/*  for homcont  - itwist=0, istart=2, nrev=0,nfixed=0,npsi=0 */
