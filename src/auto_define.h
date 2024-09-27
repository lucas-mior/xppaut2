#include "integers.h"
struct {
    int32 ndim;
    int32 ips;
    int32 irs;
    int32 ilp;
    int32 icp[20];
    double par[20];
} blbcn_;

#define blbcn_1 blbcn_
struct {
    int32 ntst;
    int32 ncol;
    int32 iad;
    int32 isp;
    int32 isw;
    int32 iplt;
    int32 nbc;
    int32 nint;
} blcde_;

#define blcde_1 blcde_

struct {
    double ds;
    double dsmin;
    double dsmax;
    int32 iads;
} bldls_;

#define bldls_1 bldls_

struct {
    int32 nmx;
    int32 nuzr;
    double rl0;
    double rl1;
    double a0;
    double a1;
} bllim_;

#define bllim_1 bllim_

struct {
    int32 npr;
    int32 mxbf;
    int32 iid;
    int32 itmx;
    int32 itnw;
    int32 nwtn;
    int32 jac;
} blmax_;

#define blmax_1 blmax_

struct {
    double epsl[20];
    double epsu;
    double epss;
} bleps_;

#define bleps_1 bleps_

#define EPSU bleps_1.epsu
#define EPSS bleps_1.epss
#define EPSL(a) bleps_1.epsl[(a)]

#define IRS blbcn_1.irs
#define NDIM blbcn_1.ndim
#define IPS blbcn_1.ips
#define ILP blbcn_1.ilp

#define NTST blcde_1.ntst
#define NCOL blcde_1.ncol
#define IAD blcde_1.iad
#define ISP blcde_1.isp
#define ISW blcde_1.isw
#define NBC blcde_1.nbc
#define NIC blcde_1.nint

#define DS bldls_1.ds
#define DSMAX bldls_1.dsmax
#define DSMIN bldls_1.dsmin

#define NMX bllim_1.nmx
#define NUZR bllim_1.nuzr
#define RL0 bllim_1.rl0
#define RL1 bllim_1.rl1
#define AUTO_A0 bllim_1.a0
#define AUTO_A1 bllim_1.a1
