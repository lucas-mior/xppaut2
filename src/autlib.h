typedef struct GlobalScratch {
    double *dfu, *dfp, *uu1, *uu2, *ff1, *ff2;
} GlobalScratch;

typedef struct GlobalRotations {
    int64 irtn;
    int64 *nrtn;
} GlobalRotations;

typedef struct GlobalParameters {
    rap_type *rav;
    iap_type *iav;
    double *dtv;
} GlobalParameters;

extern GlobalScratch global_scratch;
extern GlobalRotations global_rotations;
extern GlobalParameters global_parameters;

extern int32 fp8_is_open;
