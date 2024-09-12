#ifndef _load_eqn_h_
#define _load_eqn_h_

#include <stdio.h>
#include "integers.h"

/*
The acutual max filename length is determined by the
FILENAME_MAX (see <stdio.h>), and usually 4096 -- but
this is huge and usually overkill.  On the otherhand
the old Xpp default string buffer size of 100 is a bit
restricitive for lengths of filenames. You could also
set this define in the Makefile or at compile time to
override the below definition.
*/

#ifndef XPP_MAX_NAME
#define XPP_MAX_NAME 300
#if (XPP_MAX_NAME > FILENAME_MAX)
#define XPP_MAX_NAME FILENAME_MAX
#endif
#endif

/*
Options are set accroding to an order of precedence

command line < mfile < .xpprc < default.opt

Add any options here that you might want to track.
*/
typedef struct {
    int32 BIG_FONT_NAME;
    int32 SMALL_FONT_NAME;
    int32 BACKGROUND;
    int32 IXPLT;
    int32 IYPLT;
    int32 IZPLT;
    int32 AXES;
    int32 NMESH;
    int32 METHOD;
    int32 TIMEPLOT;
    int32 MAXSTOR;
    int32 TEND;
    int32 DT;
    int32 T0;
    int32 TRANS;
    int32 BOUND;
    int32 TOLER;
    int32 DELAY;
    int32 XLO;
    int32 XHI;
    int32 YLO;
    int32 YHI;
    int32 UserBlack;
    int32 UserWhite;
    int32 UserMainWinColor;
    int32 UserDrawWinColor;
    int32 UserGradients;
    int32 UserBGBitmap;
    int32 UserMinWidth;
    int32 UserMinHeight;
    int32 YNullColor;
    int32 XNullColor;
    int32 StableManifoldColor;
    int32 UnstableManifoldColor;
    int32 START_LINE_TYPE;
    int32 RandSeed;
    int32 PaperWhite;
    int32 COLORMAP;
    int32 NPLOT;
    int32 DLL_LIB;
    int32 DLL_FUN;
    int32 XP;
    int32 YP;
    int32 ZP;
    int32 NOUT;
    int32 VMAXPTS;
    int32 TOR_PER;
    int32 JAC_EPS;
    int32 NEWT_TOL;
    int32 NEWT_ITER;
    int32 FOLD;
    int32 DTMIN;
    int32 DTMAX;
    int32 ATOL;
    int32 TOL;
    int32 BANDUP;
    int32 BANDLO;
    int32 PHI;
    int32 THETA;
    int32 XMIN;
    int32 XMAX;
    int32 YMIN;
    int32 YMAX;
    int32 ZMIN;
    int32 ZMAX;
    int32 POIVAR;
    int32 OUTPUT;
    int32 POISGN;
    int32 POIEXT;
    int32 POISTOP;
    int32 STOCH;
    int32 POIPLN;
    int32 POIMAP;
    int32 RANGEOVER;
    int32 RANGESTEP;
    int32 RANGELOW;
    int32 RANGEHIGH;
    int32 RANGERESET;
    int32 RANGEOLDIC;
    int32 RANGE;
    int32 NTST;
    int32 NMAX;
    int32 NPR;
    int32 NCOL;
    int32 DSMIN;
    int32 DSMAX;
    int32 DS;
    int32 PARMAX;
    int32 NORMMIN;
    int32 NORMMAX;
    int32 EPSL;
    int32 EPSU;
    int32 EPSS;
    int32 RUNNOW;
    int32 SEC;
    int32 UEC;
    int32 SPC;
    int32 UPC;
    int32 AUTOEVAL;
    int32 AUTOXMAX;
    int32 AUTOYMAX;
    int32 AUTOXMIN;
    int32 AUTOYMIN;
    int32 AUTOVAR;
    int32 PS_FONT;
    int32 PS_LW;
    int32 PS_FSIZE;
    int32 PS_COLOR;
    int32 FOREVER;
    int32 BVP_TOL;
    int32 BVP_EPS;
    int32 BVP_MAXIT;
    int32 BVP_FLAG;
    int32 SOS;
    int32 FFT;
    int32 HIST;
    int32 PltFmtFlag;
    int32 ATOLER;
    int32 MaxEulIter;
    int32 EulTol;
    int32 EVEC_ITER;
    int32 EVEC_ERR;
    int32 NULL_ERR;
    int32 NEWT_ERR;
    int32 NULL_HERE;
    int32 TUTORIAL;
    int32 SLIDER1;
    int32 SLIDER2;
    int32 SLIDER3;
    int32 SLIDER1LO;
    int32 SLIDER2LO;
    int32 SLIDER3LO;
    int32 SLIDER1HI;
    int32 SLIDER2HI;
    int32 SLIDER3HI;
    int32 POSTPROCESS;
    int32 HISTCOL;
    int32 HISTLO;
    int32 HISTHI;
    int32 HISTBINS;
    int32 SPECCOL;
    int32 SPECCOL2;
    int32 SPECWIDTH;
    int32 SPECWIN;
    int32 PLOTFORMAT;
    int32 DFGRID;
    int32 DFBATCH;
    int32 NCBATCH;
    int32 COLORVIA;
    int32 COLORIZE;
    int32 COLORLO;
    int32 COLORHI;
    int32 HISTCOL2;
    int32 HISTLO2;
    int32 HISTHI2;
    int32 HISTBINS2;

} OptionsSet;

void dump_torus(FILE *fp, int32 f);
void load_eqn(void);
void set_X_vals(void);
void set_all_vals(void);
void read_defaults(FILE *fp);
void fil_flt(FILE *fpt, double *val);
void fil_int(FILE *fpt, int32 *val);
void add_intern_set(char *name, char *does);
void extract_action(char *ptr);
void extract_internset(int32 j);
void do_intern_set(char *name1, char *value);
int32 msc(char *s1, char *s2);
void set_internopts(OptionsSet *mask);
void set_internopts_xpprc_and_comline(void);
void split_apart(char *bob, char *name, char *value);
void check_for_xpprc(void);
void stor_internopts(char *s1);
void set_option(char *s1, char *s2, int32 force, OptionsSet *mask);

#endif
