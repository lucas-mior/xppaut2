/*
 *  Copyright (C) 2002-2017  Bard Ermentrout & Daniel Dougherty

 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.

 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

 *  The author can be contacted at
 *   bard@pitt.edu
 */

#include <X11/Xatom.h>
#include <X11/Xlib.h>
#include <X11/Xos.h>
#include <X11/Xproto.h>
#include <X11/Xutil.h>
#include <dirent.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "auto_nox.h"
#include "functions.h"
#include "integers.h"
#include "pp.bitmap"
#include "read_dir.h"

#define FIX_SIZE 3
#define FIX_MIN_SIZE 2
#define FIX_MAX_SIZE 1
#define LOWBIT(x) ((x) & (~(x) + 1))

#include <X11/cursorfont.h>

#include "myfonts.h"

int32 allwinvis = 0;
int32 use_ani_file = 0;
char anifile[XPP_MAX_NAME];

double xppvermaj;
double xppvermin;

int32 Xup;
int32 TipsFlag = 1;
Atom deleteWindowAtom = 0;
int32 XPPBatch = 0;
int32 batch_range = 0;
int32 BatchEquil = -1;
char batchout[XPP_MAX_NAME];
char UserOUTFILE[XPP_MAX_NAME];
int32 DisplayHeight;
int32 DisplayWidth;
int32 TrueColorFlag;
char big_font_name[100];
char small_font_name[100];
char PlotFormat[10];

int32 PaperWhite = -1;

Window draw_win;
Window main_win;
Window command_pop;
Window info_pop;
GC gc;
GC gc_graph;
GC small_gc;
GC font_gc;
char UserBlack[8];
char UserWhite[8];
char UserMainWinColor[8];
char UserDrawWinColor[8];
char UserBGBitmap[XPP_MAX_NAME];

int32 UserGradients = -1;
int32 UserMinWidth = 0;
int32 UserMinHeight = 0;
uint32 MyBackColor;
uint32 MyForeColor;
uint32 MyMainWinColor;
uint32 MyDrawWinColor;
uint32 GrFore;
uint32 GrBack;
int32 SCALEY;
Display *display;
int32 screen;
int32 periodic = 0;
int32 DCURYb;
int32 CURY_OFFb;
int32 DCURYs;
int32 DCURXs;
int32 CURY_OFFs;
int32 DCURY;
int32 DCURX;
int32 CURY_OFF;
FILE *logfile;
int32 XPPVERBOSE = 1;
int32 OVERRIDE_QUIET = 0;
int32 OVERRIDE_LOGFILE = 0;
int32 tfBell;

int32 SLIDER1 = -1;
int32 SLIDER2 = -1;
int32 SLIDER3 = -1;
char SLIDER1VAR[20];
char SLIDER2VAR[20];
char SLIDER3VAR[20];
double SLIDER1LO = 0.0;
double SLIDER2LO = 0.0;
double SLIDER3LO = 0.0;
double SLIDER1HI = 1.0;
double SLIDER2HI = 1.0;
double SLIDER3HI = 1.0;

/* Set this to 1 if you want the tutorial to come up at start-up
 * as default behavior */
int32 DoTutorial = 0;

OptionsSet notAlreadySet;
XFontStruct *small_font;

static int32 DCURXb;
static int32 SCALEX;
static Window TopButton[6];
static XFontStruct *big_font;

static int32 main_get_x_colors(XWindowAttributes *win_info, XColor **colors);
static void main_get_gc(GC *gc);
static void main_top_button_events(XEvent report);
static void main_top_button_cross(Window window, int32 b);
static void main_xpp_events(XEvent report, int32 min_wid, int32 min_hgt);
static void main_init_x(void);
static void main_check_for_quiet(int32 argc, char **argv);

void
main_plot_command(int32 nit, int32 icount, int32 cwidth) {
    double dx;

    if (nit == 0)
        return;

    dx = (double)icount*(double)cwidth / (double)nit;
    XDrawPoint(display, command_pop, gc, (int32)dx, 5);
    return;
}

void *
XMALLOC(usize size, const char *function, int32 line) {
    void *p;

    if (!(p = malloc(size))) {
        fprintf(stderr, "Error allocating %zu bytes in function %s, line %d.\n",
                size, function, line);
        exit(EXIT_FAILURE);
    }

    return p;
}

#ifndef MALLOC_DEBUG
void *
xmalloc(usize size) {
    void *p;

    if (!(p = malloc(size))) {
        fprintf(stderr, "Error allocating %zu bytes.\n", size);
        exit(EXIT_FAILURE);
    }

    return p;
}
#endif

int32
main_my_abort(void) {
    int32 ch;

    while (XPending(display) > 0) {
        XEvent event;
        XNextEvent(display, &event);

        if (ani_check_pause(event) == 27)
            return 27;

        switch (event.type) {
        case Expose:
            many_pops_do_expose(event);
            break;
        case ButtonPress:
            break;
        case KeyPress:
            ch = ggets_get_key_press(&event);
            return ch;
        default:
            break;
        }

        return 0;
    }

    return 64;
}

void
do_main(int32 argc, char **argv) {
    char myfile[XPP_MAX_NAME];
    char pptitle[sizeof(this_file) + 40];
    uint32 min_wid = 450, min_hgt = 360;
    OptionsSet *tempNS;
    /* Track which options have not been set already */
    notAlreadySet.BIG_FONT_NAME = 1;
    notAlreadySet.SMALL_FONT_NAME = 1;
    notAlreadySet.BACKGROUND = 1;
    notAlreadySet.IXPLT = 1;
    notAlreadySet.IYPLT = 1;
    notAlreadySet.IZPLT = 1;
    notAlreadySet.AXES = 1;
    notAlreadySet.NMESH = 1;
    notAlreadySet.METHOD = 1;
    notAlreadySet.TIMEPLOT = 1;
    notAlreadySet.MAXSTOR = 1;
    notAlreadySet.TEND = 1;
    notAlreadySet.DT = 1;
    notAlreadySet.T0 = 1;
    notAlreadySet.TRANS = 1;
    notAlreadySet.BOUND = 1;
    notAlreadySet.TOLER = 1;
    notAlreadySet.DELAY = 1;
    notAlreadySet.XLO = 1;
    notAlreadySet.XHI = 1;
    notAlreadySet.YLO = 1;
    notAlreadySet.YHI = 1;
    notAlreadySet.UserBlack = 1;
    notAlreadySet.UserWhite = 1;
    notAlreadySet.UserMainWinColor = 1;
    notAlreadySet.UserDrawWinColor = 1;
    notAlreadySet.UserGradients = 1;
    notAlreadySet.UserBGBitmap = 1;
    notAlreadySet.UserMinWidth = 1;
    notAlreadySet.UserMinHeight = 1;
    notAlreadySet.YNullColor = 1;
    notAlreadySet.XNullColor = 1;
    notAlreadySet.StableManifoldColor = 1;
    notAlreadySet.UnstableManifoldColor = 1;
    notAlreadySet.START_LINE_TYPE = 1;
    notAlreadySet.RandSeed = 1;
    notAlreadySet.PaperWhite = 1;
    notAlreadySet.COLORMAP = 1;
    notAlreadySet.NPLOT = 1;
    notAlreadySet.DLL_LIB = 1;
    notAlreadySet.DLL_FUN = 1;
    notAlreadySet.XP = 1;
    notAlreadySet.YP = 1;
    notAlreadySet.ZP = 1;
    notAlreadySet.NOUT = 1;
    notAlreadySet.VMAXPTS = 1;
    notAlreadySet.TOR_PER = 1;
    notAlreadySet.JAC_EPS = 1;
    notAlreadySet.NEWT_TOL = 1;
    notAlreadySet.NEWT_ITER = 1;
    notAlreadySet.FOLD = 1;
    notAlreadySet.DTMIN = 1;
    notAlreadySet.DTMAX = 1;
    notAlreadySet.ATOL = 1;
    notAlreadySet.TOL = 1;
    notAlreadySet.BANDUP = 1;
    notAlreadySet.BANDLO = 1;
    notAlreadySet.PHI = 1;
    notAlreadySet.THETA = 1;
    notAlreadySet.XMIN = 1;
    notAlreadySet.XMAX = 1;
    notAlreadySet.YMIN = 1;
    notAlreadySet.YMAX = 1;
    notAlreadySet.ZMIN = 1;
    notAlreadySet.ZMAX = 1;
    notAlreadySet.POIVAR = 1;
    notAlreadySet.OUTPUT = 1;
    notAlreadySet.POISGN = 1;
    notAlreadySet.POISTOP = 1;
    notAlreadySet.STOCH = 1;
    notAlreadySet.POIPLN = 1;
    notAlreadySet.POIMAP = 1;
    notAlreadySet.RANGEOVER = 1;
    notAlreadySet.RANGESTEP = 1;
    notAlreadySet.RANGELOW = 1;
    notAlreadySet.RANGEHIGH = 1;
    notAlreadySet.RANGERESET = 1;
    notAlreadySet.RANGEOLDIC = 1;
    notAlreadySet.RANGE = 1;
    notAlreadySet.NTST = 1;
    notAlreadySet.NMAX = 1;
    notAlreadySet.NPR = 1;
    notAlreadySet.NCOL = 1;
    notAlreadySet.DSMIN = 1;
    notAlreadySet.DSMAX = 1;
    notAlreadySet.DS = 1;
    notAlreadySet.PARMAX = 1;
    notAlreadySet.NORMMIN = 1;
    notAlreadySet.NORMMAX = 1;
    notAlreadySet.EPSL = 1;
    notAlreadySet.EPSU = 1;
    notAlreadySet.EPSS = 1;
    notAlreadySet.RUNNOW = 1;
    notAlreadySet.SEC = 1;
    notAlreadySet.UEC = 1;
    notAlreadySet.SPC = 1;
    notAlreadySet.UPC = 1;
    notAlreadySet.AUTOEVAL = 1;
    notAlreadySet.AUTOXMAX = 1;
    notAlreadySet.AUTOYMAX = 1;
    notAlreadySet.AUTOXMIN = 1;
    notAlreadySet.AUTOYMIN = 1;
    notAlreadySet.AUTOVAR = 1;
    notAlreadySet.PS_FONT = 1;
    notAlreadySet.PS_LW = 1;
    notAlreadySet.PS_FSIZE = 1;
    notAlreadySet.PS_COLOR = 1;
    notAlreadySet.FOREVER = 1;
    notAlreadySet.BVP_TOL = 1;
    notAlreadySet.BVP_EPS = 1;
    notAlreadySet.BVP_MAXIT = 1;
    notAlreadySet.BVP_FLAG = 1;
    notAlreadySet.SOS = 1;
    notAlreadySet.FFT = 1;
    notAlreadySet.HIST = 1;
    notAlreadySet.PltFmtFlag = 1;
    notAlreadySet.ATOLER = 1;
    notAlreadySet.euler_max_iter = 1;
    notAlreadySet.euler_tol = 1;
    notAlreadySet.EVEC_ITER = 1;
    notAlreadySet.EVEC_ERR = 1;
    notAlreadySet.NULL_ERR = 1;
    notAlreadySet.NEWT_ERR = 1;
    notAlreadySet.NULL_HERE = 1;
    notAlreadySet.TUTORIAL = 1;
    notAlreadySet.SLIDER1 = 1;
    notAlreadySet.SLIDER2 = 1;
    notAlreadySet.SLIDER3 = 1;
    notAlreadySet.SLIDER1LO = 1;
    notAlreadySet.SLIDER2LO = 1;
    notAlreadySet.SLIDER3LO = 1;
    notAlreadySet.SLIDER1HI = 1;
    notAlreadySet.SLIDER2HI = 1;
    notAlreadySet.SLIDER3HI = 1;
    notAlreadySet.POSTPROCESS = 1;
    notAlreadySet.HISTCOL = 1;
    notAlreadySet.HISTLO = 1;
    notAlreadySet.HISTHI = 1;
    notAlreadySet.HISTBINS = 1;
    notAlreadySet.HISTCOL2 = 1;
    notAlreadySet.HISTLO2 = 1;
    notAlreadySet.HISTHI2 = 1;
    notAlreadySet.HISTBINS2 = 1;
    notAlreadySet.SPECCOL = 1;
    notAlreadySet.SPECCOL2 = 1;
    notAlreadySet.SPECWIDTH = 1;
    notAlreadySet.SPECWIN = 1;
    notAlreadySet.PLOTFORMAT = 1;
    notAlreadySet.DFGRID = 1;
    notAlreadySet.DFBATCH = 1;
    notAlreadySet.NCBATCH = 1;
    notAlreadySet.COLORVIA = 1;
    notAlreadySet.COLORIZE = 1;
    notAlreadySet.COLORHI = 1;
    notAlreadySet.COLORLO = 1;

    get_directory(myfile);

    SCALEX = 640;
    SCALEY = 480;

    Xup = 0;
    snprintf(batchout, sizeof(batchout), "output.dat");
    snprintf(PlotFormat, sizeof(PlotFormat), "ps");

    /* Read visualization environement variables here
       since some may be overridden by command line */
    logfile = stdout;
    main_check_for_quiet(argc, argv);

    comline_do(argc, argv);

    /* We need to init_X here if there is no file on command line
     * so that a file browser can be opened.  */
    if (!XPPBatch) {
        /* Swap out the current options for a temporary place holder */
        tempNS = xmalloc(sizeof(OptionsSet));
        *tempNS = notAlreadySet;
        /* Initialize what's needed to open a browser based on
         * the current options.  */
        main_do_vis_env();
        load_eqn_set_all_vals();
        main_init_x();
        /*       XSynchronize(display,1); */
        /* Now swap back the options for proper precedence ordering of options.
         */
        notAlreadySet = *tempNS;
        free(tempNS);
    }

    load_eqn();

    tempNS = xmalloc(sizeof(*tempNS));
    *tempNS = notAlreadySet;
    load_eqn_set_internopts(tempNS);
    free(tempNS);

    storage_init_alloc_info();
    main_do_vis_env();
    load_eqn_set_all_vals();

    storage_init_alloc_info();
    dae_fun_set_init_guess();
    simplenet_update_all_ffts();

#ifdef AUTO
    auto_nox_init_win();
#endif

    if (form_ode_idsc(this_file))
        METHOD = 0;
    xppvermaj = (double)MAJOR_VERSION;
    xppvermin = (double)MINOR_VERSION;
    if (strlen(this_file) < 60) {
        snprintf(pptitle, sizeof(pptitle), "XPP Ver %g.%g >> %s", xppvermaj,
                 xppvermin, this_file);
    } else {
        snprintf(pptitle, sizeof(pptitle), "XPP Version %g.%g", xppvermaj,
                 xppvermin);
    }
    numerics_do_meth();

    numerics_set_delay();
    rhs_function = my_rhs;
    do_fit_init_info();
    form_ode_strip_saveqn();
    form_ode_create_plot_list();
    extra_auto_load_dll();

    if (XPPBatch) {
        color_map_make();
        init_browser();
        graphics_init_all();
        comline_if_needed_load_set();
        comline_if_needed_load_par();
        comline_if_needed_load_ic();
        comline_if_needed_select_sets();
        comline_if_needed_load_ext_options();
        graphics_set_extra();
        nullcline_set_colorization_stuff();
        integrate_batch();
        if (NCBatch > 0)
            silent_nullclines();
        if (DFBatch > 0)
            nullcline_silent_dfields();
        integrate_silent_equilibria();
        exit(0);
    }

    many_pops_gtitle_text(pptitle, main_win);
    Xup = 1;
    color_map_make();

    XMapWindow(display, main_win);

    {
        /* main make pops */
        int32 x;
        int32 y;
        uint32 h, w, bw, d;
        Window wn;
        XGetGeometry(display, main_win, &wn, &x, &y, &w, &h, &bw, &d);
        create_the_menus(main_win);
        command_pop =
            XCreateSimpleWindow(display, main_win, 0, DCURYs + 4, w - 2,
                                (uint)DCURY + 4, 2, MyForeColor, MyBackColor);
        info_pop = XCreateSimpleWindow(display, main_win, 0,
                                       (int32)h - DCURY - 4, w - 2, (uint)DCURY,
                                       2, MyForeColor, MyBackColor);
        XCreateFontCursor(display, XC_hand2);
        XSelectInput(display, command_pop,
                     KeyPressMask | ButtonPressMask | ExposureMask);
        XSelectInput(display, info_pop, ExposureMask);
        XMapWindow(display, info_pop);
        XMapWindow(display, command_pop);
        many_pops_init_grafs(16*DCURX + 6, DCURYs + DCURYb + 6,
                             (int32)w - 16 - 16*DCURX,
                             (int32)h - 6*DCURY - 16);
        init_conds_create_par_sliders(main_win, 0, (int32)h - 5*DCURY + 8);
        graphics_get_draw_area();
    }

    {
        /* main make top buttons */
        int32 x1 = 2, x2 = 6*DCURXs + 5, dx = DCURXs;
        TopButton[0] = make_fancy_window(main_win, x1, 1, x2, DCURYs, 1);
        x1 += x2 + dx;
        TopButton[1] = make_fancy_window(main_win, x1, 1, x2, DCURYs, 1);
        x1 += x2 + dx;

        TopButton[2] = make_fancy_window(main_win, x1, 1, x2, DCURYs, 1);
        x1 += x2 + dx;

        TopButton[3] = make_fancy_window(main_win, x1, 1, x2, DCURYs, 1);
        x1 += x2 + dx;

        TopButton[4] = make_fancy_window(main_win, x1, 1, x2, DCURYs, 1);
        x1 += x2 + dx;

        TopButton[5] = make_fancy_window(main_win, x1, 1, x2, DCURYs, 1);
        x1 += x2 + dx;
    }

    init_conds_initialize_box();

    init_browser();
    if (allwinvis == 1) {
        init_conds_make_new_ic_box();
        init_conds_make_new_bc_box();
        init_conds_make_new_delay_box();
        init_conds_make_new_param_box();
        make_new_browser();
        eig_list_create_eq_list();
    }

    Xup = 1;
    ani_zero();
    graphics_set_extra();
    nullcline_set_colorization_stuff();

    pop_list_make_scrbox_lists();

    /*          MAIN LOOP             */
    {
        /* main test color info */
        XColor *colors;
        XWindowAttributes xwa;
        TrueColorFlag = 0;

        XGetWindowAttributes(display, main_win, &xwa);
        main_get_x_colors(&xwa, &colors);

        if (colors)
            free((char *)colors);
    }
    comline_if_needed_load_set();
    comline_if_needed_load_par();
    comline_if_needed_load_ic();
    comline_if_needed_load_ext_options();
    if (use_ani_file) {
        ani_new_vcr();
        ani_get_file(anifile);
    }

    if (DoTutorial == 1)
        menudrive_do_tutorial();

    graf_par_default_window();

    main_do_events(min_wid, min_hgt);
}

void
main_check_for_quiet(int32 argc, char **argv) {
    /* First scan, check for any QUIET option set... */
    int32 i = 0;
    /* Allow for multiple calls to the QUIET and LOGFILE options
     * on the command line. The last setting is the one that will stick.
     * Settings of logfile and quiet in the xpprc file will be ignored
     * if they are set on the command line.  */
    int32 quiet_specified_once = 0;
    int32 logfile_specified_once = 0;

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-quiet") == 0) {
            load_eqn_set_option("QUIET", argv[i + 1], 1, NULL);
            quiet_specified_once = 1;
            i++;
        } else if (strcmp(argv[i], "-logfile") == 0) {
            load_eqn_set_option("LOGFILE", argv[i + 1], 1, NULL);
            logfile_specified_once = 1;
            i++;
        }
    }
    /* If -quiet or -logfile were specified at least once on the command line
     * we lock those in now...  */
    if (quiet_specified_once == 1)
        OVERRIDE_QUIET = 1;
    if (logfile_specified_once == 1)
        OVERRIDE_LOGFILE = 1;
    return;
}

void
main_do_vis_env(void) {
    load_eqn_set_x_vals();
    load_eqn_check_for_xpprc();
    load_eqn_set_internopts_xpprc_and_comline();
    return;
}

void
main_init_x(void) {
    char *icon_name = "xpp";
    char *win_name = "XPPAUT";
    int32 x = 0, y = 0;
    uint32 min_wid = 450, min_hgt = 360;

    static uint32 Black;
    static uint32 White;

    char teststr[] = "The Quick Brown Fox Jumped Over The Lazy Dog?";

    if (UserMinWidth > 0) {
        min_wid = (uint32)UserMinWidth;
        SCALEX = (int32)min_wid;
    }

    if (UserMinHeight > 0) {
        min_hgt = (uint32)UserMinHeight;
        SCALEY = (int32)min_hgt;
    }

    if (PaperWhite == 0) {
        GrFore = White;
        GrBack = Black;
    }

    main_win =
        main_init_win(4, icon_name, win_name, x, y, min_wid, min_hgt, 0, NULL);

    /* Set up foreground and background colors */

    Black = (uint32)BlackPixel(display, screen);
    White = (uint32)WhitePixel(display, screen);

    if (strlen(UserBlack) != 0) {
        XColor user_col;

        XParseColor(display, DefaultColormap(display, screen), UserBlack,
                    &user_col);
        XAllocColor(display, DefaultColormap(display, screen), &user_col);

        MyForeColor = GrFore = (uint32)user_col.pixel;
        Black = MyForeColor;
    }

    if (strlen(UserWhite) != 0) {
        XColor user_col;

        XParseColor(display, DefaultColormap(display, screen), UserWhite,
                    &user_col);
        XAllocColor(display, DefaultColormap(display, screen), &user_col);

        MyBackColor = GrBack = (uint32)user_col.pixel;
        White = MyBackColor;
    }

    /*  Switch for reversed video  */
    MyForeColor = GrFore = Black;
    MyBackColor = GrBack = White;

    if (PaperWhite == 1) {
        /* Respect the swapping implied by the -white option. */
        char swapcol[8];
        printf("Doing swap!\n");
        strcpy(swapcol, UserWhite);
        strcpy(UserWhite, UserBlack);
        strcpy(UserBlack, swapcol);

        MyForeColor = GrFore = White;
        MyBackColor = GrBack = Black;
    }

    if (strlen(UserMainWinColor) != 0) {
        XColor main_win_col;

        XParseColor(display, DefaultColormap(display, screen), UserMainWinColor,
                    &main_win_col);
        XAllocColor(display, DefaultColormap(display, screen), &main_win_col);

        MyMainWinColor = (uint32)main_win_col.pixel;
    } else {
        MyMainWinColor = MyBackColor;
    }

    XSetWindowBorder(display, main_win, MyForeColor);
    XSetWindowBackground(display, main_win, MyMainWinColor);

    if (strlen(UserDrawWinColor) != 0) {
        XColor draw_win_col;
        XParseColor(display, DefaultColormap(display, screen), UserDrawWinColor,
                    &draw_win_col);
        XAllocColor(display, DefaultColormap(display, screen), &draw_win_col);

        MyDrawWinColor = (uint32)draw_win_col.pixel;
    } else {
        MyDrawWinColor = MyBackColor;
    }

    main_fix_window_size(main_win, SCALEX, SCALEY, FIX_MIN_SIZE);
    periodic = 1;
    if (DefaultDepth(display, screen) >= 8) {
        COLOR = 1;
    } else {
        COLOR = 0;
    }

    XSelectInput(display, main_win,
                 ExposureMask | KeyPressMask | ButtonPressMask |
                     StructureNotifyMask | ButtonReleaseMask |
                     ButtonMotionMask);

    /* main load fonts */
    if ((big_font = XLoadQueryFont(display, big_font_name)) == NULL) {
        ggets_plintf("X Error: Failed to load big font: %s\n", big_font_name);
        exit(-1);
    }

    if ((small_font = XLoadQueryFont(display, small_font_name)) == NULL) {
        ggets_plintf("X Error: Failed to load small font: %s\n",
                     small_font_name);
        exit(-1);
    }

    for (int32 i = 0; i < 5; i++) {
        if ((symfonts[i] = XLoadQueryFont(display, symbolfonts[i])) == NULL) {
            if (i == 0 || i == 1)
                symfonts[i] = small_font;
            else
                symfonts[i] = big_font;
            avsymfonts[i] = 1;
        } else {
            avsymfonts[i] = 1;
            ggets_plintf(" sym %d loaded ..", i);
        }

        if ((romfonts[i] = XLoadQueryFont(display, timesfonts[i])) == NULL) {
            if (i == 0 || i == 1)
                romfonts[i] = small_font;
            else
                romfonts[i] = big_font;
            avromfonts[i] = 1;
        } else {
            avromfonts[i] = 1;
            ggets_plintf(" times %d loaded ..", i);
        }
    }
    ggets_plintf("\n");

    /* BETTER SUPPORT FOR VARIABLE WIDTH FONTS
     * Use a statistical average to get average spacing. Some fonts don't
     * or are not able to report this accurately so this is reliable way to
     * get the information. If person only has variable width font on their
     * system they can get by.
     * The average spacing will be too small for some int16 strings having
     * capital letters (for example "GO"). Thus, we divide by the string
     * length of our test string minus 2 for a little more wiggle room. */

    DCURXb = XTextWidth(big_font, teststr, (int)strlen(teststr)) /
             (int)(strlen(teststr) - 2);

    DCURYb = big_font->ascent + big_font->descent;
    CURY_OFFb = big_font->ascent - 1;

    DCURXs = XTextWidth(small_font, teststr, (int)strlen(teststr)) /
             (int)(strlen(teststr) - 2);

    DCURYs = small_font->ascent + small_font->descent;
    CURY_OFFs = small_font->ascent - 1;

    main_get_gc(&gc);
    main_get_gc(&gc_graph);
    main_get_gc(&small_gc);
    main_get_gc(&font_gc);

    if (strlen(UserBGBitmap) != 0) {
        uint32 width_return, height_return;
        int32 x_hot, y_hot;
        uchar *pixdata;

        int32 success =
            XReadBitmapFileData(UserBGBitmap, &width_return, &height_return,
                                &pixdata, &x_hot, &y_hot);

        if (success != BitmapSuccess) {
            if (success == BitmapOpenFailed) {
                ggets_plintf(
                    "Problem reading bitmap file %s -> BitmapOpenFailed\n",
                    UserBGBitmap);
            } else if (success == BitmapFileInvalid) {
                ggets_plintf(
                    "Problem reading bitmap file %s -> BitmapFileInvalid\n",
                    UserBGBitmap);
            } else if (success == BitmapNoMemory) {
                ggets_plintf(
                    "Problem reading bitmap file %s -> BitmapNoMemory\n",
                    UserBGBitmap);
            }
        } else {
            Pixmap pmap = XCreatePixmapFromBitmapData(
                display, main_win, (char *)pixdata, width_return, height_return,
                MyForeColor, MyMainWinColor,
                (uint)DefaultDepth(display, DefaultScreen(display)));
            XSetWindowBackgroundPixmap(display, main_win, pmap);
            XFreePixmap(display, pmap);
            XFree(pixdata);
        }
    }

    if (COLOR)
        color_map_make();

    /* main set big font */
    DCURX = DCURXb;
    DCURY = DCURYb;
    CURY_OFF = CURY_OFFb;
    XSetFont(display, gc, big_font->fid);

    XSetFont(display, small_gc, small_font->fid);

    /* If the user didn't specify specifically heights and widths
     * we try to set the initial size to fit everything nicely especially
     * if they are using wacky fonts...  */
    if (UserMinWidth <= 0)
        SCALEX = 10 + 36*2*DCURXs + 32*DCURXs;

    if (UserMinHeight <= 0)
        SCALEY = 25*DCURYb + 7*DCURYs;

    XResizeWindow(display, main_win, (uint)SCALEX, (uint)SCALEY);
    return;
}

/* not sure what to do with this - but it works pretty well!
 * it allows you to create a KB script and send it as
 * fake presses to the X11 event handler

 * special keypresses are
 * #t tab
 * #e escape
 * #b backspace
 * #d delete
 * #r return
 *so for example to change a parameter to some numbert and then
 * run the integration, you would script
 * "piapp#r#b#b#b#b.12#r#rig"

 * p calls parameter prompt
 * iapp#r  types in iapp with a return
 * #b#b#b#b deletes the current value (assuming no more than 4 numbers)
 * .12# types in the number
 * #r gets out of the parameter picker
 * ig  runs XPP

*/
void
main_xpp_events(XEvent report, int32 min_wid, int32 min_hgt) {
    char ch;

    int32 used = 0;

    array_plot_do_events(report);
    txt_view_events(report);
    ani_do_events(report);
    main_top_button_events(report);
    switch (report.type) {
    case ConfigureNotify: /* this needs to be fixed!!! */
        init_conds_resize_par_box(report.xany.window);
        resize_my_browser(report.xany.window);
        eig_list_resize_eq_list(report.xany.window);
        auto_x11_resize_window(report);
        if (report.xconfigure.window == main_win) {
            SCALEX = report.xconfigure.width;
            SCALEY = report.xconfigure.height;
            if ((SCALEX < min_wid) || (SCALEY < min_hgt)) {
                SCALEX = min_wid;
                SCALEY = min_hgt;
            } else {
                XResizeWindow(display, command_pop, (uint)SCALEX - 4,
                              (uint)DCURY + 1);
                XMoveResizeWindow(display, info_pop, 0, SCALEY - DCURY - 4,
                                  (uint)SCALEX - 4, (uint)DCURY);
                init_conds_resize_par_slides(SCALEY - 3*DCURYs - 1*DCURYb -
                                             13);
                many_pops_resize_all(SCALEX, SCALEY);
                main_redraw_all();
            }
        }

        break;
    case Expose:
    case MapNotify:
        if (report.xany.window == command_pop)
            ggets_put_command("Command:");
        many_pops_do_expose(report);

        break;
    case KeyPress:
        used = 0;
        init_conds_box_keypress(report, &used);
        if (used)
            break;
        eig_list_eq_list_keypress(report, &used);
        if (used)
            break;
        my_browse_keypress(report, &used);
        if (used)
            break;
#ifdef AUTO
        auto_x11_keypress(report, &used);
        if (used)
            break;
#endif
        ch = (char)ggets_get_key_press(&report);
        main_commander(ch);

        break;
    case EnterNotify:
        eig_list_enter_eq_stuff(report.xcrossing.window, 2);
        enter_my_browser(report, 1);
        init_conds_enter_slides(report.xcrossing.window, 1);
        init_conds_box_enter_events(report.xcrossing.window, 1);
        menu_crossing(report.xcrossing.window, 1);
#ifdef AUTO
        auto_x11_enter(report.xcrossing.window, 2);
#endif
        break;
    case LeaveNotify:
        eig_list_enter_eq_stuff(report.xcrossing.window, 1);
        enter_my_browser(report, 0);
        init_conds_enter_slides(report.xcrossing.window, 0);
        init_conds_box_enter_events(report.xcrossing.window, 0);
        menu_crossing(report.xcrossing.window, 0);
#ifdef AUTO
        auto_x11_enter(report.xcrossing.window, 1);
#endif
        break;
    case MotionNotify:
        many_pops_do_motion_events(report);
        break;
    case ButtonRelease:
        init_conds_slide_release(report.xbutton.window);

        break;
    case ButtonPress:
        if (!many_pops_rotate_3dcheck(report)) {
            menu_button(report.xbutton.window);
            init_conds_box_buttons(report.xbutton.window);

            init_conds_slide_button_press(report.xbutton.window);
            eig_list_eq_list_button(report);
            my_browse_button(report);
#ifdef AUTO
            auto_x11_button(report);
#endif

            ggets_show_position(report);
        }
        break;
    default:
        break;
    }
    return;
}

void
main_do_events(uint32 min_wid, uint32 min_hgt) {
    XEvent report;

    ggets_blank_screen(main_win);
    menu_help();
    if (RunImmediately == 1) {
        menudrive_run_the_commands(4);
        RunImmediately = 0;
    }
    while (true) {
        XNextEvent(display, &report);
        main_xpp_events(report, (int32)min_wid, (int32)min_hgt);
    }
}

void
main_bye_bye(void) {
    int32 i;
    auto_nox_yes_reset();
    XUnloadFont(display, big_font->fid);
    XUnloadFont(display, small_font->fid);
    for (i = 0; i < 5; i++) {
        if (avsymfonts[i])
            XUnloadFont(display, symfonts[i]->fid);
        if (avromfonts[i])
            XUnloadFont(display, romfonts[i]->fid);
    }
    XFreeGC(display, gc);
    XCloseDisplay(display);
    exit(1);
}

void
main_clr_scrn(void) {
    ggets_blank_screen(draw_win);
    many_pops_restore_off();
    axes2_do();
    return;
}

void
main_redraw_all(void) {
    if (manual_expose == 0) {
        nullcline_redraw_dfield();
        integrate_restore(0, my_browser.maxrow);
        many_pops_draw_label(draw_win);
        graf_par_draw_freeze(draw_win);
        many_pops_restore_on();
    }
    return;
}

void
main_commander(int32 ch) {
    switch (help_menu) {
    case MAIN_MENU: {
        switch (ch) {
        case 'i':
            menudrive_ini_data_menu();
            break;
        case 'c':
            integrate_cont_integ();
            break;
        case 'n':
            menudrive_new_clines();
            break;
        case 'd':
            menudrive_direct_field();
            break;
        case 'w':
            menudrive_window_zoom();
            break;
        case 'a':
            menudrive_do_torus();
            break;
        case 'k':
            menudrive_do_movie();
            break;
        case 'g':
            menudrive_add_a_curve();
            break;
        case 'u':
            menu_help_num();
            break;
        case 'f':
            menu_help_file();
            break;
        case 'p':
            menudrive_new_param();
            break;
        case 'e':
            menudrive_clear_screens();
            break;
        case 'h':
        case 'm':
            menudrive_do_windows();
            break;
        case 't':
            menudrive_do_gr_objs();
            break;
        case 's':
            menudrive_find_equilibrium();
            break;
        case 'v':
            menudrive_change_view();
            break;
        case 'b':
            menudrive_find_bvp();
            break;

        case 'x':
            menudrive_x_vs_t();
            break;
        case 'r':
            menudrive_redraw_them_all();
            break;
        case '3':
            menudrive_get_3d_par();
            break;
        case 'y':
            graphics_draw_many_lines();
            break;
        default:
            break;
        }
        break;
    }
    case NUM_MENU: {
        numerics_get_num_par(ch);
        break;
    }
    case FILE_MENU: {
        switch (ch) {
        case 't':
            adj2_do_transpose();
            break;
        case 'g':
            many_pops_get_intern_set();
            break;
        case 'i':
            TipsFlag = 1 - TipsFlag;
            break;
        case 'p':
            txt_make_view();
            break;
        case 'w':
            do_lunch(0);
            break;
        case 's':
            lunch_file_inf();
            break;
        case 'a':
#ifdef AUTO
            auto_nox_win();
#endif
            break;
        case 'c':
            calc_q_calc();
            break;
        case 'r':
            do_lunch(1);
            break;
        case 'e':
            edit_rhs_menu();
            break;
        case 'b':
            tfBell = 1 - tfBell;
            break;
        case 'h':
            menudrive_xpp_hlp();
            break;
        case 'q':
            if (pop_list_yes_no_box())
                main_bye_bye();
            break;
        case 'l':
            init_conds_clone_ode();
            break;
        case 'x':
            menudrive_edit_xpprc();
            break;
        case 'u':
            menudrive_do_tutorial();
            break;
        default:
            break;
        }
        menu_help();
        break;
    }
    default:
        break;
    }
    return;
}

Window
main_init_win(uint32 bw, char *icon_name, char *win_name, int32 x, int32 y,
              uint32 min_wid, uint32 min_hgt, int32 argc, char **argv) {
    Window wine;
    int32 count;
    int32 dp_h;
    int32 dp_w;
    Pixmap icon_map;
    XIconSize *size_list;
    XSizeHints size_hints;
    char *display_name = NULL;

    if ((display = XOpenDisplay(display_name)) == NULL) {
        ggets_plintf(" Failed to open X-Display \n");
        exit(-1);
    }
    screen = DefaultScreen(display);
    if (!deleteWindowAtom)
        deleteWindowAtom = XInternAtom(display, "WM_DELETE_WINDOW", 0);
    dp_w = DisplayWidth(display, screen);
    dp_h = DisplayHeight(display, screen);
    DisplayWidth = dp_w;
    DisplayHeight = dp_h;
    if (SCALEX > dp_w)
        SCALEX = dp_w;
    if (SCALEY > dp_h)
        SCALEY = dp_h;
    wine = XCreateSimpleWindow(display, RootWindow(display, screen), x, y,
                               (uint)SCALEX, (uint)SCALEY, bw, MyForeColor,
                               MyBackColor);
    XGetIconSizes(display, RootWindow(display, screen), &size_list, &count);
    icon_map = XCreateBitmapFromData(display, wine, (char *)pp_bits, pp_width,
                                     pp_height);

#ifdef X11R3
    size_hints.flags = PPosition | PSize | PMinsize;
    size_hints.x = x;
    size_hints.y = y;
    size_hints.width = width;
    size_hints.height = height;
    size_hints.min_width = min_wid;
    size_hints.min_height = min_hgt;
#else

    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.min_width = (int)min_wid;
    size_hints.min_height = (int)min_hgt;
#endif

#ifdef X11R3
    XSetStandardProperties(display, wine, win_name, icon_name, icon_map, argv,
                           argc, &size_hints);
#else
    {
        XWMHints wm_hints;
        XClassHint class_hints;
        XTextProperty winname, iconname;
        if (XStringListToTextProperty(&icon_name, 1, &iconname) == 0) {
            ggets_plintf("X error: failure for iconname\n");
            exit(-1);
        }
        if (XStringListToTextProperty(&win_name, 1, &winname) == 0) {
            ggets_plintf("X error: failure for winname\n");
            exit(-1);
        }

        wm_hints.initial_state = NormalState;
        wm_hints.input = True;
        wm_hints.icon_pixmap = icon_map;
        wm_hints.flags = StateHint | IconPixmapHint | InputHint;
        class_hints.res_name = "base";
        class_hints.res_class = win_name;

        XSetWMProperties(display, wine, &winname, &iconname, argv, argc,
                         &size_hints, &wm_hints, &class_hints);
        XSetWMProtocols(display, wine, &deleteWindowAtom, 1);
    }
#endif
    return wine;
}

void
main_top_button_draw(Window window) {
    if (window == TopButton[0])
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "ICs  ", 5);
    if (window == TopButton[1])
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "BCs  ", 5);
    if (window == TopButton[2])
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Delay", 5);
    if (window == TopButton[3])
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Param", 5);
    if (window == TopButton[4])
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Eqns ", 5);
    if (window == TopButton[5])
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Data ", 5);
    return;
}

void
main_top_button_cross(Window window, int32 b) {
    int32 i;
    for (i = 0; i < 6; i++)
        if (window == TopButton[i]) {
            XSetWindowBorderWidth(display, window, (uint)b);
            return;
        }
    return;
}

void
main_top_button_events(XEvent report) {
    Window window = report.xbutton.window;
    switch (report.type) {
    case Expose:
    case MapNotify:
        main_top_button_draw(report.xany.window);
        break;
    case EnterNotify:
        main_top_button_cross(report.xcrossing.window, 2);
        break;
    case LeaveNotify:
        main_top_button_cross(report.xcrossing.window, 1);
        break;
    case ButtonPress:
        /* main top button press */
        if (window == TopButton[0])
            init_conds_make_new_ic_box();
        if (window == TopButton[1])
            init_conds_make_new_bc_box();
        if (window == TopButton[2])
            init_conds_make_new_delay_box();
        if (window == TopButton[3])
            init_conds_make_new_param_box();
        if (window == TopButton[4])
            eig_list_create_eq_list();
        if (window == TopButton[5])
            make_new_browser();
        break;
    default:
        break;
    }
    return;
}

void
main_get_gc(GC *gc2) {
    uint32 valuemask = 0;
    XGCValues values;
    *gc2 = XCreateGC(display, main_win, valuemask, &values);
    XSetForeground(display, *gc2, MyForeColor);
    return;
}

void
main_fix_window_size(Window window, int32 width, int32 height, int32 flag) {
    XSizeHints size_hints;
    switch (flag) {
    case FIX_SIZE:
        size_hints.flags = PSize | PMinSize | PMaxSize;
        size_hints.width = width;
        size_hints.min_width = width;
        size_hints.max_width = width;
        size_hints.height = height;
        size_hints.min_height = height;
        size_hints.max_height = height;
        break;
    case FIX_MIN_SIZE:
        size_hints.flags = PMinSize;
        size_hints.min_width = width;

        size_hints.min_height = height;

        break;
    case FIX_MAX_SIZE:
        size_hints.flags = PMaxSize;
        size_hints.max_width = width;
        size_hints.max_height = height;
        break;
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    XSetWMProperties(display, window, NULL, NULL, NULL, 0, &size_hints, NULL,
                     NULL);
    return;
}

int32
main_get_x_colors(XWindowAttributes *win_info, XColor **colors) {
    int32 ncolors;

    *colors = (XColor *)NULL;
    TrueColorFlag = 0;
    if (win_info->visual->class == TrueColor) {
        TrueColorFlag = 1;
        ggets_plintf("TrueColor visual:  no colormap needed\n");
        return 0;
    }

    else if (!win_info->colormap) {
        ggets_plintf("no colormap associated with window\n");
        return 0;
    }

    ncolors = win_info->visual->map_entries;
    ggets_plintf("%d entries in colormap\n", ncolors);

    *colors = xmalloc(sizeof(XColor)*(uint)ncolors);
    xorfix = 0;

    if (win_info->visual->class == DirectColor) {
        int32 red, green, blue, red1, green1, blue1;

        ggets_plintf("DirectColor visual\n");

        red = green = blue = 0;
        red1 = (int32)LOWBIT(win_info->visual->red_mask);
        green1 = (int32)LOWBIT(win_info->visual->green_mask);
        blue1 = (int32)LOWBIT(win_info->visual->blue_mask);
        for (int32 i = 0; i < ncolors; i++) {
            (*colors)[i].pixel = (ulong)(red | green | blue);
            (*colors)[i].pad = 0;
            red += red1;
            if (red > (int32)win_info->visual->red_mask)
                red = 0;
            green += green1;
            if (green > (int32)win_info->visual->green_mask)
                green = 0;
            blue += blue1;
            if (blue > (int32)win_info->visual->blue_mask)
                blue = 0;
        }
    } else {
        for (int32 i = 0; i < ncolors; i++) {
            (*colors)[i].pixel = (ulong)i;
            (*colors)[i].pad = 0;
        }
    }

    XQueryColors(display, win_info->colormap, *colors, ncolors);

    return ncolors;
}

int32
main_get_command_width(void) {
    int32 x;
    int32 y;
    uint32 w, h, bw, de;

    Window root;
    XGetGeometry(display, command_pop, &root, &x, &y, &w, &h, &bw, &de);
    XClearWindow(display, command_pop);
    return (int32)w;
}
