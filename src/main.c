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
#include "xmalloc.h"
#include "integers.h"
#include "pp.bitmap"
#include "read_dir.h"

#define LOWBIT(x) ((x) & (~(x) + 1))

#include <X11/cursorfont.h>

#include "myfonts.h"

int32 all_win_vis = 0;
int32 use_ani_file = 0;
char anifile[XPP_MAX_NAME];

double xpp_version_maj;
double xpp_version_min;

int32 Xup;
int32 flag_tips = 1;
Atom atom_delete_window = 0;
int32 xpp_batch = 0;
int32 batch_range = 0;
int32 batch_equil = -1;
char batch_out[XPP_MAX_NAME];
char user_out_file[XPP_MAX_NAME];
int32 display_height;
int32 display_width;
int32 flag_true_color;
char font_name_big[XPP_MAX_FONT_NAME];
char font_name_small[XPP_MAX_FONT_NAME];
char plot_format[10];

int32 paper_white = -1;

Window draw_win;
Window main_win;
Window command_pop;
Window info_pop;
GC gc;
GC gc_graph;
GC small_gc;
GC font_gc;
char user_black[8];
char user_white[8];
char user_main_win_color[8];
char user_draw_win_color[8];
char user_bg_bitmap[XPP_MAX_NAME];

int32 user_gradients = -1;
int32 user_min_width = 0;
int32 user_min_height = 0;
uint32 my_back_color;
uint32 my_fore_color;
uint32 my_main_win_color;
uint32 my_draw_win_color;
uint32 gr_fore;
uint32 gr_back;
int32 scale_y;
Display *display;
int32 screen;
int32 periodic = 0;
int32 dcur_yb;
int32 cury_offb;
int32 dcur_ys;
int32 dcur_xs;
int32 cury_offs;
int32 dcur_y;
int32 dcur_x;
int32 cury_off;
FILE *logfile;
int32 xpp_verbose = 1;
int32 override_quiet = 0;
int32 override_logfile = 0;
int32 tfBell;

int32 slider1 = -1;
int32 slider2 = -1;
int32 slider3 = -1;
char slider1var[20];
char slider2var[20];
char slider3var[20];
double slider1lo = 0.0;
double slider2lo = 0.0;
double slider3lo = 0.0;
double slider1hi = 1.0;
double slider2hi = 1.0;
double slider3hi = 1.0;

/* Set this to 1 if you want the tutorial to come up at start-up
 * as default behavior */
int32 do_tutorial = 0;

OptionsSet not_already_set;
XFontStruct *font_small;

static int32 dcurx_b;
static int32 SCALEX;
static Window top_button[6];
static XFontStruct *font_big;

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

    if (nit == 0) {
        return;
    }

    dx = (double)icount*(double)cwidth / (double)nit;
    XDrawPoint(display, command_pop, gc, (int32)dx, 5);
    return;
}

int32
main_my_abort(void) {
    int32 ch;

    while (XPending(display) > 0) {
        XEvent event;
        XNextEvent(display, &event);

        if (ani_check_pause(event) == PAUSE_NUMBER) {
            return PAUSE_NUMBER;
        }

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
    uint32 min_wid = 450;
    uint32 min_hgt = 360;
    OptionsSet *tempNS;

    // Track which options have not been set already
    not_already_set.BIG_FONT_NAME = 1;
    not_already_set.SMALL_FONT_NAME = 1;
    not_already_set.BACKGROUND = 1;
    not_already_set.ix_plt = 1;
    not_already_set.iy_plt = 1;
    not_already_set.iz_plt = 1;
    not_already_set.axes = 1;
    not_already_set.nmesh = 1;
    not_already_set.METHOD = 1;
    not_already_set.TIMEPLOT = 1;
    not_already_set.max_stor = 1;
    not_already_set.TEND = 1;
    not_already_set.DT = 1;
    not_already_set.T0 = 1;
    not_already_set.TRANS = 1;
    not_already_set.bound = 1;
    not_already_set.TOLER = 1;
    not_already_set.delay = 1;
    not_already_set.XLO = 1;
    not_already_set.XHI = 1;
    not_already_set.YLO = 1;
    not_already_set.YHI = 1;
    not_already_set.user_black = 1;
    not_already_set.user_white = 1;
    not_already_set.user_main_win_color = 1;
    not_already_set.user_draw_win_color = 1;
    not_already_set.user_gradients = 1;
    not_already_set.user_bg_bitmap = 1;
    not_already_set.user_min_width = 1;
    not_already_set.user_min_height = 1;
    not_already_set.YNullColor = 1;
    not_already_set.XNullColor = 1;
    not_already_set.stable_manifold_color = 1;
    not_already_set.UnstableManifoldColor = 1;
    not_already_set.start_line_type = 1;
    not_already_set.rand_seed = 1;
    not_already_set.paper_white = 1;
    not_already_set.COLORMAP = 1;
    not_already_set.NPLOT = 1;
    not_already_set.DLL_LIB = 1;
    not_already_set.DLL_FUN = 1;
    not_already_set.XP = 1;
    not_already_set.YP = 1;
    not_already_set.ZP = 1;
    not_already_set.NOUT = 1;
    not_already_set.VMAXPTS = 1;
    not_already_set.TOR_PER = 1;
    not_already_set.JAC_EPS = 1;
    not_already_set.NEWT_TOL = 1;
    not_already_set.NEWT_ITER = 1;
    not_already_set.FOLD = 1;
    not_already_set.DTMIN = 1;
    not_already_set.DTMAX = 1;
    not_already_set.ATOL = 1;
    not_already_set.TOL = 1;
    not_already_set.BANDUP = 1;
    not_already_set.BANDLO = 1;
    not_already_set.PHI = 1;
    not_already_set.THETA = 1;
    not_already_set.XMIN = 1;
    not_already_set.XMAX = 1;
    not_already_set.YMIN = 1;
    not_already_set.YMAX = 1;
    not_already_set.ZMIN = 1;
    not_already_set.ZMAX = 1;
    not_already_set.POIVAR = 1;
    not_already_set.OUTPUT = 1;
    not_already_set.POISGN = 1;
    not_already_set.POISTOP = 1;
    not_already_set.STOCH = 1;
    not_already_set.POIPLN = 1;
    not_already_set.POIMAP = 1;
    not_already_set.RANGEOVER = 1;
    not_already_set.RANGESTEP = 1;
    not_already_set.RANGELOW = 1;
    not_already_set.RANGEHIGH = 1;
    not_already_set.RANGERESET = 1;
    not_already_set.RANGEOLDIC = 1;
    not_already_set.RANGE = 1;
    not_already_set.NTST = 1;
    not_already_set.NMAX = 1;
    not_already_set.NPR = 1;
    not_already_set.NCOL = 1;
    not_already_set.DSMIN = 1;
    not_already_set.DSMAX = 1;
    not_already_set.DS = 1;
    not_already_set.PARMAX = 1;
    not_already_set.NORMMIN = 1;
    not_already_set.NORMMAX = 1;
    not_already_set.EPSL = 1;
    not_already_set.EPSU = 1;
    not_already_set.EPSS = 1;
    not_already_set.RUNNOW = 1;
    not_already_set.SEC = 1;
    not_already_set.UEC = 1;
    not_already_set.SPC = 1;
    not_already_set.UPC = 1;
    not_already_set.AUTOEVAL = 1;
    not_already_set.AUTOXMAX = 1;
    not_already_set.AUTOYMAX = 1;
    not_already_set.AUTOXMIN = 1;
    not_already_set.AUTOYMIN = 1;
    not_already_set.AUTOVAR = 1;
    not_already_set.ps_font = 1;
    not_already_set.ps_lw = 1;
    not_already_set.PS_FSIZE = 1;
    not_already_set.PS_COLOR = 1;
    not_already_set.forever = 1;
    not_already_set.bvp_tol = 1;
    not_already_set.bvp_eps = 1;
    not_already_set.bpv_maxit = 1;
    not_already_set.bvp_flag = 1;
    not_already_set.SOS = 1;
    not_already_set.FFT = 1;
    not_already_set.hist = 1;
    not_already_set.plt_fmt_flag = 1;
    not_already_set.atoler = 1;
    not_already_set.euler_max_iter = 1;
    not_already_set.euler_tol = 1;
    not_already_set.evec_iter = 1;
    not_already_set.evec_err = 1;
    not_already_set.NULL_ERR = 1;
    not_already_set.newt_err = 1;
    not_already_set.NULL_HERE = 1;
    not_already_set.TUTORIAL = 1;
    not_already_set.slider1 = 1;
    not_already_set.slider2 = 1;
    not_already_set.slider3 = 1;
    not_already_set.slider1lo = 1;
    not_already_set.slider2lo = 1;
    not_already_set.slider3lo = 1;
    not_already_set.slider1hi = 1;
    not_already_set.slider2hi = 1;
    not_already_set.slider3hi = 1;
    not_already_set.POSTPROCESS = 1;
    not_already_set.HISTCOL = 1;
    not_already_set.HISTLO = 1;
    not_already_set.HISTHI = 1;
    not_already_set.HISTBINS = 1;
    not_already_set.HISTCOL2 = 1;
    not_already_set.HISTLO2 = 1;
    not_already_set.HISTHI2 = 1;
    not_already_set.HISTBINS2 = 1;
    not_already_set.SPECCOL = 1;
    not_already_set.SPECCOL2 = 1;
    not_already_set.SPECWIDTH = 1;
    not_already_set.SPECWIN = 1;
    not_already_set.PLOTFORMAT = 1;
    not_already_set.DFGRID = 1;
    not_already_set.DFBATCH = 1;
    not_already_set.NCBATCH = 1;
    not_already_set.COLORVIA = 1;
    not_already_set.COLORIZE = 1;
    not_already_set.COLORHI = 1;
    not_already_set.COLORLO = 1;

    read_dir_get_directory(myfile);

    SCALEX = 640;
    scale_y = 480;

    Xup = 0;
    snprintf(batch_out, sizeof(batch_out), "output.dat");
    snprintf(plot_format, sizeof(plot_format), "ps");

    /* Read visualization environement variables here
       since some may be overridden by command line */
    logfile = stdout;
    main_check_for_quiet(argc, argv);

    cli_do(argc, argv);

    /* We need to init_X here if there is no file on command line
     * so that a file browser can be opened.  */
    if (!xpp_batch) {
        // Swap out the current options for a temporary place holder
        tempNS = xmalloc(sizeof(OptionsSet));
        *tempNS = not_already_set;
        /* Initialize what's needed to open a browser based on
         * the current options.  */
        main_do_vis_env();
        load_eqn_set_all_vals();
        main_init_x();
        /* Now swap back the options for proper precedence ordering of options.
         */
        not_already_set = *tempNS;
        free(tempNS);
    }

    load_eqn();

    tempNS = xmalloc(sizeof(*tempNS));
    *tempNS = not_already_set;
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

    if (form_ode_idsc(this_file)) {
        METHOD = 0;
    }
    xpp_version_maj = (double)MAJOR_VERSION;
    xpp_version_min = (double)MINOR_VERSION;
    if (strlen(this_file) < 60) {
        snprintf(pptitle, sizeof(pptitle), "XPP Ver %g.%g >> %s", xpp_version_maj, xpp_version_min,
                 this_file);
    } else {
        snprintf(pptitle, sizeof(pptitle), "XPP Version %g.%g", xpp_version_maj, xpp_version_min);
    }
    numerics_do_meth();

    numerics_set_delay();
    rhs_function = main_rhs;
    do_fit_init_info();
    form_ode_strip_saveqn();
    form_ode_create_plot_list();
    extra_auto_load_dll();

    if (xpp_batch) {
        color_map_make();
        browser_init();
        graphics_init_all();
        cli_if_needed_load_set();
        cli_if_needed_load_par();
        cli_if_needed_load_ic();
        cli_if_needed_select_sets();
        cli_if_needed_load_ext_options();
        graphics_set_extra();
        nullcline_set_colorization_stuff();
        integrate_batch();
        if (NCBatch > 0) {
            silent_nullclines();
        }
        if (df_batch > 0) {
            nullcline_silent_dfields();
        }
        integrate_silent_equilibria();
        exit(0);
    }

    many_pops_gtitle_text(pptitle, main_win);
    Xup = 1;
    color_map_make();

    XMapWindow(display, main_win);

    {
        // main make pops
        int32 x;
        int32 y;
        uint32 h;
        uint32 w;
        uint32 bw;
        uint32 d;
        Window wn;
        XGetGeometry(display, main_win, &wn, &x, &y, &w, &h, &bw, &d);
        menu_create_them(main_win);
        command_pop = XCreateSimpleWindow(display, main_win, 0, dcur_ys + 4, w - 2,
                                          (uint)dcur_y + 4, 2, my_fore_color, my_back_color);
        info_pop = XCreateSimpleWindow(display, main_win, 0, (int32)h - dcur_y - 4, w - 2,
                                       (uint)dcur_y, 2, my_fore_color, my_back_color);
        XCreateFontCursor(display, XC_hand2);
        XSelectInput(display, command_pop, KeyPressMask | ButtonPressMask | ExposureMask);
        XSelectInput(display, info_pop, ExposureMask);
        XMapWindow(display, info_pop);
        XMapWindow(display, command_pop);
        many_pops_init_grafs(16*dcur_x + 6, dcur_ys + dcur_yb + 6, (int32)w - 16 - 16*dcur_x,
                             (int32)h - 6*dcur_y - 16);
        init_conds_create_par_sliders(main_win, 0, (int32)h - 5*dcur_y + 8);
        graphics_get_draw_area();
    }

    {
        // main make top buttons
        int32 x1 = 2, x2 = 6*dcur_xs + 5, dx = dcur_xs;
        top_button[0] = pop_list_make_fancy_window(main_win, x1, 1, x2, dcur_ys, 1);
        x1 += x2 + dx;
        top_button[1] = pop_list_make_fancy_window(main_win, x1, 1, x2, dcur_ys, 1);
        x1 += x2 + dx;

        top_button[2] = pop_list_make_fancy_window(main_win, x1, 1, x2, dcur_ys, 1);
        x1 += x2 + dx;

        top_button[3] = pop_list_make_fancy_window(main_win, x1, 1, x2, dcur_ys, 1);
        x1 += x2 + dx;

        top_button[4] = pop_list_make_fancy_window(main_win, x1, 1, x2, dcur_ys, 1);
        x1 += x2 + dx;

        top_button[5] = pop_list_make_fancy_window(main_win, x1, 1, x2, dcur_ys, 1);
        x1 += x2 + dx;
    }

    init_conds_initialize_box();

    browser_init();
    if (all_win_vis == 1) {
        init_conds_make_new_ic_box();
        init_conds_make_new_bc_box();
        init_conds_make_new_delay_box();
        init_conds_make_new_param_box();
        browser_make_new();
        eig_list_create_eq_list();
    }

    Xup = 1;
    ani_zero();
    graphics_set_extra();
    nullcline_set_colorization_stuff();

    pop_list_make_scrbox_lists();

    //          MAIN LOOP
    {
        // main test color info
        XColor *colors;
        XWindowAttributes xwa;
        flag_true_color = 0;

        XGetWindowAttributes(display, main_win, &xwa);
        main_get_x_colors(&xwa, &colors);

        if (colors) {
            free((char *)colors);
        }
    }
    cli_if_needed_load_set();
    cli_if_needed_load_par();
    cli_if_needed_load_ic();
    cli_if_needed_load_ext_options();
    if (use_ani_file) {
        ani_new_vcr();
        ani_get_file(anifile);
    }

    if (do_tutorial == 1) {
        menudrive_do_tutorial();
    }

    graf_par_default_window();

    main_do_events(min_wid, min_hgt);
}

void
main_check_for_quiet(int32 argc, char **argv) {
    // First scan, check for any QUIET option set...
    /* Allow for multiple calls to the QUIET and LOGFILE options
     * on the command line. The last setting is the one that will stick.
     * Settings of logfile and quiet in the xpprc file will be ignored
     * if they are set on the command line.  */
    int32 quiet_specified_once = 0;
    int32 logfile_specified_once = 0;

    for (int32 i = 1; i < argc; i++) {
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
    if (quiet_specified_once == 1) {
        override_quiet = 1;
    }
    if (logfile_specified_once == 1) {
        override_logfile = 1;
    }
    return;
}

void
main_do_vis_env(void) {
    load_eqn_set_x_vals();
    load_eqn_check_for_xpprc();
    load_eqn_set_internopts_xpprc_and_cli();
    return;
}

void
main_init_x(void) {
    char *icon_name = "xpp";
    char *win_name = "XPPAUT";
    int32 x = 0;
    int32 y = 0;
    uint32 min_wid = 450;
    uint32 min_hgt = 360;

    static uint32 Black;
    static uint32 White;

    char teststr[] = "The Quick Brown Fox Jumped Over The Lazy Dog?";

    if (user_min_width > 0) {
        min_wid = (uint32)user_min_width;
        SCALEX = (int32)min_wid;
    }

    if (user_min_height > 0) {
        min_hgt = (uint32)user_min_height;
        scale_y = (int32)min_hgt;
    }

    if (paper_white == 0) {
        gr_fore = White;
        gr_back = Black;
    }

    main_win = main_init_win(4, icon_name, win_name, x, y, min_wid, min_hgt, 0, NULL);

    // Set up foreground and background colors

    Black = (uint32)BlackPixel(display, screen);
    White = (uint32)WhitePixel(display, screen);

    if (strlen(user_black) != 0) {
        XColor user_col;

        XParseColor(display, DefaultColormap(display, screen), user_black, &user_col);
        XAllocColor(display, DefaultColormap(display, screen), &user_col);

        my_fore_color = gr_fore = (uint32)user_col.pixel;
        Black = my_fore_color;
    }

    if (strlen(user_white) != 0) {
        XColor user_col;

        XParseColor(display, DefaultColormap(display, screen), user_white, &user_col);
        XAllocColor(display, DefaultColormap(display, screen), &user_col);

        my_back_color = gr_back = (uint32)user_col.pixel;
        White = my_back_color;
    }

    //  Switch for reversed video
    my_fore_color = gr_fore = Black;
    my_back_color = gr_back = White;

    if (paper_white == 1) {
        // Respect the swapping implied by the -white option.
        char swapcol[8];
        printf("Doing swap!\n");
        strcpy(swapcol, user_white);
        strcpy(user_white, user_black);
        strcpy(user_black, swapcol);

        my_fore_color = gr_fore = White;
        my_back_color = gr_back = Black;
    }

    if (strlen(user_main_win_color) != 0) {
        XColor main_win_col;

        XParseColor(display, DefaultColormap(display, screen), user_main_win_color, &main_win_col);
        XAllocColor(display, DefaultColormap(display, screen), &main_win_col);

        my_main_win_color = (uint32)main_win_col.pixel;
    } else {
        my_main_win_color = my_back_color;
    }

    XSetWindowBorder(display, main_win, my_fore_color);
    XSetWindowBackground(display, main_win, my_main_win_color);

    if (strlen(user_draw_win_color) != 0) {
        XColor draw_win_col;
        XParseColor(display, DefaultColormap(display, screen), user_draw_win_color, &draw_win_col);
        XAllocColor(display, DefaultColormap(display, screen), &draw_win_col);

        my_draw_win_color = (uint32)draw_win_col.pixel;
    } else {
        my_draw_win_color = my_back_color;
    }

    main_fix_window_size(main_win, SCALEX, scale_y, FIX_MIN_SIZE);
    periodic = 1;
    if (DefaultDepth(display, screen) >= 8) {
        COLOR = 1;
    } else {
        COLOR = 0;
    }

    XSelectInput(display, main_win,
                 ExposureMask | KeyPressMask | ButtonPressMask | StructureNotifyMask |
                     ButtonReleaseMask | ButtonMotionMask);

    // main load fonts
    if ((font_big = XLoadQueryFont(display, font_name_big)) == NULL) {
        ggets_plintf("X Error: Failed to load big font: %s\n", font_name_big);
        exit(-1);
    }

    if ((font_small = XLoadQueryFont(display, font_name_small)) == NULL) {
        ggets_plintf("X Error: Failed to load small font: %s\n", font_name_small);
        exit(-1);
    }

    for (int32 i = 0; i < 5; i++) {
        if ((symfonts[i] = XLoadQueryFont(display, symbolfonts[i])) == NULL) {
            if (i == 0 || i == 1) {
                symfonts[i] = font_small;
            } else {
                symfonts[i] = font_big;
            }
            avsymfonts[i] = 1;
        } else {
            avsymfonts[i] = 1;
            ggets_plintf(" sym %d loaded ..", i);
        }

        if ((romfonts[i] = XLoadQueryFont(display, timesfonts[i])) == NULL) {
            if (i == 0 || i == 1) {
                romfonts[i] = font_small;
            } else {
                romfonts[i] = font_big;
            }
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

    dcurx_b = XTextWidth(font_big, teststr, (int)strlen(teststr)) / (int)(strlen(teststr) - 2);

    dcur_yb = font_big->ascent + font_big->descent;
    cury_offb = font_big->ascent - 1;

    dcur_xs = XTextWidth(font_small, teststr, (int)strlen(teststr)) / (int)(strlen(teststr) - 2);

    dcur_ys = font_small->ascent + font_small->descent;
    cury_offs = font_small->ascent - 1;

    main_get_gc(&gc);
    main_get_gc(&gc_graph);
    main_get_gc(&small_gc);
    main_get_gc(&font_gc);

    if (strlen(user_bg_bitmap) != 0) {
        uint32 width_return;
        uint32 height_return;
        int32 x_hot;
        int32 y_hot;
        uchar *pixdata;

        int32 success = XReadBitmapFileData(user_bg_bitmap, &width_return, &height_return, &pixdata,
                                            &x_hot, &y_hot);

        if (success != BitmapSuccess) {
            if (success == BitmapOpenFailed) {
                ggets_plintf("Problem reading bitmap file %s -> BitmapOpenFailed\n",
                             user_bg_bitmap);
            } else if (success == BitmapFileInvalid) {
                ggets_plintf("Problem reading bitmap file %s -> BitmapFileInvalid\n",
                             user_bg_bitmap);
            } else if (success == BitmapNoMemory) {
                ggets_plintf("Problem reading bitmap file %s -> BitmapNoMemory\n", user_bg_bitmap);
            }
        } else {
            Pixmap pmap = XCreatePixmapFromBitmapData(
                display, main_win, (char *)pixdata, width_return, height_return, my_fore_color,
                my_main_win_color, (uint)DefaultDepth(display, DefaultScreen(display)));
            XSetWindowBackgroundPixmap(display, main_win, pmap);
            XFreePixmap(display, pmap);
            XFree(pixdata);
        }
    }

    if (COLOR) {
        color_map_make();
    }

    // main set big font
    dcur_x = dcurx_b;
    dcur_y = dcur_yb;
    cury_off = cury_offb;
    XSetFont(display, gc, font_big->fid);

    XSetFont(display, small_gc, font_small->fid);

    /* If the user didn't specify specifically heights and widths
     * we try to set the initial size to fit everything nicely especially
     * if they are using wacky fonts...  */
    if (user_min_width <= 0) {
        SCALEX = 10 + 36*2*dcur_xs + 32*dcur_xs;
    }

    if (user_min_height <= 0) {
        scale_y = 25*dcur_yb + 7*dcur_ys;
    }

    XResizeWindow(display, main_win, (uint)SCALEX, (uint)scale_y);
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
 * so for example to change a parameter to some numbert and then
 * run the integration, you would script
 * "piapp#r#b#b#b#b.12#r#rig"

 * p calls parameter prompt
 * iapp#r  types in iapp with a return
 * #b#b#b#b deletes the current value (assuming no more than 4 numbers)
 * .12# types in the number
 * #r gets out of the parameter picker
 * ig  runs XPP */
void
main_xpp_events(XEvent report, int32 min_wid, int32 min_hgt) {
    char ch;

    int32 used = 0;

    array_plot_do_events(report);
    txt_view_events(report);
    ani_do_events(report);
    main_top_button_events(report);
    switch (report.type) {
    case ConfigureNotify:  // this needs to be fixed!!!
        init_conds_resize_par_box(report.xany.window);
        browser_my_resize(report.xany.window);
        eig_list_resize_eq_list(report.xany.window);
        auto_x11_resize_window(report);
        if (report.xconfigure.window == main_win) {
            SCALEX = report.xconfigure.width;
            scale_y = report.xconfigure.height;
            if ((SCALEX < min_wid) || (scale_y < min_hgt)) {
                SCALEX = min_wid;
                scale_y = min_hgt;
            } else {
                XResizeWindow(display, command_pop, (uint)SCALEX - 4, (uint)dcur_y + 1);
                XMoveResizeWindow(display, info_pop, 0, scale_y - dcur_y - 4, (uint)SCALEX - 4,
                                  (uint)dcur_y);
                init_conds_resize_par_slides(scale_y - 3*dcur_ys - 1*dcur_yb - 13);
                many_pops_resize_all(SCALEX, scale_y);
                main_redraw_all();
            }
        }

        break;
    case Expose:
    case MapNotify:
        if (report.xany.window == command_pop) {
            ggets_put_command("Command:");
        }
        many_pops_do_expose(report);

        break;
    case KeyPress:
        used = 0;
        init_conds_box_keypress(report, &used);
        if (used) {
            break;
        }
        eig_list_eq_list_keypress(report, &used);
        if (used) {
            break;
        }
        browser_my_keypress(report, &used);
        if (used) {
            break;
        }
#ifdef AUTO
        auto_x11_keypress(report, &used);
        if (used) {
            break;
        }
#endif
        ch = (char)ggets_get_key_press(&report);
        main_commander(ch);

        break;
    case EnterNotify:
        eig_list_enter_eq_stuff(report.xcrossing.window, 2);
        browser_my_enter(report, 1);
        init_conds_enter_slides(report.xcrossing.window, 1);
        init_conds_box_enter_events(report.xcrossing.window, 1);
        menu_crossing(report.xcrossing.window, 1);
#ifdef AUTO
        auto_x11_enter(report.xcrossing.window, 2);
#endif
        break;
    case LeaveNotify:
        eig_list_enter_eq_stuff(report.xcrossing.window, 1);
        browser_my_enter(report, 0);
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
            browser_my_button(report);
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
    if (run_immediately == 1) {
        menudrive_run_the_commands(4);
        run_immediately = 0;
    }
    while (true) {
        XNextEvent(display, &report);
        main_xpp_events(report, (int32)min_wid, (int32)min_hgt);
    }
}

void
main_bye_bye(void) {
    auto_nox_yes_reset();
    XUnloadFont(display, font_big->fid);
    XUnloadFont(display, font_small->fid);

    for (int32 i = 0; i < 5; i++) {
        if (avsymfonts[i]) {
            XUnloadFont(display, symfonts[i]->fid);
        }
        if (avromfonts[i]) {
            XUnloadFont(display, romfonts[i]->fid);
        }
    }

    XFreeGC(display, gc);
    XCloseDisplay(display);
    exit(1);
}

void
main_clr_scrn(void) {
    ggets_blank_screen(draw_win);
    many_pops_restore_off();
    axes_do();
    return;
}

void
main_redraw_all(void) {
    if (manual_expose == 0) {
        nullcline_redraw_dfield();
        integrate_restore(0, browser_my.maxrow);
        many_pops_draw_label(draw_win);
        graf_par_draw_freeze(draw_win);
        many_pops_restore_on();
    }
    return;
}

void
main_commander(int32 ch) {
    switch (help_menu) {
    case MENU_MAIN: {
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
    case MENU_NUM: {
        numerics_get_num_par(ch);
        break;
    }
    case MENU_FILE: {
        switch (ch) {
        case 't':
            adjoints_do_transpose();
            break;
        case 'g':
            many_pops_get_intern_set();
            break;
        case 'i':
            flag_tips = 1 - flag_tips;
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
            if (pop_list_yes_no_box()) {
                main_bye_bye();
            }
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
main_init_win(uint32 bw, char *icon_name, char *win_name, int32 x, int32 y, uint32 min_wid,
              uint32 min_hgt, int32 argc, char **argv) {
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
    if (!atom_delete_window) {
        atom_delete_window = XInternAtom(display, "WM_DELETE_WINDOW", 0);
    }
    dp_w = DisplayWidth(display, screen);
    dp_h = DisplayHeight(display, screen);
    display_width = dp_w;
    display_height = dp_h;
    if (SCALEX > dp_w) {
        SCALEX = dp_w;
    }
    if (scale_y > dp_h) {
        scale_y = dp_h;
    }
    wine = XCreateSimpleWindow(display, RootWindow(display, screen), x, y, (uint)SCALEX,
                               (uint)scale_y, bw, my_fore_color, my_back_color);
    XGetIconSizes(display, RootWindow(display, screen), &size_list, &count);
    icon_map = XCreateBitmapFromData(display, wine, (char *)pp_bits, pp_width, pp_height);

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
    XSetStandardProperties(display, wine, win_name, icon_name, icon_map, argv, argc, &size_hints);
#else
    {
        XWMHints wm_hints;
        XClassHint class_hints;
        XTextProperty winname;
        XTextProperty iconname;
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

        XSetWMProperties(display, wine, &winname, &iconname, argv, argc, &size_hints, &wm_hints,
                         &class_hints);
        XSetWMProtocols(display, wine, &atom_delete_window, 1);
    }
#endif
    return wine;
}

void
main_top_button_draw(Window window) {
    if (window == top_button[0]) {
        XDrawString(display, window, small_gc, 5, cury_offs, "ICs  ", 5);
    }
    if (window == top_button[1]) {
        XDrawString(display, window, small_gc, 5, cury_offs, "BCs  ", 5);
    }
    if (window == top_button[2]) {
        XDrawString(display, window, small_gc, 5, cury_offs, "Delay", 5);
    }
    if (window == top_button[3]) {
        XDrawString(display, window, small_gc, 5, cury_offs, "Param", 5);
    }
    if (window == top_button[4]) {
        XDrawString(display, window, small_gc, 5, cury_offs, "Eqns ", 5);
    }
    if (window == top_button[5]) {
        XDrawString(display, window, small_gc, 5, cury_offs, "Data ", 5);
    }
    return;
}

void
main_top_button_cross(Window window, int32 b) {
    for (int32 i = 0; i < 6; i++) {
        if (window == top_button[i]) {
            XSetWindowBorderWidth(display, window, (uint)b);
            return;
        }
    }
    return;
}

void
main_top_button_events(XEvent report) {
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
    case ButtonPress: {
        Window window = report.xbutton.window;
        // main top button press
        if (window == top_button[0]) {
            init_conds_make_new_ic_box();
        }
        if (window == top_button[1]) {
            init_conds_make_new_bc_box();
        }
        if (window == top_button[2]) {
            init_conds_make_new_delay_box();
        }
        if (window == top_button[3]) {
            init_conds_make_new_param_box();
        }
        if (window == top_button[4]) {
            eig_list_create_eq_list();
        }
        if (window == top_button[5]) {
            browser_make_new();
        }
        break;
    }
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
    XSetForeground(display, *gc2, my_fore_color);
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
    XSetWMProperties(display, window, NULL, NULL, NULL, 0, &size_hints, NULL, NULL);
    return;
}

int32
main_get_x_colors(XWindowAttributes *win_info, XColor **colors) {
    int32 ncolors;

    *colors = (XColor *)NULL;
    flag_true_color = 0;
    if (win_info->visual->class == TrueColor) {
        flag_true_color = 1;
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
        int32 red;
        int32 green;
        int32 blue;
        int32 red1;
        int32 green1;
        int32 blue1;

        ggets_plintf("DirectColor visual\n");

        red = green = blue = 0;
        red1 = (int32)LOWBIT(win_info->visual->red_mask);
        green1 = (int32)LOWBIT(win_info->visual->green_mask);
        blue1 = (int32)LOWBIT(win_info->visual->blue_mask);
        for (int32 i = 0; i < ncolors; i++) {
            (*colors)[i].pixel = (ulong)(red | green | blue);
            (*colors)[i].pad = 0;
            red += red1;
            if (red > (int32)win_info->visual->red_mask) {
                red = 0;
            }
            green += green1;
            if (green > (int32)win_info->visual->green_mask) {
                green = 0;
            }
            blue += blue1;
            if (blue > (int32)win_info->visual->blue_mask) {
                blue = 0;
            }
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
    uint32 w;
    uint32 h;
    uint32 bw;
    uint32 de;

    Window root;
    XGetGeometry(display, command_pop, &root, &x, &y, &w, &h, &bw, &de);
    XClearWindow(display, command_pop);
    return (int32)w;
}
