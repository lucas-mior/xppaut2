#include "functions.h"
#include "integers.h"

#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>

int32 help_menu;

/* end CLONE change */
static struct MenuDef {
    Window base;
    Window title;
    Window window[25];
    char key[25];
    char **names;
    char **hints;
    int32 n;
    int32 visible;
} my_menus[3];

#include <X11/Xlib.h>

#define MAIN_ENTRIES 20
/* CLONE */
#define FILE_ENTRIES 16
#define NUM_ENTRIES 18
static char *main_menu[] = {
    "XPP",        "Initialconds", "Continue",      "Nullcline", "Dir.field/flow", "Window/zoom",
    "phAsespace", "Kinescope",    "Graphic stuff", "nUmerics",  "File",           "Parameters",
    "Erase",      "Makewindow",   "Text,etc",      "Sing pts",  "Viewaxes",       "Xi vs t",
    "Restore",    "3d-params",    "Bndryval"};

static char *num_menu[] = {
    "NUMERICS",    "Total",  "Start time", "tRansient", "Dt",         "Ncline ctrl", "sIng pt ctrl",
    "nOutput",     "Bounds", "Method",     "dElay",     "Color code", "stocHast",    "Poincare map",
    "rUelle plot", "looKup", "bndVal",     "Averaging", "[Esc]-exit"};
/* CLONE change */

static char *fileon_menu[] = {"FILE",       "Prt src",   "Write set", "Read set",    "Auto",
                              "Calculator", "Edit",      "Save info", "Bell off",    "Help",
                              "Quit",       "Transpose", "tIps",      "Get par set", "cLone",
                              ".Xpprc",     "tUtorial"};

static char *fileoff_menu[] = {"FILE",       "Prt src",   "Write set", "Read set",    "Auto",
                               "Calculator", "Edit",      "Save info", "Bell on",     "Help",
                               "Quit",       "Transpose", "tIps",      "Get par set", "cLone",
                               ".Xpprc",     "tUtorial"};

/* hints for the main menus */
static char *main_hint[] = {"Integrate the equations",
                            "Continue integration for specified time",
                            "Draw nullclines",
                            "Direction fields and flows of the phaseplane",
                            "Change the size of two-dimensional view",
                            "Set up periodic/torus phase space",
                            "Take snapshots of the screen",
                            "Adding graphs,hard copy, etc",
                            "Numerics options",
                            "Quit, save stuff, etc",
                            "Change problem parameters",
                            "Clear screen",
                            "Create other windows",
                            "Add fancy text and lines,arrows",
                            "Find fixed points and stability",
                            "Change 2 or 3d views",
                            "Plot variable vs time",
                            "Redraw the graph ",
                            "Set parameters for 3D view",
                            "Run boundary value solver"};

static char *file_hint[] = {"Display source and active comments",
                            "Save information for restart",
                            "Read information for restart",
                            "Run AUTO, the bifurcation package",
                            "A little calculator -- press ESC to exit",
                            "Edit right-hand sides or functions or auxiliaries",
                            "Save info about simulation in human readable format",
                            "Turn bell on/off",
                            "Browser help",
                            "Duh!",
                            "Transpose storage",
                            "Turn off these silly tips",
                            "Set predefined parameters",
                            "Clone the ode file",
                            "Edit your .xpprc preferences file",
                            "Run a quick tutorial on XPPAUT"};

static char *num_hint[] = {"Total time to integrate eqns",
                           "Starting time -- T0",
                           "Time to integrate before storing",
                           "Time step to use",
                           "Mesh for nullclines",
                           "Numerical parameters for fixed points",
                           "Number of steps per plotted point",
                           "Maximum allowed size of any variable",
                           "Integration method",
                           "Maximum delay and delay related stuff",
                           "Color trajectories according to velocity,etc",
                           "Curve fitting, FFT, mean, variance, seed, etc",
                           "Define Poincare map parameters",
                           "Define shifted plots",
                           "Modify lookup tables",
                           "Numerical setup for boundary value solver",
                           "Compute adjoint and averaged functions",
                           "Return to main menu"};

/* other hints  */

char *null_hint[] = {"Compute new nullclines",
                     "Redraw last nullclines",
                     "Set automatic redraw -- X redraws when needed",
                     "Only redraw when asked (Redraw)",
                     "Freeze multiple nullclines",
                     "Save nullcline values to a file"};

char *null_freeze[] = {"Freeze current clines", "Delete all frozen clines",
                       "Range freeze a bunch of clines", "Animate nullclines"};
char *ic_hint[] = {"Integrate over a range of parameters, init data, etc",
                   "Integrate over range of 2 parameters,init data, etc",
                   "Pick up from last step of previous solution",
                   "Use current initial data",
                   "Use current initial data",
                   "Specify initial data with mouse",
                   "Pick up from last step and shift T to latest time",
                   "Input new initial data",
                   "Use the initial data from last shooting from fixed pt",
                   "Read init data from file",
                   "Type in function of 't' for t^th equation",
                   "Repeated ICs using mouse ",
                   "Guess new values for DAE",
                   "Integrate backwards"};

char *wind_hint[] = {"Manually choose 2D view", "Zoom into with mouse",
                     "Zoom out with mouse",     "Let XPP automatically choose window",
                     "Reset to default view",   "Scroll around the view"};

char *flow_hint[] = {"Draw vector field for 2D section", "Draw regular series of trajectories", " ",
                     "Color the PP on a grid", "Draw only directions"};

char *phas_hint[] = {"Each variable is on a circle", "No variable on circle",
                     "Choose circle variables"};

char *kin_hint[] = {"Grab a screen shot",
                    "Clear all screen shots",
                    "Manually cycle thru screenshots",
                    "Continuously cycle through screen shots",
                    "Dump the screen shots to disk",
                    "Make animated gif file from screenshots"};

char *graf_hint[] = {"Add another curve to the current plot",
                     "Delete last added plot",
                     "Remove all the added plots except the main one",
                     "Edit parameters for a plot",
                     "Create a postscript file of current plot",
                     "Create a styleable svg file of current plot",
                     "Options for permanently saving curve",
                     "Axes label sizes and zero axes for postscript",
                     "Export the numbers used in the graphs on the screen",
                     "Change colormap"};

char *cmap_hint[] = {
    " blue-green-red", "red-...-violet-red", "black-red-yellow-white",     "pale cyan->pale yellow",
    "blue-violet-red", "black-white",        "helical luminence corrected"};
char *frz_hint[] = {
    "Permanently keep main curve -- even after reintegrating",
    "Delete specified frozen curve",
    "Edit specified frozen curve",
    "Remove all frozen curves",
    "Toggle key on/off",
    "Import bifurcation data",
    "Clear imported bifurcation curve",
    "Automatically freeze after each integration",
};

char *stoch_hint[] = {"Seed random number generator",
                      "Run many simulations to get average trajectory",
                      "Get data from last simulation",
                      "Load average trajectory from many runs",
                      "Load variance of from many runs",
                      "Compute histogram",
                      "Reload last histogram",
                      "Fourier series of trajectory",
                      "Power spectrum/phase",
                      "Fit data from file to trajectory",
                      "Get mean/variance of single trajectory",
                      "Compute maximal Liapunov exponent",
                      "Compute spike-time autocorrel",
                      "Compute correlations - subtracting mean",
                      "Compute windowed spectral density",
                      "Compute two-variable histograms"};

char *bvp_hint[] = {"Solve BVP over range of parameters", "Don't show any but final step",
                    "Show each step of iteration", "Solve BVP with periodic conditions",
                    "Set up special homoclinic stuff"};

char *adj_hint[] = {"Compute a new adjoint function",
                    "Compute averaging interaction function",
                    "Load computed adjoint",
                    "Load computed orbit",
                    "Load computed interaction function",
                    "Adjoint numerical parameters",
                    "Range over stuff to computte many adjoints"};

char *map_hint[] = {
    "Turn off Poincare map",
    "Define section for Poincare map",
    "Compute Poincare map on maximum/minimum of variable",
    "Compute map based on period between events",
};

char *view_hint[] = {"Two-dimensional view settings", "Three-dimensional view settings",
                     "Plot array ", "Animation window"};

char *half_hint[] = {"Create new window",
                     "Delete all but main window",
                     "Delete last window",
                     "Place current window at bottom",
                     "Automatically redraw",
                     "Redraw only when requested",
                     "Plot all graphs simultaneously -- slows you down"};

char *text_hint[] = {"Create text labels in different fonts ",
                     "Add arrows to trajectories",
                     "Create lines with arrowheads",
                     "Add squares, circles, etc to plot",
                     "Change properties, delete, or move existing add-on",
                     "Delete all text, markers, arrows",
                     "Create many markers based on browser data"};

char *edit_hint[] = {"Move the selected item", "Change properties of selected item",
                     "Delete selected item"};

char *sing_hint[] = {"Find fixed points over range of parameter", " ",
                     "Use mouse to guess fixed point", "Monte carlo search for fixed points"};

char *meth_hint[] = {"Discrete time -- difference equations",
                     "Euler method",
                     "Heun method -- 2nd order Euler",
                     "4th order Runge-Kutta",
                     "4th order predictor-corrector",
                     "Gear's method -- for stiff systems",
                     "Integrator for Volterra equations",
                     "Implicit Euler scheme",
                     "4th order Runge-Kutta with adaptive steps <- RECOMMENDED",
                     "Another stiff method if Gear fails",
                     "Stiff - bad for discontinuous systems",
                     "Dormand-Prince5",
                     "Dormand-Prince8(3)",
                     "Rosenbrock(2,3) - good with discontinuties",
                     "Symplectic - x''=F(x)"};

char *color_hint[] = {" ", "Color according to magnitude of derivative",
                      "Color according to height of Z-axis"};

char *tab_hint[] = {"Edit the lookup tables", "View a table in the data browser"};

char *edrh_hint[] = {"Edit right-hand sides and auxiliaries", "Edit function definitions",
                     "Save current file with new defs", "Load external C right-hand sides"};

char *auto_hint[] = {"Tell AUTO the parameters you may vary",
                     "What will be plotted on the axes and what parameter(s)",
                     "Tell AUTO range, direction, and tolerance",
                     "Run the continuation",
                     "Grab a point to continue from",
                     "Tell AUTO to save at values of period or parameter",
                     "Erase the screen",
                     "Redraw the diagram",
                     "Save and output options"};

char *no_hint[] = {" ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " "};

char *aaxes_hint[] = {"Plot maximum of variable vs parameter",
                      "Plot norm of solution vs parameter",
                      "Plot max/min of variable vs parameter",
                      "Plot period of orbit vs parameter",
                      "Set up two-parameter bifurcation",
                      "Zoom in with mouse",
                      "Zoom out with mouse",
                      "Recall last 1 parameter plot",
                      "Recall last 2 parameter plot",
                      "Let XPP determine the bounds of the plot",
                      "Plot frequency vs parameter",
                      "Plot average of orbit vs parameter",
                      "Return to default bounds of the plot",
                      "Scroll around the plot"};

char *afile_hint[] = {
    "Load a computed orbit into XPP",
    "Write diagram info to file for reuse",
    "Load previously saved file for restart",
    "Create postscript file of picture",
    "Create SVG file of picture",
    "Delete all points of diagram and associated files",
    "Clear grab point to allow start from new point",
    "Write the x-y values of the current diagram to file",
    "Write all the info for the whole diagram!",
    "Save initial data for whole diagram",
    "Toggle automatic redraw",
    "Range over a marked branch",
    "Select a point in 2 parameter diagram",
    "Draw orbits of labeled points automatically",
    "Put all data from branch into browser",
};

char *aspecial_hint[] = {
    "Bifurcation or branch point",
    "Endpoint of a branch",
    "Hopf bifurcation point",
    "Limit point or turning point of a branch",
    "Failure to converge",
    "Period doubling bifurcation",
    "Torus bifurcation from a periodic",
    "User defined function",
};

char *arun_hint[] = {
    "Start at fixed point",
    "Start at periodic orbit",
    "Start at solution to boundary value problem",
    "Start at homoclinics",
    "Start at a heteroclinic",
};

char *browse_hint[] = {"Find closest data point to given value",
                       "Scroll up",
                       "Scroll down",
                       "Scroll up a page",
                       "Scroll down a page",
                       "Scroll left",
                       "Scroll right",
                       "First plotted point",
                       "Last plotted point",
                       "Mark first point for plotting",
                       "Mark last point for plotting",
                       "Redraw data",
                       "Write data to ascii file",
                       "Load first line of Browser to initial data",
                       "Replace column by formula",
                       "Unreplace last replacement",
                       "Write a column of data in tabular format",
                       "Load data from a file into Browser",
                       " ",
                       "Add a new column to Browser",
                       "Delete a column from Browser"};

static void unshow_menu(int32 j);
static void show_menu(int32 j);

void
menu_add(Window base, int32 j, int32 n, char **names, char *key, char **hint) {
    Window window;
    Cursor cursor;
    cursor = XCreateFontCursor(display, XC_hand2);
    window = pop_list_make_plain_unmapped_window(base, 0, dcur_ys + dcur_yb + 10, 16*dcur_x,
                                                 21*(dcur_y + 2) - 3, 1);
    my_menus[j].base = window;
    XDefineCursor(display, window, cursor);
    my_menus[j].names = names;
    my_menus[j].n = n;
    my_menus[j].hints = hint;
    strcpy(my_menus[j].key, key);
    my_menus[j].title = pop_list_make_unmapped_window(window, 0, 0, 16*dcur_x, dcur_y, 1);
    for (int32 i = 0; i < n; i++) {
        my_menus[j].window[i] = pop_list_make_unmapped_window(window, 0, (i + 1)*(dcur_y + 2),
                                                              16*dcur_x, dcur_y, 0);
    }
    my_menus[j].visible = 0;
    XMapRaised(display, my_menus[j].base);
    XMapSubwindows(display, my_menus[j].base);
    return;
}

void
menu_create_them(Window base) {
    char key[30];
    strcpy(key, "icndwakgufpemtsvxr3b");
    menu_add(base, MENU_MAIN, MAIN_ENTRIES, main_menu, key, main_hint);
    strcpy(key, "tsrdniobmechpukva");
    key[17] = 27;
    key[18] = 0;
    menu_add(base, MENU_NUM, NUM_ENTRIES, num_menu, key, num_hint);
    // CLONE
    strcpy(key, "pwracesbhqtiglxu");
    menu_add(base, MENU_FILE, FILE_ENTRIES, fileon_menu, key, file_hint);
    help_menu = -1;
    return;
}

void
show_menu(int32 j) {
    XRaiseWindow(display, my_menus[j].base);
    my_menus[j].visible = 1;
    help_menu = j;
    return;
}

void
unshow_menu(int32 j) {
    if (j < 0) {
        return;
    }
    my_menus[j].visible = 0;
    return;
}

void
menu_help(void) {
    unshow_menu(help_menu);
    show_menu(MENU_MAIN);
    return;
}

void
menu_help_num(void) {
    unshow_menu(help_menu);
    show_menu(MENU_NUM);
    return;
}

void
menu_help_file(void) {
    if (tfBell) {
        my_menus[MENU_FILE].names = fileon_menu;
    } else {
        my_menus[MENU_FILE].names = fileoff_menu;
    }
    unshow_menu(help_menu);
    show_menu(MENU_FILE);
    return;
}

void
menu_crossing(Window win, int32 yn) {
    int32 n;
    int32 j = help_menu;
    char **z;
    if (j < 0) {
        return;
    }
    if (my_menus[j].visible == 0) {
        return;
    }
    n = my_menus[j].n;
    z = my_menus[j].hints;
    for (int32 i = 0; i < n; i++) {
        if (win == my_menus[j].window[i]) {
            XSetWindowBorderWidth(display, win, (uint)yn);
            if (yn && flag_tips) {
                ggets_bottom_msg(z[i]);
            }
            return;
        }
    }
    return;
}

void
menu_expose(Window win) {
    int32 n;
    int32 j = help_menu;
    char **z;
    if (j < 0) {
        return;
    }
    if (my_menus[j].visible == 0) {
        return;
    }
    n = my_menus[j].n;
    z = my_menus[j].names;
    if (win == my_menus[j].title) {
        ggets_set_fore();
        ggets_bar(0, 0, 16*dcur_x, dcur_y, win);
        ggets_set_back();
        XDrawString(display, win, gc, dcur_x / 2 + 5, cury_off, z[0], (int)strlen(z[0]));
        ggets_set_fore();
        return;
    }
    for (int32 i = 0; i < n; i++) {
        if (win == my_menus[j].window[i]) {
            many_pops_base_col();
            XDrawString(display, win, gc, 5, cury_off, z[i + 1], (int)strlen(z[i + 1]));
            return;
        }
    }
    return;
}

void
menu_button(Window win) {
    int32 n;
    int32 j = help_menu;

    if (j < 0) {
        return;
    }
    if (my_menus[j].visible == 0) {
        return;
    }
    n = my_menus[j].n;
    for (int32 i = 0; i < n; i++) {
        if (win == my_menus[j].window[i]) {
            XSetWindowBorderWidth(display, win, 0);
            main_commander(my_menus[j].key[i]);
            return;
        }
    }
    return;
}

void
menu_draw_help(void) {
    int32 j = help_menu;
    int32 n;
    if (j < 0) {
        return;
    }
    if (my_menus[j].visible == 0) {
        return;
    }
    n = my_menus[j].n;
    menu_expose(my_menus[j].title);
    for (int32 i = 0; i < n; i++) {
        menu_expose(my_menus[j].window[i]);
    }
}
