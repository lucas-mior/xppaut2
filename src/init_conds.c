#include "functions.h"
#include "parserslow.h"
#include "integers.h"
#include <stdbool.h>

#include <stdlib.h>
#include <string.h>
#include <strings.h>

/*    This makes a big box with windows that have the names of the
       variables and their current initial data, parameters, BCs
       etc

         This also has the slider boxes

        This also has the file selector gadget
This also has the clone gadget
*/

#include <dirent.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#include <time.h>
#include "xpplim.h"
#include "ic.bitmap"
#include "param.bitmap"
#include "delay.bitmap"

#include "filebrowse.bitmap"
#include "pageup.bitmap"
#include "pagedn.bitmap"
#include "lineup.bitmap"
#include "linedn.bitmap"
#include "home.bitmap"
#include "start.bitmap"

#include "bc.bitmap"

#include "read_dir.h"

#include "mykeydef.h"
#define HOTWILD 2
#define HOTFILE 1

#define READEM 1
#define WRITEM 0

/*extern char UserBGBitmap[100];*/
#define PARAMBOX 1
#define ICBOX 2
#define DELAYBOX 3
#define BCBOX 4
#define BOXEVENT                                                               \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     LeaveWindowMask | EnterWindowMask)

#define EDIT_WAIT 0
#define EDIT_NEXT 1
#define EDIT_ESC 2
#define EDIT_DONE 3

#define FILESELNWIN 10
static struct FileSel {
    int32 n, n0, here;
    Window base, cancel, ok, up, dn, pgup, pgdn, file, wild;
    Window window[FILESELNWIN];
    Window dir, home, start;
    Window fw;
    Window ww;
    char wildtxt[256];
    char filetxt[256];
    int32 nwin, minwid, minhgt;
    int32 off, pos, hot;
    char title[256];
} filesel;

static void display_file_sel(struct FileSel f, Window window);

static struct ParSlider {
    int32 use, pos, l;
    char parname[20];
    double lo, hi, val;
    int32 hgt;
    int32 type;
    int32 index;
    Window left, right, top, main, slide, go;
} my_par_slide[3];

static void do_slide_button(Window window, struct ParSlider *p);
static void redraw_slide(struct ParSlider *p);
static void set_slide_pos(struct ParSlider *p);
static void do_slide_release(Window window, struct ParSlider *p);
static void do_slide_motion(Window window, int32 x, struct ParSlider *p,
                            int32 state);
static void enter_slider(Window window, struct ParSlider *p, int32 val);
static void expose_slider(Window window, struct ParSlider *p);
static void load_entire_box(BoxList *b);
static void set_value_from_box(BoxList *b, int32 i);
static int32 to_float(char *s, double *z);
static void add_editval(BoxList *b, int32 i, char *string);
static void add_edit_float(BoxList *b, int32 i, double z);
static int32 edit_bitem(BoxList *b, int32 i, int32 ch);
static void put_edit_cursor(Window window, int32 pos);
static void draw_editable(Window win, char *string, int32 off, int32 cursor,
                          int32 mc);
static void set_default_params(void);
static void do_box_key(BoxList *b, XEvent event, int32 *used);
static void box_list_scroll(BoxList *b, int32 i);
static void do_box_button(BoxList *b, Window window);
static void redraw_entire_box(BoxList *b);
static void box_enter(BoxList b, Window window, int32 val);
static void display_box(BoxList b, Window window);
static void justify_string(Window w1, char *s1);
static void make_box_list_window(BoxList *b, int32 type);
static void destroy_box(BoxList *b);
static void get_nrow_from_hgt(int32 h, int32 *n, int32 *w);
static void make_par_slider(Window base, int32 x, int32 y, int32 width,
                            int32 index);
static void draw_slider(Window window, int32 x, int32 hgt, int32 l);
static int32 selector_key(XEvent event);
static void stringintersect(char *target, char *sother);
static void create_file_selector(char *title, char *file, char *wild);
static void crossing_selector(Window window, int32 c);
static int32 button_selector(Window window);
static void fs_scroll(int32 i);
static void redraw_fs_text(char *string, Window window, int32 flag);
static void redraw_file_list(void);
static void redraw_directory(void);
static void expose_selector(Window window);

static BoxList *HotBox;
static int32 HotBoxItem = -1;
static BoxList ICBox;
BoxList ParamBox;
static BoxList DelayBox;
static BoxList BCBox;

/* CLONE */
void
init_conds_clone_ode(void) {
    int32 i, j, x, y;
    FILE *fp;
    char clone[256];

    char *s;
    time_t ttt;
    double z;
    clone[0] = 0;
    if (!init_conds_file_selector("Clone ODE file", clone, "*.ode"))
        return;
    if ((fp = fopen(clone, "w")) == NULL) {
        ggets_err_msg(" Cant open clone file");
        return;
    }
    ttt = time(0);
    fprintf(fp, "# clone of %s on %s", this_file, ctime(&ttt));
    for (i = 0; i < NLINES; i++) {
        s = save_eqn[i];

        if (s[0] == 'p' || s[0] == 'P' || s[0] == 'b' || s[0] == 'B') {
            x = form_ode_find_char(s, "'", 0, &j);
            y = form_ode_find_char(s, "=", 0, &j);

            if (x != 0 || y != 0) {
                fprintf(fp, "# original\n# %s\n", s);
                continue;
            }
        }
        if (strncasecmp("done", s, 4) == 0)
            continue;
        fprintf(fp, "%s\n", s);
    }
    fprintf(fp, "# Cloned parameters etc here\n");
    /* now we do parameters boundary conds and ICs */
    j = 0;
    fprintf(fp, "init ");
    for (i = 0; i < (NODE + NMarkov); i++) {
        if (j == 8) {
            fprintf(fp, "\ninit ");
            j = 0;
        }

        fprintf(fp, "%s=%g ", uvar_names[i], last_ic[i]);
        j++;
    }
    fprintf(fp, "\n");

    /* BDRY conds */
    if (my_bc[0].string[0] != '0') {
        for (i = 0; i < NODE; i++)
            fprintf(fp, "bdry %s\n", my_bc[i].string);
    }
    j = 0;
    if (NUPAR > 0) {
        fprintf(fp, "par ");
        for (i = 0; i < NUPAR; i++) {
            if (j == 8) {
                fprintf(fp, "\npar ");
                j = 0;
            }
            get_val(upar_names[i], &z);
            fprintf(fp, "%s=%g ", upar_names[i], z);
            j++;
        }
    }
    fprintf(fp, "\n");
    fprintf(fp, "done \n");
    fclose(fp);
    return;
}

int32
init_conds_find_user_name(int32 type, char *oname) {
    char name[25];
    int32 k = 0, i = -1;
    for (usize j = 0; j < strlen(oname); j++) {
        if (!isspace(oname[j])) {
            name[k] = oname[j];
            k++;
        }
    }
    name[k] = 0;

    for (i = 0; i < NUPAR; i++)
        if ((type == PARAMBOX) && (strcasecmp(upar_names[i], name) == 0))
            break;
    if (i < NUPAR)
        return i;
    for (i = 0; i < NEQ; i++)
        if ((type == ICBOX) && (strcasecmp(uvar_names[i], name) == 0))
            break;
    if (i < NEQ)
        return i;
    return -1;
}

void
init_conds_create_par_sliders(Window base, int32 x0, int32 h0) {
    int32 i;
    for (i = 0; i < 3; i++)
        make_par_slider(base, x0 + i*36*DCURXs, h0, 100, i);
    return;
}

void
init_conds_resize_par_slides(int32 h) {
    for (int32 i = 0; i < 3; i++) {
        XMoveResizeWindow(display, my_par_slide[i].main, 10 + 36*i*DCURXs,
                          h, 32*(uint)DCURXs, 3*(uint)(DCURYs + 2));
    }
    return;
}

void
init_conds_slide_button_press(Window window) {
    for (int32 i = 0; i < 3; i++)
        do_slide_button(window, &my_par_slide[i]);
    return;
}

void
do_slide_button(Window window, struct ParSlider *p) {
    static char *n[] = {"*3Par/Var", "Value", "Low", "High"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    int32 status;
    double lo, hi, val;
    if (window == p->go && p->use == 1)
        integrate_run_now();

    if (window != p->top)
        return;
    strcpy(values[0], p->parname);
    sprintf(values[1], "%.16g", p->val);
    sprintf(values[2], "%.16g", p->lo);
    sprintf(values[3], "%.16g", p->hi);
    status = do_string_box(4, 4, 1, "Set Sliders", n, values, 35);
    if (status == 0)
        return;
    if (strlen(values[0]) == 0) { /* empty string cancels */
        p->use = 0;
        return;
    }
    status = init_conds_find_user_name(PARAMBOX, values[0]);
    if (status == -1) {
        status = init_conds_find_user_name(ICBOX, values[0]);
        if (status == -1) {
            ggets_err_msg("Not a parameter or variable !");
            return;
        }
        p->type = ICBOX;
        p->index = status;
    } else {
        p->type = PARAMBOX;
        p->index = status;
    }
    lo = atof(values[2]);
    hi = atof(values[3]);
    val = atof(values[1]);
    if (val < lo || val > hi || hi <= lo) {
        ggets_err_msg(" low <= value <= high ");
        return;
    }
    p->val = val;
    p->hi = hi;
    p->lo = lo;
    strcpy(p->parname, values[0]);
    set_val(p->parname, val);
    if (p->type == ICBOX)
        last_ic[p->index] = val;
    init_conds_redraw_params();
    init_conds_redraw_ics();
    p->use = 1;
    set_slide_pos(p);
    redraw_slide(p);
    return;
}

void
expose_selector(Window window) {
    display_file_sel(filesel, window);
    return;
}

/* this is rather lazy and slow but hey it works */

void
redraw_directory(void) {
    XClearWindow(display, filesel.dir);
    expose_selector(filesel.dir);
    return;
}

void
redraw_file_list(void) {
    for (int32 i = 0; i < filesel.nwin; i++) {
        XClearWindow(display, filesel.window[i]);
        expose_selector(filesel.window[i]);
    }
    return;
}

void
redraw_fs_text(char *string, Window window, int32 flag) {
    XClearWindow(display, window);
    filesel.off = 0;
    if (flag)
        filesel.pos = (int32)strlen(string);
    XDrawString(display, window, small_gc, 0, CURY_OFF, string,
                (int)strlen(string));
    if (flag)
        put_edit_cursor(window, DCURXs*(int32)strlen(string));
    return;
}

void
display_file_sel(struct FileSel f, Window window) {
    int32 i0;
    Window root;
    int32 xloc;
    int32 yloc;

    uint32 cwid;
    uint32 chgt;
    uint32 cbwid;
    uint32 cdepth;

    char t[sizeof(f.title) + sizeof(cur_dir) + sizeof(f.wildtxt) + 1];
    int32 hgt = DCURYs + 4;

    XGetGeometry(display, f.base, &root, &xloc, &yloc, &cwid, &chgt, &cbwid,
                 &cdepth);
    XResizeWindow(display, f.wild, cwid - 7*(uint)DCURXs - 5, (uint)DCURYs);
    XResizeWindow(display, f.file, cwid - 7*(uint)DCURXs - 5, (uint)DCURYs);
    for (int32 i = 0; i < f.nwin; i++) {
        XResizeWindow(display, f.window[i], cwid - 6*(uint)DCURXs - 10,
                      (uint)DCURYs);
    }
    XMoveResizeWindow(display, f.ok, (int)cwid / 2 - 7*DCURXs - 3,
                      (int32)chgt - hgt, 7*(uint)DCURXs, (uint)DCURYs);
    XMoveResizeWindow(display, f.cancel, cwid / 2 + 3, (int)chgt - hgt,
                      7*(uint)DCURXs, (uint)DCURYs);

    if (f.here != 1)
        return;
    if (f.ok == window)
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Ok", 2);
    if (f.cancel == window)
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Cancel", 6);
    if (f.up == window) {
        /*XDrawString(display,w,small_gc,5+DCURX/2,CURY_OFFs,"^",1);
         */
    }
    if (f.dn == window) {
        /*XDrawString(display,w,small_gc,5+DCURX/2,CURY_OFFs,"vv",1);
         */
    }
    if (f.pgup == window) {
        /*XDrawString(display,w,small_gc,5,CURY_OFFs,"^^",2);
         */
    }
    if (f.pgdn == window) {
        /* XDrawString(display,w,small_gc,5,CURY_OFFs,"vv",2);
         */
    }
    if (f.file == window) {
        XClearWindow(display, window);
        XDrawString(display, window, small_gc, 2, CURY_OFFs, f.filetxt,
                    (int)strlen(f.filetxt));
    }
    if (f.wild == window) {
        XClearWindow(display, window);
        XDrawString(display, window, small_gc, 2, CURY_OFFs, f.wildtxt,
                    (int)strlen(f.wildtxt));
    }
    if (f.fw == window)
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "File: ", 6);
    if (f.ww == window)
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Wild: ", 6);
    if (f.dir == window) {
        XTextProperty windowName;
        snprintf(t, sizeof(t), " %s", f.title);
        XDrawString(display, window, small_gc, 0, CURY_OFFs, t, (int)strlen(t));
        snprintf(t, sizeof(t), "%s - %s", f.wildtxt, cur_dir);
        {
            char *nameit[] = {t};
            XStringListToTextProperty(nameit, 1, &windowName);
            XSetWMName(display, f.base, &windowName);
        }
    }
    for (int32 i = 0; i < f.nwin; i++) {
        if (window == f.window[i]) {
            i0 = i + f.n0;
            if (i0 >= f.n)
                XDrawString(display, window, small_gc, 5, CURY_OFFs, " ", 1);
            else {
                if (i0 < my_ff.ndirs)
                    sprintf(t, "<>%s", my_ff.dirnames[i0]);
                else
                    sprintf(t, "%s", my_ff.filenames[i0 - my_ff.ndirs]);
                XDrawString(display, window, small_gc, 5, CURY_OFFs, t,
                            (int)strlen(t));
            }
        }
    }
    return;
}

void
fs_scroll(int32 i) {
    int32 n0 = filesel.n0;
    int32 new;
    int32 nend;
    int32 nw = filesel.nwin, n = filesel.n;
    if (n <= nw)
        return;
    new = n0 - i;
    nend = new + nw;
    if (new < 0)
        new = 0;
    if (nend > n)
        new = n - nw;
    filesel.n0 = new;
    redraw_file_list();
}

int32
button_selector(Window window) {
    int32 i0;
    int32 k, n = filesel.n;
    if (window == filesel.ok)
        return 1;
    if (window == filesel.cancel)
        return 2;
    if (window == filesel.up)
        fs_scroll(1);
    if (window == filesel.dn)
        fs_scroll(-1);
    if (window == filesel.pgup)
        fs_scroll(filesel.nwin);
    if (window == filesel.pgdn)
        fs_scroll(-filesel.nwin);
    if (window == filesel.home) {
        char *HOMEDIR = getenv("KEY_HOME");
        int32 m;
        if ((HOMEDIR == NULL) || (strlen(HOMEDIR) == 0)) {
            ggets_plintf("User's KEY_HOME environment variable not set.\n");
            return 0;
        }
        change_directory(HOMEDIR);

        get_directory(cur_dir);
        redraw_directory();
        free_finfo(&my_ff); /* delete the old file info */
        filesel.n0 = 0;     /* back to the top of the list */
        get_fileinfo(filesel.wildtxt, cur_dir, &my_ff);
        filesel.n = my_ff.ndirs + my_ff.nfiles;

        strcpy(filesel.filetxt, cur_dir);

        m = (int32)strlen(filesel.filetxt);
        if (filesel.filetxt[m - 1] != '/') {
            strcat(filesel.filetxt, "/");
        }

        redraw_fs_text(filesel.filetxt, filesel.file, 1);
        redraw_file_list();
        XFlush(display);

        return 0;
    }
    if (window == filesel.start) {
        char *START = getenv("XPPSTART");
        int32 m;

        if ((START == NULL) || (strlen(START) == 0)) {
            ggets_plintf("User's XPPSTART environment variable not set.\n");
            return 0;
        }

        change_directory(START);

        get_directory(cur_dir);
        redraw_directory();
        free_finfo(&my_ff); /* delete the old file info */
        filesel.n0 = 0;     /* back to the top of the list */
        get_fileinfo(filesel.wildtxt, cur_dir, &my_ff);
        filesel.n = my_ff.ndirs + my_ff.nfiles;

        strcpy(filesel.filetxt, cur_dir);

        m = (int32)strlen(filesel.filetxt);
        if (filesel.filetxt[m - 1] != '/') {
            strcat(filesel.filetxt, "/");
        }

        redraw_fs_text(filesel.filetxt, filesel.file, 1);
        redraw_file_list();
        XFlush(display);

        return 0;
    }
    if (window == filesel.file) { /* selected the file text */
        if (filesel.hot != HOTFILE)
            filesel.pos = (int32)strlen(filesel.filetxt);

        filesel.hot = HOTFILE;
        redraw_fs_text(filesel.filetxt, filesel.file, 1);
        redraw_fs_text(filesel.wildtxt, filesel.wild, 0);
        /* set up text stuff */
        return 0;
    }
    if (window == filesel.wild) {
        if (filesel.hot != HOTWILD)
            filesel.pos = (int32)strlen(filesel.wildtxt);
        filesel.hot = HOTWILD;
        redraw_fs_text(filesel.filetxt, filesel.file, 0);
        redraw_fs_text(filesel.wildtxt, filesel.wild, 1);
        return 0;
    }
    i0 = -1;
    for (int32 i = 0; i < filesel.nwin; i++) {
        if (window == filesel.window[i])
            i0 = i;
    }
    if (i0 > -1) { /* clicked on a file or directory */
        k = i0 + filesel.n0;
        if (k < my_ff.ndirs) { /* it is a directory so we should reset */
            int32 m;
            change_directory(my_ff.dirnames[k]);
            get_directory(cur_dir);
            redraw_directory();
            free_finfo(&my_ff); /* delete the old file info */
            filesel.n0 = 0;     /* back to the top of the list */
            get_fileinfo(filesel.wildtxt, cur_dir, &my_ff);
            filesel.n = my_ff.ndirs + my_ff.nfiles;

            strcpy(filesel.filetxt, cur_dir);

            m = (int32)strlen(filesel.filetxt);
            if (filesel.filetxt[m - 1] != '/') {
                strcat(filesel.filetxt, "/");
            }

            redraw_fs_text(filesel.filetxt, filesel.file, 1);
            redraw_file_list();
            XFlush(display);

            return 0;
        }
        if (k < n) {
            k = k - my_ff.ndirs;
            strcpy(filesel.filetxt, my_ff.filenames[k]);
            return 1; /* got a file */
        }
    }
    return 0;
}

void
crossing_selector(Window window, int32 c) {
    int32 t1 = 1, t2 = 2;
    Window w = window;
    if (c == 1) {
        t1 = 0;
        t2 = 1;
    }
    for (int32 i = 0; i < filesel.nwin; i++) {
        if (w == filesel.window[i]) {
            XSetWindowBorderWidth(display, w, (uint)t1);
            return;
        }
    }
    if (w == filesel.ok || w == filesel.cancel || w == filesel.pgup ||
        w == filesel.pgdn || w == filesel.up || w == filesel.dn ||
        w == filesel.file || w == filesel.wild || w == filesel.home ||
        w == filesel.start)
        XSetWindowBorderWidth(display, w, (uint)t2);
    return;
}

void
create_file_selector(char *title, char *file, char *wild) {
    int32 n = my_ff.ndirs + my_ff.nfiles;
    int32 nwin = FILESELNWIN;
    int32 hgt;
    int32 i;
    int32 width, height;

    Window base;
    XTextProperty winname;
    XSizeHints size_hints;
    filesel.n = n;
    filesel.n0 = 0;
    filesel.nwin = nwin;
    strcpy(filesel.title, title);
    strcpy(filesel.wildtxt, wild);
    strcpy(filesel.filetxt, file);
    width = 80*DCURXs;
    /*wid=30*DCURXs;*/
    hgt = DCURYs + 4;
    height = (5 + nwin)*hgt;
    filesel.minwid = width;
    filesel.minhgt = height;
    /* printf("Here now 23!\n");*/
    base =
        make_plain_window(RootWindow(display, screen), 0, 0, width, height, 4);
    /* printf("Here now 23!\n"); */
    filesel.base = base;
    XStringListToTextProperty(&title, 1, &winname);
    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.x = 0;
    size_hints.y = 0;
    size_hints.width = width;
    size_hints.height = height;
    size_hints.min_width = width;
    size_hints.min_height = height;
    size_hints.max_width = width;
    size_hints.max_height = height;

    many_pops_make_icon((char *)filebrowse_bits, filebrowse_width,
                        filebrowse_height, base);

    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";

        XSetWMProperties(display, base, &winname, NULL, NULL, 0, &size_hints,
                         NULL, &class_hints);
    }

    /*
    filesel.up=make_window(base,DCURXs,2+4*hgt,3*DCURXs+5,DCURYs,1);
    filesel.dn=make_window(base,DCURXs,2+5*hgt,3*DCURXs+5,DCURYs,1);
    filesel.pgup=make_plain_window(base,DCURXs,2+7*hgt,3*DCURXs+5,DCURYs,1);
    filesel.pgdn=make_window(base,DCURXs,2+8*hgt,3*DCURXs+5,DCURYs,1);
    */

    filesel.up = make_icon_window(base, DCURXs, 2 + 3*hgt + 72 + 15, 32, 24,
                                  1, lineup_bits);
    filesel.dn = make_icon_window(base, DCURXs, 2 + 3*hgt + 96 + 18, 32, 24,
                                  1, linedn_bits);
    filesel.pgup = make_icon_window(base, DCURXs, 2 + 3*hgt + 48 + 12, 32, 24,
                                    1, pageup_bits);
    filesel.pgdn = make_icon_window(base, DCURXs, 2 + 3*hgt + 120 + 21, 32,
                                    24, 1, pagedn_bits);
    filesel.home =
        make_icon_window(base, DCURXs, 2 + 3*hgt, 32, 24, 1, home_bits);
    filesel.start = make_icon_window(base, DCURXs, 2 + 3*hgt + 24 + 3, 32, 24,
                                     1, start_bits);

    filesel.dir = make_plain_window(base, 7*DCURXs, 2, width - 7*DCURXs - 5,
                                    DCURYs, 0);
    filesel.wild = make_plain_window(base, 7*DCURXs, 2 + hgt,
                                     width - 7*DCURXs - 5, DCURYs, 1);
    filesel.ww = make_window(base, 2, 2 + hgt, 6*DCURXs + 2, DCURYs, 0);
    filesel.file = make_plain_window(base, 7*DCURXs, 2 + 2*hgt,
                                     width - 7*DCURXs - 5, DCURYs, 1);
    filesel.fw = make_window(base, 2, 2 + 2*hgt, 6*DCURXs + 2, DCURYs, 0);
    for (i = 0; i < nwin; i++) {
        filesel.window[i] =
            make_plain_window(base, 6*DCURXs + 5, 2 + (3 + i)*hgt,
                              width - 6*DCURXs - 10, DCURYs, 0);
    }

    filesel.ok = make_window(base, width / 2 - 7*DCURXs - 3, height - hgt,
                             7*DCURXs, DCURYs, 1);
    filesel.cancel =
        make_window(base, width / 2 + 3, height - hgt, 7*DCURXs, DCURYs, 1);
    /*  XSelectInput(display,filesel.wild,BOXEVENT);
        XSelectInput(display,filesel.file,BOXEVENT); */
    filesel.here = 1;
    filesel.hot = HOTFILE;
    filesel.pos = (int32)strlen(filesel.filetxt);
    filesel.off = 0;
    return;
}

void
stringintersect(char *target, char *sother) {
    int32 m = (int32)strlen(target);
    int32 n = (int32)strlen(sother);
    int32 j = 0;
    if (n < m)
        m = n;
    while (j < m) {
        if (target[j] != sother[j]) {
            break;
        }
        j++;
    }
    target[j] = '\0';
    return;
}

static int32 edit_fitem(int32 ch, char *string, Window window, int32 *off1,
                        int32 *pos1, int32 mc);
int32
edit_fitem(int32 ch, char *string, Window window, int32 *off1, int32 *pos1,
           int32 mc) {
    int32 l = (int32)strlen(string), cp;
    int32 off = *off1, pos = *pos1, wpos = pos - off;
    switch (ch) {
    case KEY_LEFT:
        if (pos > 0) {
            pos--;
            wpos--;
            if (wpos < 0) {
                off = off - 4;
                if (off < 0)
                    off = 0;
                wpos = pos - off;
            }
        } else
            ggets_ping();
        break;
    case KEY_RIGHT:
        if (pos < l) {
            pos++;
            wpos++;
            if (wpos > mc) {
                off = off + 4;
                if (off + mc > l)
                    off = l - mc;
                wpos = pos - off;
            }
        } else
            ggets_ping();
        break;
    case KEY_HOME:
        pos = 0;
        wpos = 0;
        break;
    case KEY_END:
        pos = l;
        wpos = mc;
        break;
    case KEY_BADKEY:
        return 0;

    case KEY_DOWN:
        fs_scroll(-1);
        return 0;
    case KEY_UP:
        fs_scroll(1);
        return 0;
    case KEY_PGUP:
        fs_scroll(filesel.nwin);
        return 0;
    case KEY_PGDN:
        fs_scroll(-filesel.nwin);
        return 0; /* junk key  */
    case KEY_ESC:
        return EDIT_ESC;
    case KEY_FINE:
        return EDIT_DONE;
    case KEY_BKSP:
        /*
        if(pos<l){
          ggets_mem_mov(&string[pos],&string[pos+1],l-pos);
          l--;
        }
        else
         ggets_ping();
         break; */
    case KEY_DEL:

        if (pos > 0) {
            ggets_mem_mov(&string[pos - 1], &string[pos], l - pos + 1);
            pos--;
            wpos--;
            if (wpos < 0) {
                off = off - 4;
                if (off < 0)
                    off = 0;
                wpos = pos - off;
            }
            l--;
        } else
            ggets_ping();
        break;
    case KEY_TAB: /*KEY_TAB completion of file names */
    {
        struct dirent *dp;
        /*char ft[100];
        char ftpath[100];
        */

        char ft[XPP_MAX_NAME];
        char ftpath[XPP_MAX_NAME];
        int32 m;

        /*User may have typed ahead (maybe they remember the path they want)"*/
        /*Try to change to that new directory if it is one.*/
        if ((dp = (struct dirent *)opendir(filesel.filetxt)) != NULL) {
            if (strcmp(cur_dir, filesel.filetxt) != 0) {
                change_directory(filesel.filetxt);
                get_directory(cur_dir);
                redraw_directory();
                free_finfo(&my_ff); /* delete the old file info */
                filesel.n0 = 0;     /* back to the top of the list */
                get_fileinfo(filesel.wildtxt, cur_dir, &my_ff);
                filesel.n = my_ff.ndirs + my_ff.nfiles;
                strcpy(filesel.filetxt, cur_dir);
                m = (int32)strlen(filesel.filetxt);
                if (filesel.filetxt[m - 1] != '/') {
                    strcat(filesel.filetxt, "/");
                }
                redraw_fs_text(filesel.filetxt, filesel.file, 1);
                redraw_file_list();
                XFlush(display);
            }
            return EDIT_WAIT; /*Wait for further instruction...*/
        }

        m = (int32)strlen(filesel.filetxt) + 1;
        // strcpy(ft,filesel.filetxt);

        if (m > 1) {
            strcpy(ft, filesel.filetxt);
        } else {
            /* We are already at the root level of file system! */
            strcpy(filesel.filetxt, "/");
            strcpy(ft, filesel.filetxt);
        }

        while ((dp = (struct dirent *)opendir(ft)) == NULL) {
            /* This WHILE is perhaps a bit clunky but since we can't be sure of
             * path separator user will type in the box a trial-by-error
             * approach may be more robust. */

            ft[m] = '\0';
            m--;
            if ((ft[m] != '/') & (ft[m] != '\\')) {
                ft[m] = '\0';
                m--;
            }

            if (m == 0) {
                break;
            }
        }

        ft[0] = '\0';
        if (m > (int32)strlen(filesel.filetxt)) {
            return EDIT_WAIT;
        }
        for (usize n = 0; n < strlen(filesel.filetxt) - (usize)m; n++) {
            ft[n] = filesel.filetxt[(usize)m + n + 1];
            ft[n + 1] = '\0';
        }
        strcat(ft, "*");
        strcpy(ftpath, filesel.filetxt);
        ftpath[m + 1] = '\0';
        /*Make sure we are in the correct directory now
        since user could have moved cursor back up _several_
        branches in directory tree before hitting tab key.*/
        change_directory(ftpath);
        free_finfo(&my_ff);
        filesel.n0 = 0;
        get_fileinfo_tab(ft, ftpath, &my_ff, filesel.wildtxt);
        filesel.n = my_ff.ndirs + my_ff.nfiles;
        if ((my_ff.ndirs + my_ff.nfiles) == 1) {
            if (my_ff.ndirs == 1) /*Only possible directory -- take it.*/
            {
                int32 m2;
                change_directory(my_ff.dirnames[0]);
                get_directory(cur_dir);
                redraw_directory();
                free_finfo(&my_ff); /* delete the old file info */
                filesel.n0 = 0;     /* back to the top of the list */
                get_fileinfo(filesel.wildtxt, cur_dir, &my_ff);
                filesel.n = my_ff.ndirs + my_ff.nfiles;
                strcpy(filesel.filetxt, cur_dir);
                m2 = (int32)strlen(filesel.filetxt);
                if (filesel.filetxt[m2 - 1] != '/') {
                    strcat(filesel.filetxt, "/");
                }
                redraw_fs_text(filesel.filetxt, filesel.file, 1);
            } else /*Must be that (my_ff.nfiles == 1) -- but best wait for user
                      to confim they actually want this file by clicking on
                      it.*/
            {
                /*Copy the only remaining choice into the file box and wait for
                user to make final choice.
                */
                strcpy(filesel.filetxt, cur_dir);
                strcat(filesel.filetxt, "/");
                strcat(filesel.filetxt, my_ff.filenames[0]);
                redraw_fs_text(filesel.filetxt, filesel.file, 1);
            }
        } else if (filesel.n > 1) /*Expand the file text to most generic
                                     currently represented*/
        {
            int32 j = 0;
            char U[256];
            if (my_ff.ndirs > 0) {
                strcpy(U, my_ff.dirnames[0]);
            } else {
                strcpy(U, my_ff.filenames[0]);
            }

            for (j = 0; j < filesel.n; j++) {
                if (j < my_ff.ndirs) {
                    stringintersect(U, my_ff.dirnames[j]);
                } else {
                    stringintersect(U, my_ff.filenames[j - my_ff.ndirs]);
                }

                if (strlen(U) == 0) /*No common substring*/
                {
                    break;
                }
            }
            strcpy(filesel.filetxt, ftpath);
            strcat(filesel.filetxt, U);
            /*Actually don't want to force appending of path separator here
            since we want user to decide between matching directory and files.
            */
            redraw_fs_text(filesel.filetxt, filesel.file, 1);
        }
        redraw_file_list();
        XFlush(display);
        return 0;
    }
    default:
        if ((ch >= ' ') && (ch <= '~')) {
            if (strlen(string) >= 256)
                ggets_ping();
            else {
                ggets_mov_mem(&string[pos + 1], &string[pos], l - pos + 1);
                string[pos] = (char)ch;
                pos = pos + 1;
                wpos++;
                l++;
                if (wpos > mc) {
                    off = off + 4;
                    if (off + mc > l)
                        off = l - mc;
                    wpos = pos - off;
                }
            }
        }
        break;
    }
    /* all done lets save everything */
    off = pos - wpos;
    *off1 = off;
    *pos1 = pos;
    XClearWindow(display, window);
    XDrawString(display, window, small_gc, 0, CURY_OFF, string + off,
                (int32)strlen(string) - off);
    cp = DCURXs*(pos - off);
    put_edit_cursor(window, cp);
    return 0;
}

int32
selector_key(XEvent event) {
    char ch;
    int32 flag;
    ch = (char)ggets_get_key_press(&event);
    switch (filesel.hot) {
    case HOTFILE:
        flag = edit_fitem(ch, filesel.filetxt, filesel.file, &(filesel.off),
                          &(filesel.pos), 29);
        if (flag == EDIT_DONE)
            return 1;
        if (flag == EDIT_ESC)
            return 2;
        return 0;
    case HOTWILD:
        flag = edit_fitem(ch, filesel.wildtxt, filesel.wild, &(filesel.off),
                          &(filesel.pos), 29);
        if (flag == EDIT_DONE) {
            /* new wild */
            free_finfo(&my_ff); /* delete the old file info */
            filesel.n0 = 0;     /* back to the top of the list */
            get_fileinfo(filesel.wildtxt, cur_dir, &my_ff);
            filesel.n = my_ff.ndirs + my_ff.nfiles;
            redraw_file_list();
            XFlush(display);
            return 0;
        }
        if (flag == EDIT_ESC)
            return 2;
        return 0;
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
}

int32
init_conds_file_selector(char *title, char *file, char *wild) {
    int32 selected;
    int32 done;
    if (!get_directory(cur_dir))
        return 0;
    if (!get_fileinfo(wild, cur_dir, &my_ff))
        return 0;

    create_file_selector(title, file, wild);

    /* do file select events */
    while (true) {
        XEvent event;
        XNextEvent(display, &event);
        switch (event.type) {
        case ConfigureNotify:
        case Expose:
        case MapNotify:
            if (Xup)
                many_pops_do_expose(event);
            expose_selector(event.xany.window);
            break;
        case ButtonPress:
            done = button_selector(event.xbutton.window);
            if (done == 1) {
                selected = 1; /* OK made a selection */
                goto end;
            }
            if (done == 2) {
                selected = 0; /* canceled the whole thing */
                goto end;
            }
            break;
        case EnterNotify:
            crossing_selector(event.xcrossing.window, 0);
            break;
        case LeaveNotify:
            crossing_selector(event.xcrossing.window, 1);
            break;
        case KeyPress:
            done = selector_key(event);
            if (done == 2) {
                selected = 0;
                goto end;
            }
            if (done == 1) {
                selected = 1;
                goto end;
            }
            break;
        default:
            break;
        }
    }
end:

    /* destroy selector */
    filesel.here = 0;
    browse_wait_a_sec(ClickTime);
    XDestroySubwindows(display, filesel.base);
    XDestroyWindow(display, filesel.base);
    free_finfo(&my_ff);

    XFlush(display); /*Need to do this otherwise the file dialog hangs around*/
    if (selected == 0)
        return 0;
    strcpy(file, filesel.filetxt);
    return 1; /* got a file name */
}

void
init_conds_reset_sliders(void) {
    int32 i;
    double val;
    struct ParSlider *p;
    for (i = 0; i < 3; i++) {
        p = &my_par_slide[i];
        if (p->use) {
            if (p->type == ICBOX)
                val = last_ic[p->index];
            else
                get_val(p->parname, &val);
            p->val = val;
            set_slide_pos(p);
            expose_slider(p->slide, p);
            expose_slider(p->top, p);
        }
    }
    return;
}

void
redraw_slide(struct ParSlider *p) {
    expose_slider(p->slide, p);
    expose_slider(p->top, p);
    expose_slider(p->left, p);
    expose_slider(p->right, p);
    return;
}

void
set_slide_pos(struct ParSlider *p) {
    double pos;
    int32 ip;
    pos = 2. + (p->l - 4)*(p->val - p->lo) / (p->hi - p->lo);
    ip = (int32)pos;
    if (ip < 2)
        ip = 2;
    if (ip > (p->l - 2))
        ip = p->l - 2;
    p->pos = ip;
    return;
}

void
init_conds_slide_release(Window window) {
    for (int32 i = 0; i < 3; i++)
        do_slide_release(window, &my_par_slide[i]);
    return;
}

void
do_slide_release(Window window, struct ParSlider *p) {
    if (p->use == 0)
        return;
    if (p->slide == window) {
        set_val(p->parname, p->val);
        if (p->type == ICBOX)
            last_ic[p->index] = p->val;
        init_conds_redraw_ics();
        init_conds_redraw_params();
    }
    return;
}

void
init_conds_slider_motion(XEvent event) {
    int32 x;
    Window window;
    window = event.xmotion.window;
    x = event.xmotion.x;
    for (int32 i = 0; i < 3; i++)
        do_slide_motion(window, x, &my_par_slide[i],
                        (int32)event.xmotion.state);
    return;
}

void
do_slide_motion(Window window, int32 x, struct ParSlider *p, int32 s) {
    int32 sp = SuppressBounds;
    if (window == p->slide) {
        p->pos = x;
        if (x < 2)
            p->pos = 2;
        if (x > (p->l - 2))
            p->pos = p->l - 2;
        expose_slider(p->slide, p);
        if (p->use) {
            p->val = p->lo + (p->hi - p->lo)*(double)(p->pos - 2) /
                                 (double)(p->l - 4);
            expose_slider(p->top, p);
            set_val(p->parname, p->val);
            if (p->type == ICBOX)
                last_ic[p->index] = p->val;
            if (s < 300) {
                menudrive_clr_all_scrns();
                nullcline_redraw_dfield();
                nullcline_create_new_cline();
                many_pops_draw_label(draw_win);
                SuppressBounds = 1;
                integrate_run_now();
                SuppressBounds = sp;
            }
        }
    }
    return;
}

void
init_conds_enter_slides(Window window, int32 val) {
    int32 i;
    for (i = 0; i < 3; i++)
        enter_slider(window, &my_par_slide[i], val);
    return;
}

void
enter_slider(Window window, struct ParSlider *p, int32 val) {
    if (window == p->top || window == p->go)
        XSetWindowBorderWidth(display, window, (uint)val + 1);
    return;
}

void
init_conds_expose_slides(Window window) {
    int32 i;
    for (i = 0; i < 3; i++)
        expose_slider(window, &my_par_slide[i]);
    return;
}

void
expose_slider(Window window, struct ParSlider *p) {
    int32 x, len = 12*DCURXs;
    char top[256];
    if (window == p->slide) {
        draw_slider(window, p->pos, p->hgt, p->l);
        return;
    }
    if (window == p->go) {
        XDrawString(display, window, small_gc, 2,
                    (int32)(0.75*(double)CURY_OFFs), "go", 2);
        return;
    }
    if (p->use) {
        if (window == p->left) {
            sprintf(top, "%.16g", p->lo);
            x = 1;
            XClearWindow(display, window);
            XDrawString(display, window, small_gc, x, CURY_OFFs, top,
                        (int)strlen(top));
            return;
        }
        if (window == p->right) {
            sprintf(top, "%.16g", p->hi);
            x = 1;
            if (strlen(top) < 12)
                x = len - DCURXs*(int32)strlen(top) - 1;
            XClearWindow(display, window);
            XDrawString(display, window, small_gc, x, CURY_OFFs, top,
                        (int)strlen(top));
            return;
        }
        if (window == p->top) {
            sprintf(top, "%s=%.16g", p->parname, p->val);
            XClearWindow(display, window);
            XDrawString(display, window, small_gc, 2, CURY_OFFs, top,
                        (int)strlen(top));
        }
    } else {
        if (window == p->top) {
            sprintf(top, "Par/Var?");
            x = 1;
            XClearWindow(display, window);
            XDrawString(display, window, small_gc, x, CURY_OFFs, top,
                        (int)strlen(top));
        }
    }
    return;
}

void
draw_slider(Window window, int32 x, int32 hgt, int32 l) {
    int32 x0 = x - 2, i;
    if (x0 < 0)
        x0 = 0;
    if (x0 > (l - 4))
        x0 = l - 4;
    XClearWindow(display, window);
    for (i = 0; i < 4; i++)
        XDrawLine(display, window, small_gc, x0 + i, 0, x0 + i, hgt);
    return;
}

void
make_par_slider(Window base, int32 x, int32 y, int32 width, int32 index) {
    int32 mainhgt = 3*(DCURYs + 2);
    int32 mainwid = 32*DCURXs;
    int32 xs;
    Window window;
    if (mainwid < (width + 4))
        mainwid = width + 4;

    window = make_plain_window(base, x, y, mainwid, mainhgt, 1);
    my_par_slide[index].main = window;
    xs = (mainwid - width - 4) / 2;
    my_par_slide[index].slide =
        make_window(window, xs, DCURYs + 5, width + 4, DCURYs - 4, 1);
    my_par_slide[index].go = make_window(window, xs + width + 8, DCURYs + 5,
                                         3*DCURXs, DCURYs - 3, 1);
    my_par_slide[index].top = make_window(window, 2, 2, mainwid - 6, DCURYs, 1);
    my_par_slide[index].left =
        make_window(window, 2, 2*DCURYs + 3, 12*DCURXs, DCURYs, 0);
    my_par_slide[index].right =
        make_window(window, mainwid - 12*DCURXs - 4, 2*DCURYs + 3,
                    12*DCURXs, DCURYs, 0);
    my_par_slide[index].lo = 0.0;
    my_par_slide[index].hi = 1.0;
    my_par_slide[index].val = 0.5;
    my_par_slide[index].use = 0;
    my_par_slide[index].l = width + 4;
    my_par_slide[index].pos = (width + 4) / 2;
    my_par_slide[index].parname[0] = 0;
    my_par_slide[index].hgt = DCURYs - 4;

    if ((notAlreadySet.SLIDER1 == 0) && (index == 0)) {
        strcpy(my_par_slide[index].parname, SLIDER1VAR);
        my_par_slide[index].use = 1;
        my_par_slide[index].lo = SLIDER1LO;
        my_par_slide[index].hi = SLIDER1HI;
        get_val(my_par_slide[index].parname, &my_par_slide[index].val);
    }

    if ((notAlreadySet.SLIDER2 == 0) && (index == 1)) {
        strcpy(my_par_slide[index].parname, SLIDER2VAR);
        my_par_slide[index].use = 1;
        my_par_slide[index].lo = SLIDER2LO;
        my_par_slide[index].hi = SLIDER2HI;
        get_val(my_par_slide[index].parname, &my_par_slide[index].val);
    }

    if ((notAlreadySet.SLIDER3 == 0) && (index == 2)) {
        strcpy(my_par_slide[index].parname, SLIDER3VAR);
        my_par_slide[index].use = 1;
        my_par_slide[index].lo = SLIDER3LO;
        my_par_slide[index].hi = SLIDER3HI;
        get_val(my_par_slide[index].parname, &my_par_slide[index].val);
    }
    return;
}

/*     The rest of the code is good
                    |
                    V
  */

void
init_conds_make_new_ic_box(void) {
    if (ICBox.xuse) {
        XRaiseWindow(display, ICBox.base);
        return;
    }
    make_box_list_window(&ICBox, ICBOX);
    many_pops_make_icon((char *)ic_bits, ic_width, ic_height, ICBox.base);
    return;
}

void
init_conds_make_new_bc_box(void) {
    if (BCBox.xuse) {
        XRaiseWindow(display, BCBox.base);
        return;
    }
    make_box_list_window(&BCBox, BCBOX);
    many_pops_make_icon((char *)bc_bits, bc_width, bc_height, BCBox.base);
    return;
}

void
init_conds_make_new_delay_box(void) {
    if (DelayBox.use == 0)
        return;
    if (DelayBox.xuse == 1) {
        XRaiseWindow(display, DelayBox.base);
        return;
    }
    make_box_list_window(&DelayBox, DELAYBOX);
    many_pops_make_icon((char *)delay_bits, delay_width, delay_height,
                        DelayBox.base);
    return;
}

void
init_conds_make_new_param_box(void) {
    if (ParamBox.use == 0)
        return;
    if (ParamBox.xuse == 1) {
        XRaiseWindow(display, ParamBox.base);
        return;
    }
    make_box_list_window(&ParamBox, PARAMBOX);
    many_pops_make_icon((char *)param_bits, param_width, param_height,
                        ParamBox.base);
    return;
}

void
init_conds_initialize_box(void) {
    make_box_list(&ICBox, "Initial Data", "ICs", NODE + NMarkov, ICBOX, 1);
    if (NUPAR > 0)
        make_box_list(&ParamBox, "Parameters", "Par", NUPAR, PARAMBOX, 1);
    else
        ParamBox.use = 0;
    if (NDELAYS > 0)
        make_box_list(&DelayBox, "Delay ICs", "Delay", NODE, DELAYBOX, 1);
    else
        DelayBox.use = 0;
    make_box_list(&BCBox, "Boundary Conds", "BCs", NODE, BCBOX, 1);

    /*  Iconify them !!   */
    /*  if(noicon==0){
    if(ICBox.xuse)XIconifyWindow(display,ICBox.base,screen);
    if(DelayBox.xuse) XIconifyWindow(display,DelayBox.base,screen);
    if(ParamBox.xuse) XIconifyWindow(display,ParamBox.base,screen);
     if(BCBox.xuse)XIconifyWindow(display,BCBox.base,screen);
     } */
    return;
}

void
init_conds_resize_par_box(Window window) {
    uint32 h;
    uint32 w;
    int32 nwin = 0;
    int32 ok = 0;
    BoxList *b = NULL;

    if (ICBox.xuse == 1 && window == ICBox.base) {
        ok = 1;
        b = &ICBox;
        eig_list_get_new_size(window, &w, &h);
        get_nrow_from_hgt((int32)h, &nwin, (int32 *)&w);
    }

    if (ParamBox.xuse == 1 && window == ParamBox.base) {
        ok = 2;
        b = &ParamBox;
        browse_wait_a_sec(ClickTime);

        eig_list_get_new_size(window, &w, &h);
        get_nrow_from_hgt((int32)h, &nwin, (int32 *)&w);
    }
    if (BCBox.xuse == 1 && window == BCBox.base) {
        ok = 3;
        b = &BCBox;
        eig_list_get_new_size(window, &w, &h);
        get_nrow_from_hgt((int32)h, &nwin, (int32 *)&w);
    }
    if (DelayBox.xuse == 1 && window == DelayBox.base) {
        ok = 4;
        b = &DelayBox;
        eig_list_get_new_size(window, &w, &h);
        get_nrow_from_hgt((int32)h, &nwin, (int32 *)&w);
    }
    if (ok == 0)
        return;
    if (nwin > b->n)
        nwin = b->n;

    if (nwin == b->nwin)
        return;

    b->nwin = nwin;

    b->n0 = (b->n - b->nwin);
    /* b->n0=0;This is a work-around for the following bug (still happening):
                i) User scrolls to bottom of parameter list
                ii) then resizes window
                iii) then Segmentation Fault
   */

    switch (ok) {
    case 1:
        destroy_box(&ICBox);
        init_conds_make_new_ic_box();
        break;
    case 2:
        destroy_box(&ParamBox);
        init_conds_make_new_param_box();
        break;
    case 3:
        destroy_box(&BCBox);
        init_conds_make_new_bc_box();
        break;
    case 4:
        destroy_box(&DelayBox);
        init_conds_make_new_delay_box();
        break;
    default:
        break;
    }
    return;
}

/* this returns the fixed width, the number of entries
    allowed
*/

void
get_nrow_from_hgt(int32 h, int32 *n, int32 *w) {
    int32 hgt = DCURYs + 4;
    *w = 28*DCURXs;
    *n = h / (hgt + 4) - 5;
    return;
}

void
destroy_box(BoxList *b) {
    /*int32 n,nrow;
     */
    if (b->xuse == 0)
        return;
    b->xuse = 0;
    XFlush(display);
    XSetInputFocus(display, main_win, RevertToParent, CurrentTime);
    if (b->use == 0)
        return;
    browse_wait_a_sec(ClickTime);

    XDestroySubwindows(display, b->base);
    XDestroyWindow(display, b->base);

    /*n=b->n;
    nrow=b->nwin;
    */
    /* now free up stuff */
    free(b->w);
    free(b->we);
    if (b->type == ICBOX) {
        free(b->ck);
        free(b->isck);
    }
    browse_wait_a_sec(200);
    XFlush(display);
    return;
}

void
make_box_list_window(BoxList *b, int32 type) {
    int32 nrow;
    int32 n;
    int32 x, y;
    int32 xb1, xb2, xb3, xb4;
    int32 i, wid1, wid2;
    int32 width, height, wid, hgt;

    Window base;
    XTextProperty winname;
    XTextProperty iconame;
    XSizeHints size_hints;
    n = b->n;
    nrow = b->nwin;

    wid1 = 16*DCURXs;
    wid2 = 22*DCURXs;
    wid = wid1 + wid2 + DCURXs;
    hgt = DCURYs + 4;
    height = (nrow + 4)*(hgt + 4) + 2*hgt;
    width = wid + 8*DCURXs;
    b->minwid = width;
    b->minhgt = height;
    base =
        make_plain_window(RootWindow(display, screen), 0, 0, width, height, 4);
    b->base = base;
    XStringListToTextProperty(&b->wname, 1, &winname);
    XStringListToTextProperty(&b->iname, 1, &iconame);
    size_hints.flags = PPosition | PSize | PMinSize;
    size_hints.x = 0;
    size_hints.y = 0;
    size_hints.width = width;
    size_hints.height = height;
    size_hints.min_width = width;
    size_hints.min_height = height;
    size_hints.max_width = width;
    size_hints.max_height = height;
    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";
        XSetWMProperties(display, base, &winname, &iconame, NULL, 0,
                         &size_hints, NULL, &class_hints);
    }
    b->w = xmalloc((usize)nrow*sizeof(*(b->w)));
    b->we = xmalloc((usize)nrow*sizeof(*(b->we)));
    if (type == ICBOX) {
        b->ck = xmalloc((usize)nrow*sizeof(*(b->ck)));
        b->isck = xmalloc((usize)n*sizeof(*(b->isck)));
        for (i = 0; i < n; i++)
            b->isck[i] = 0;
    }

    xb1 = 2 + 7*DCURXs + 14;

    xb2 = xb1 + 7*DCURXs + 14;
    xb3 = xb2 + 7*DCURXs + 14;
    xb4 = xb3 + 7*DCURXs + 14;

    b->close = make_window(base, 2, 5, 7*DCURXs + 10, DCURYs, 1);
    b->ok = make_window(base, xb1, 5, 7*DCURXs + 10, DCURYs, 1);
    b->def = make_window(base, xb2, 5, 7*DCURXs + 10, DCURYs, 1);
    b->cancel = make_window(base, xb3, 5, 7*DCURXs + 10, DCURYs, 1);
    b->go = make_window(base, xb4, 5, 7*DCURXs + 10, DCURYs, 1);
    xb1 = DCURXs + wid1 + wid2 + 12;

    b->up = make_icon_window(base, xb1, (int32)(1.75*DCURYs) + 24 + 3, 32, 24,
                             1, lineup_bits);
    b->dn = make_icon_window(base, xb1, (int32)(1.75*DCURYs) + 48 + 6, 32, 24,
                             1, linedn_bits);
    b->pgup = make_icon_window(base, xb1, (int32)(1.75*DCURYs), 32, 24, 1,
                               pageup_bits);
    b->pgdn = make_icon_window(base, xb1, (int32)(1.75*DCURYs) + 72 + 9, 32,
                               24, 1, pagedn_bits);

    for (i = 0; i < nrow; i++) {
        x = DCURXs;
        y = DCURYs + (hgt + 4)*i + (int32)(1.5*hgt);
        b->w[i] = make_plain_window(base, x, y, wid1, hgt, 0);
        b->we[i] = make_plain_window(base, x + wid1 + 2, y, wid2, hgt, 1);
        XSelectInput(display, b->w[i], BOXEVENT);
        if (type == ICBOX) {
            b->ck[i] = make_plain_window(base, 1, y, 6, DCURYs, 1);
        }
    }

    y = DCURYs + (hgt + 4)*nrow + (int32)(1.5*hgt);
    x = (width - 24) / 3;
    if (type == ICBOX) {
        b->xvt = make_window(base, x, y, 5*DCURXs, DCURYs, 1);
        b->pp = make_window(base, x + 6*DCURXs, y, 5*DCURXs, DCURYs, 1);
        b->arr = make_window(base, x + 12*DCURXs, y, 5*DCURXs, DCURYs, 1);
    }

    b->xuse = 1;
    return;
}

void
make_box_list(BoxList *b, char *wname, char *iname, int32 n, int32 type,
              int32 use) {
    int32 nrow;
    int32 i;
    char sss[256];
    double z;

    nrow = 10;

    if (n < 10)
        nrow = n;
    b->xuse = 0;
    b->use = use;
    b->mc = 21;
    b->type = type;
    b->n = n;
    b->n0 = 0;
    b->nwin = nrow;
    b->value = xmalloc((usize)n*sizeof(char *));
    b->pos = xmalloc((usize)n*sizeof(*(b->pos)));
    b->off = xmalloc((usize)n*sizeof(*(b->off)));
    b->iname = xmalloc((usize)strlen(iname) + 5);
    strcpy(b->iname, iname);
    b->wname = xmalloc(strlen(wname) + 5);
    strcpy(b->wname, wname);

    for (i = 0; i < n; i++) {
        b->value[i] = xmalloc(256);
        switch (type) {
        case PARAMBOX:
            get_val(upar_names[i], &z);
            sprintf(sss, "%.16g", z);
            init_conds_set_edit_params(b, i, sss);
            break;
        case ICBOX:
            sprintf(sss, "%.16g", last_ic[i]);
            init_conds_set_edit_params(b, i, sss);
            break;
        case BCBOX:
            init_conds_set_edit_params(b, i, my_bc[i].string);
            break;
        case DELAYBOX:
            init_conds_set_edit_params(b, i, delay_string[i]);
            break;
        default:
            fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
            exit(EXIT_FAILURE);
        }
    }
    return;
}

/* this is added to take care of making sure
    exposure of the boxes is easily taken care of
 */

void
init_conds_do_box_expose(Window window) {
    if (ICBox.xuse)
        display_box(ICBox, window);
    if (BCBox.xuse)
        display_box(BCBox, window);
    if (ParamBox.xuse)
        display_box(ParamBox, window);
    if (DelayBox.xuse)
        display_box(DelayBox, window);
    return;
}

void
justify_string(Window w1, char *s1) {
    int32 n1 = (int32)strlen(s1)*DCURXs, nt = 10*DCURXs;
    int32 i = 0;
    if (n1 < nt)
        i = nt - n1;
    XClearWindow(display, w1);
    XDrawString(display, w1, small_gc, i, CURY_OFFs, s1, (int)strlen(s1));
    return;
}

/* new code is a bit tricky here - we dont want
 * to draw it if it is not visible
 *  there are nwin windows covering indexes
 *  n0,n0+1,...n0+nwin-1
 *  if the index is beyond this dont draw it
 */

void
init_conds_draw_one_box(BoxList b, int32 index) {
    Window window;
    Window we;

    int32 n0 = b.n0;
    int32 n1 = n0 + b.nwin - 1;
    int32 i;
    if (b.xuse == 0)
        return;
    if (index < n0 || index > n1)
        return; /* don't draw the ones out of range*/
    i = index - n0;
    window = b.w[i];
    we = b.we[i];
    switch (b.type) {
    case PARAMBOX:
        draw_editable(we, b.value[index], b.off[index], b.pos[index], b.mc);
        justify_string(window, upar_names[index]);
        break;
    case BCBOX:
        justify_string(window, my_bc[index].name);
        draw_editable(we, b.value[index], b.off[index], b.pos[index], b.mc);
        break;
    case ICBOX:
        draw_editable(we, b.value[index], b.off[index], b.pos[index], b.mc);
        justify_string(window, uvar_names[index]);
        break;
    case DELAYBOX:
        justify_string(window, uvar_names[index]);
        draw_editable(we, b.value[index], b.off[index], b.pos[index], b.mc);
        break;
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    return;
}

void
init_conds_redraw_params(void) {
    int32 i;
    double z;
    derived_evaluate();
    if (ParamBox.use)
        for (i = 0; i < NUPAR; i++) {
            get_val(upar_names[i], &z);
            add_edit_float(&ParamBox, i, z);
            init_conds_draw_one_box(ParamBox, i);
        }
    init_conds_reset_sliders();
}

void
init_conds_redraw_delays(void) {
    int32 i;
    if (DelayBox.use)
        for (i = 0; i < NODE; i++)
            init_conds_draw_one_box(DelayBox, i);
    return;
}

void
init_conds_redraw_ics(void) {
    int32 i;
    int32 in;
    for (i = 0; i < NODE + NMarkov; i++) {
        add_edit_float(&ICBox, i, last_ic[i]);
        init_conds_draw_one_box(ICBox, i);
    }
    init_conds_reset_sliders();
    if (ICBox.xuse == 0)
        return;
    for (i = 0; i < ICBox.nwin; i++) {
        in = i + ICBox.n0;
        if (ICBox.isck[in])
            XDrawString(display, ICBox.ck[i], small_gc, 0, CURY_OFFs, "*", 1);
        else
            XClearWindow(display, ICBox.ck[i]);
    }
    return;
}

void
init_conds_redraw_bcs(void) {
    int32 i;
    for (i = 0; i < NODE; i++)
        init_conds_draw_one_box(BCBox, i);
    return;
}

void
display_box(BoxList b, Window window) {
    int32 i;
    int32 n0 = b.n0;
    int32 n1 = n0 + b.nwin;
    int32 index;
    if (b.xuse == 0)
        return;
    if (b.close == window)
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Close", 5);
    if (b.go == window)
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Go", 2);
    if (b.ok == window)
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Ok", 2);
    if (b.cancel == window)
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Cancel", 6);
    if (b.def == window)
        XDrawString(display, window, small_gc, 5, CURY_OFFs, "Default", 7);
    if (b.up == window) {
        /*XDrawString(display,w,small_gc,5+DCURX,CURY_OFFs,
        "^",1);
        */
    }
    if (b.dn == window) {
        /*XDrawString(display,w,small_gc,5+DCURX,CURY_OFFs,
        "v",1);*/
    }
    if (b.pgup == window) {
        /*XDrawString(display,w,small_gc,5,CURY_OFFs,
        "^^",2);*/
    }
    if (b.pgdn == window) {
        /*XDrawString(display,w,small_gc,5,CURY_OFFs,
        "vv",2);*/
    }
    if (b.type == ICBOX) {
        if (b.xvt == window)
            XDrawString(display, window, small_gc, 3, CURY_OFFs, "xvst", 4);
        if (b.pp == window)
            XDrawString(display, window, small_gc, 3, CURY_OFFs, "xvsy", 4);
        if (b.arr == window)
            XDrawString(display, window, small_gc, 3, CURY_OFFs, "arry", 4);
    }

    for (i = 0; i < b.nwin; i++)
        if (b.w[i] == window || b.we[i] == window) {
            init_conds_draw_one_box(b, i + b.n0);
            return;
        }
    if (b.type == ICBOX) {
        for (i = 0; i < b.nwin; i++) {
            index = i + b.n0;
            if (index >= n0 && index < n1) {
                if (b.ck[i] == window && b.isck[index] == 1)
                    XDrawString(display, window, small_gc, 5, CURY_OFFs, "*",
                                1);
            }
        }
    }
    return;
}

void
init_conds_box_enter_events(Window window, int32 yn) {
    int32 i;
    int32 val;
    if (yn == 1)
        val = 2;
    else
        val = 1;
    if (ICBox.xuse)
        box_enter(ICBox, window, val);
    if (BCBox.xuse)
        box_enter(BCBox, window, val);
    if (ParamBox.xuse)
        box_enter(ParamBox, window, val);
    if (DelayBox.xuse)
        box_enter(DelayBox, window, val);
    if (ICBox.xuse &&
        (window == ICBox.xvt || window == ICBox.pp || window == ICBox.arr))
        XSetWindowBorderWidth(display, window, (uint)val);
    if (ICBox.xuse == 0)
        return;
    for (i = 0; i < ICBox.nwin; i++)
        if (window == ICBox.ck[i])
            XSetWindowBorderWidth(display, window, (uint)val);
    return;
}

void
box_enter(BoxList b, Window window, int32 val) {
    Window w = window;
    if (w == b.ok || w == b.cancel || w == b.def || w == b.go || w == b.close ||
        w == b.dn || w == b.up || w == b.pgdn || w == b.pgup)
        XSetWindowBorderWidth(display, w, (uint)val);
    return;
}

void
redraw_entire_box(BoxList *b) {
    if (b->xuse == 0)
        return;
    switch (b->type) {
    case PARAMBOX:
        init_conds_redraw_params();
        return;
    case BCBOX:
        init_conds_redraw_bcs();
        return;
    case ICBOX:
        init_conds_redraw_ics();
        return;
    case DELAYBOX:
        init_conds_redraw_delays();
        return;
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
}

void
do_box_button(BoxList *b, Window window) {
    int32 n = b->nwin;
    if (b->xuse == 0)
        return;
    if (window == b->close) {
        destroy_box(b);
        return;
    }
    if (window == b->ok || window == b->go)
        load_entire_box(b);
    if (window == b->cancel)
        redraw_entire_box(b);
    if (window == b->go)
        integrate_run_now();
    if (window == b->def && b->type == PARAMBOX)
        set_default_params();
    if (window == b->def && b->type == ICBOX)

        /* set default ics */
        for (int32 i2 = 0; i2 < NODE + NMarkov; i2++)
            last_ic[i2] = default_ic[i2];
    init_conds_redraw_ics();

    /* now for the "scrolling"

     */
    if (window == b->up)
        box_list_scroll(b, 1);
    if (window == b->pgup)
        box_list_scroll(b, b->nwin);
    if (window == b->dn)
        box_list_scroll(b, -1);
    if (window == b->pgdn)
        box_list_scroll(b, -b->nwin);

    for (int32 i = 0; i < n; i++) {
        if (window == b->we[i]) {
            XSetInputFocus(display, window, RevertToParent, CurrentTime);
            do {
                /* check box cursor */
                int32 n0;
                if (HotBoxItem < 0 || HotBox->xuse == 0)
                    break;
                n0 = HotBox->n0;
                draw_editable(HotBox->we[HotBoxItem],
                              HotBox->value[HotBoxItem + n0],
                              HotBox->off[HotBoxItem], HotBox->pos[HotBoxItem],
                              HotBox->mc);
                HotBoxItem = -1;
            } while (0);
            HotBoxItem = i;
            HotBox = b;
            draw_editable(window, b->value[i + b->n0], b->off[i + b->n0],
                          b->pos[i + b->n0], b->mc);
        }
    }

    if (b->type == ICBOX) {
        for (int32 i = 0; i < b->nwin; i++) {
            if (window == b->ck[i]) {
                b->isck[i + b->n0] = 1 - b->isck[i + b->n0];
                if (b->isck[i + b->n0])
                    XDrawString(display, window, small_gc, 0, CURY_OFFs, "*",
                                1);
                else
                    XClearWindow(display, window);
            }
        }
        if (window == b->xvt) {
            /* set up xvt */
            int32 plot_list[10];
            int32 n2 = 0;
            for (int32 i = 0; i < ICBox.n; i++)
                if (ICBox.isck[i]) {
                    if (n2 < 10) {
                        plot_list[n2] = i + 1;
                        n2++;
                    }
                    ICBox.isck[i] = 0;
                }
            for (int32 i = 0; i < ICBox.nwin; i++)
                XClearWindow(display, ICBox.ck[i]);
            if (n2 > 0)
                graf_par_graph_all(plot_list, n2, 0);
        }
        if (window == b->pp) {
            /* set up pp */
            int32 plot_list[3], n2 = 0;

            for (int32 i = 0; i < ICBox.n; i++)
                if (ICBox.isck[i]) {
                    if (n2 < 3) {
                        plot_list[n2] = i + 1;
                        n2++;
                    }
                    ICBox.isck[i] = 0;
                }
            for (int32 i = 0; i < ICBox.nwin; i++)
                XClearWindow(display, ICBox.ck[i]);
            if (n2 > 1)
                graf_par_graph_all(plot_list, n2, 1);
        }

        if (window == b->arr) {
            /* set up arry */
            int32 plot_list[2], n2 = 0;

            for (int32 i = 0; i < ICBox.n; i++)
                if (ICBox.isck[i]) {
                    if (n2 < 2) {
                        plot_list[n2] = i + 1;
                        n2++;
                    }
                    ICBox.isck[i] = 0;
                }
            for (int32 i = 0; i < ICBox.nwin; i++)
                XClearWindow(display, ICBox.ck[i]);
            if (n2 == 2)
                array_plot_optimize(plot_list);
        }
    }
    return;
}

void
box_list_scroll(BoxList *b, int32 i) {
    int32 n0 = b->n0;
    int32 new;
    int32 nw = b->nwin, n = b->n;
    int32 nend;
    if (n <= nw)
        return; /* do nothing - there is nothing to do */
    new = n0 - i;
    nend = new + nw;
    if (new < 0)
        new = 0;
    if (nend > n)
        new = n - nw;
    b->n0 = new;
    switch (b->type) {
    case PARAMBOX:
        load_entire_box(b);
        init_conds_redraw_params();
        init_conds_reset_sliders();
        break;
    case BCBOX:
        load_entire_box(b);
        init_conds_redraw_bcs();

        break;
    case ICBOX:
        load_entire_box(b);
        init_conds_redraw_ics();
        init_conds_reset_sliders();
        break;
    case DELAYBOX:
        load_entire_box(b);
        init_conds_redraw_delays();
        break;
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    return;
}

void
init_conds_box_buttons(Window window) {
    if (ICBox.xuse)
        do_box_button(&ICBox, window);
    if (BCBox.xuse)
        do_box_button(&BCBox, window);
    if (DelayBox.xuse)
        do_box_button(&DelayBox, window);
    if (ParamBox.xuse)
        do_box_button(&ParamBox, window);
    return;
}

void
init_conds_box_keypress(XEvent event, int32 *used) {
    if (ICBox.xuse) {
        do_box_key(&ICBox, event, used);
        if (*used)
            return;
    }
    if (BCBox.xuse) {
        do_box_key(&BCBox, event, used);
        if (*used)
            return;
    }
    if (DelayBox.xuse) {
        do_box_key(&DelayBox, event, used);
        if (*used)
            return;
    }
    if (ParamBox.xuse) {
        do_box_key(&ParamBox, event, used);
        if (*used)
            return;
    }
    return;
}

void
do_box_key(BoxList *b, XEvent event, int32 *used) {
    Window window = event.xkey.window;
    char ch;
    Window focus;
    int32 rev, n = b->nwin, i, j, flag;
    *used = 0;
    if (b->xuse == 0)
        return;
    for (i = 0; i < n; i++) {
        if (b->we[i] == window) {
            XGetInputFocus(display, &focus, &rev);
            if (window == focus) {
                *used = 1;
                ch = (char)ggets_get_key_press(&event);
                flag = edit_bitem(b, i, ch);
                if (flag == EDIT_NEXT && n > 1) {
                    j = i + 1;
                    if (j == n)
                        j = 0;
                    XSetInputFocus(display, b->we[j], RevertToParent,
                                   CurrentTime);
                    set_value_from_box(b, i);

                    HotBoxItem = j;
                    draw_editable(b->we[i], b->value[i + b->n0],
                                  b->off[i + b->n0], b->pos[i + b->n0], b->mc);
                    draw_editable(b->we[j], b->value[j + b->n0],
                                  b->off[j + b->n0], b->pos[j + b->n0], b->mc);
                    if (b->type == PARAMBOX || b->type == ICBOX)
                        init_conds_reset_sliders();
                }

                if (flag == EDIT_DONE) {
                    HotBoxItem = -1;
                    XSetInputFocus(display, main_win, RevertToParent,
                                   CurrentTime);
                    load_entire_box(b);
                }
                if (flag == EDIT_ESC) {
                    HotBoxItem = -1;
                    XSetInputFocus(display, main_win, RevertToParent,
                                   CurrentTime);
                }
            }
        }
    }
    return;
}

void
init_conds_man_ic(void) {
    int32 done, index = 0;
    double z;
    char name[256], junk[256];
    while (true) {
        sprintf(name, "%s :", uvar_names[index]);
        z = last_ic[index];
        done = ggets_new_float(name, &z);
        if (done == 0) {
            last_ic[index] = z;
            sprintf(junk, "%.16g", z);
            init_conds_set_edit_params(&ICBox, index, junk);
            init_conds_draw_one_box(ICBox, index);
            index++;
            if (index >= NODE + NMarkov)
                break;
        }
        if (done == -1)
            break;
    }
    return;
}

void
init_conds_new_parameter(void) {
    int32 done;
    int32 index;
    double z;
    char name[256];
    char value[sizeof(name) + 2];
    char junk[256];
    while (true) {
        name[0] = 0;
        done = ggets_new_string("Parameter:", name);
        if (strlen(name) == 0 || done == 0) {
            init_conds_redo_stuff();
            return;
        }
        if (strncasecmp(name, "DEFAULT", 7) == 0) {
            set_default_params();
            continue;
        }

        if (strncasecmp(name, "!LOAD", 5) == 0) {
            lunch_io_parameter_file(name, READEM);
            continue;
        }
        if (strncasecmp(name, "!SAVE", 5) == 0) {
            lunch_io_parameter_file(name, WRITEM);
            continue;
        }

        else {
            index = init_conds_find_user_name(PARAMBOX, name);
            if (index >= 0) {
                get_val(upar_names[index], &z);
                snprintf(value, sizeof(value), "%s :", name);
                done = ggets_new_float(value, &z);
                if (done == 0) {
                    set_val(upar_names[index], z);
                    sprintf(junk, "%.16g", z);
                    init_conds_set_edit_params(&ParamBox, index, junk);
                    init_conds_draw_one_box(ParamBox, index);
                    init_conds_reset_sliders();
                }
                if (done == -1) {
                    init_conds_redo_stuff();
                    break;
                }
            }
        }
    }
    return;
}

void
init_conds_redo_stuff(void) {
    derived_evaluate();
    volterra_re_evaluate_kernels();
    tabular_redo_all_fun_tables();
    derived_evaluate();
    return;
}

void
set_default_params(void) {
    int32 i;
    char junk[256];
    for (i = 0; i < NUPAR; i++) {
        set_val(upar_names[i], default_val[i]);
        sprintf(junk, "%.16g", default_val[i]);
        init_conds_set_edit_params(&ParamBox, i, junk);
    }

    init_conds_redraw_params();
    volterra_re_evaluate_kernels();
    tabular_redo_all_fun_tables();
}

void
draw_editable(Window window, char *string, int32 off, int32 cursor, int32 mc)
/* cursor position in letters to the left */
/* first character of string is off */
{

    int32 l = (int32)strlen(string) - off, rev, cp;
    Window focus;
    if (l > mc)
        l = mc;
    XClearWindow(display, window);
    XDrawString(display, window, small_gc, 0, CURY_OFF, string + off, l);
    XGetInputFocus(display, &focus, &rev);
    if (focus == window) {
        cp = DCURXs*(cursor - off); /* must be fixed */
        put_edit_cursor(window, cp);
    }
    return;
}

void
put_edit_cursor(Window window, int32 pos) {
    int32 x1 = pos;
    int32 x2 = x1 + 1;
    XDrawLine(display, window, small_gc, x1, 1, x1, DCURYs - 1);
    XDrawLine(display, window, small_gc, x2, 1, x2, DCURYs - 1);
    return;
}

int32
edit_bitem(BoxList *b, int32 i, int32 ch) {
    Window window = b->we[i];
    int32 i0 = i + b->n0;
    char *string = b->value[i0];
    int32 off = b->off[i0];
    int32 pos = b->pos[i0];
    int32 mc = b->mc;
    int32 l = (int32)strlen(string), wpos = pos - off;

    switch (ch) {
    case KEY_LEFT:
        if (pos > 0) {
            pos--;
            wpos--;
            if (wpos < 0) {
                off = off - 4;
                if (off < 0)
                    off = 0;
                wpos = pos - off;
            }
        } else
            ggets_ping();
        break;
    case KEY_RIGHT:
        if (pos < l) {
            pos++;
            wpos++;
            if (wpos > mc) {
                off = off + 4;
                if (off + mc > l)
                    off = l - mc;
                wpos = pos - off;
            }
        } else
            ggets_ping();
        break;
    case KEY_HOME:
        pos = 0;
        wpos = 0;
        break;
    case KEY_END:
        pos = l;
        wpos = mc;
        break;
    case KEY_BADKEY:
        return 0;

    case KEY_DOWN:
        box_list_scroll(b, -1);
        return 0;
    case KEY_UP:
        box_list_scroll(b, 1);
        return 0;
    case KEY_PGUP:
        box_list_scroll(b, b->nwin);
        return 0;
    case KEY_PGDN:
        box_list_scroll(b, -b->nwin);
        return 0; /* junk key  */
    case KEY_ESC:
        return EDIT_ESC;
    case KEY_FINE:
        return EDIT_NEXT;
    case KEY_BKSP:
        /*
        if(pos<l){
          ggets_mem_mov(&string[pos],&string[pos+1],l-pos);
          l--;
        }
        else
         ggets_ping();
         break; */
    case KEY_DEL:

        if (pos > 0) {
            ggets_mem_mov(&string[pos - 1], &string[pos], l - pos + 1);
            pos--;
            wpos--;
            if (wpos < 0) {
                off = off - 4;
                if (off < 0)
                    off = 0;
                wpos = pos - off;
            }
            l--;
        } else
            ggets_ping();
        break;
    case KEY_TAB:
        return EDIT_DONE;
    default:
        if ((ch >= ' ') && (ch <= '~')) {
            if (strlen(string) >= 256)
                ggets_ping();
            else {
                ggets_mov_mem(&string[pos + 1], &string[pos], l - pos + 1);
                string[pos] = (char)ch;
                pos = pos + 1;
                wpos++;
                l++;
                if (wpos > mc) {
                    off = off + 4;
                    if (off + mc > l)
                        off = l - mc;
                    wpos = pos - off;
                }
            }
        }
        break;
    }
    /* all done lets save everything */
    off = pos - wpos;

    b->off[i0] = off;
    b->pos[i0] = pos;
    draw_editable(window, string, off, pos, mc);
    return 0;
}

void
add_edit_float(BoxList *b, int32 i, double z) {
    char junk[256];
    sprintf(junk, "%.16g", z);
    add_editval(b, i, junk);
    return;
}

void
init_conds_set_edit_params(BoxList *b, int32 i, char *string) {
    int32 l = (int32)strlen(string);
    strcpy(b->value[i], string);
    b->off[i] = 0;
    if (l > b->mc)
        b->pos[i] = b->mc;
    else
        b->pos[i] = l;
    return;
}

void
add_editval(BoxList *b, int32 i, char *string) {
    int32 n0 = b->n0, n1 = b->n0 + b->nwin - 1;
    int32 iw;
    init_conds_set_edit_params(b, i, string);
    if (i < n0 || i > n1)
        return;
    iw = i - n0;
    if (b->xuse)
        draw_editable(b->we[iw], string, b->off[i], b->pos[i], b->mc);
    return;
}

int32
to_float(char *s, double *z) {
    int32 flag;
    *z = 0.0;
    if (s[0] == '%') {
        flag = calc_do_calc(&s[1], z);
        if (flag == -1)
            return -1;
        return 0;
    }
    *z = atof(s);
    return 0;
}

void
set_value_from_box(BoxList *b, int32 i) {
    char *s;
    double z;
    s = b->value[i];
    switch (b->type) {
    case ICBOX:
        if (to_float(s, &z) == -1)
            return;
        last_ic[i] = z;
        add_edit_float(b, i, z);
        break;
    case PARAMBOX:
        if (to_float(s, &z) == -1)
            return;
        set_val(upar_names[i], z);
        add_edit_float(b, i, z);
        break;

    case BCBOX:
        strcpy(my_bc[i].string, s);
        add_editval(b, i, s);
        break;
    case DELAYBOX:
        strcpy(delay_string[i], s);
        add_editval(b, i, s);
        break;
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    return;
}

void
load_entire_box(BoxList *b) {
    int32 n = b->n;

    for (int32 i = 0; i < n; i++)
        set_value_from_box(b, i);
    if (b->type == PARAMBOX) {
        volterra_re_evaluate_kernels();
        tabular_redo_all_fun_tables();
        init_conds_reset_sliders();
    }
    if (b->type == DELAYBOX)
        delay_handle_do_init_delay(DELAY);
    return;
}
