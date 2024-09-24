#include "xpplim.h"
#include "functions.h"
#include "parserslow.h"
#include "integers.h"
#include <stdbool.h>

#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include "graph.bitmap"
#include "struct.h"

#define MAXLAB 50
#define MAXGROB 400
#define POINTER 0
#define ARROW 1
#define MARKER 2 /* markers start at 2  there are several of them */

#define WDMARK .001
#define HTMARK .0016

static struct Grob {
    double xs, ys, xe, ye;
    double size;
    int16 use;
    Window window;
    int32 type;
    int32 color;
} grob[MAXGROB];

static struct Label {
    Window window;
    double x;
    double y;
    char s[MAXCHAR];
    int16 use;
    int32 font;
    int32 size;
} lb[MAXLAB];

static struct MarkInfo {
    int32 type;
    int32 color;
    int32 number, start, skip;
    double size;
} markinfo = {2, 0, 1, 0, 1, 1.0};

int32 manual_expose = 0;

Graph graph[MAXPOP];
Curve frz[MAXFRZ];
Graph *MyGraph;

int32 SimulPlotFlag = 0;
int32 current_pop;
int32 num_pops;
static int32 MINI_H = 300;
static int32 MINI_W = 450;

int32 ActiveWinList[MAXPOP];

static void select_sym(Window window);
static void lo_lite(Window wi);
static void set_gr_back(void);
static void set_gr_fore(void);
static void select_window(Window window);
static void kill_all_pops(void);
static void destroy_a_pop(void);
static void set_restore(int32 flag);
static void add_pntarr(int32 type);
static void add_markers(void);
static void add_marker(void);
static int32 get_markers_info(void);
static int32 get_marker_info(void);
static int32 select_marker_type(int32 *type);
static void destroy_label(Window window);
static void destroy_grob(Window window);
static void arrow_head(double xs, double ys, double xe, double ye, double size);
static void draw_grob(int32 i);
static void draw_marker(double x, double y, double size, int32 type);

int32
many_pops_select_table(void) {
    int32 i;
    int32 j;
    Window temp = main_win;
    char *n[MAX_TAB], key[MAX_TAB], ch;
    for (i = 0; i < NTable; i++) {
        n[i] = xmalloc(25);
        key[i] = 'a' + (char)i;
        sprintf(n[i], "%c: %s", key[i], my_table[i].name);
    }
    key[NTable] = 0;
    ch = (char)pop_up_list(&temp, "Table", n, key, NTable, 12, 0, 10, 0,
                           no_hint, info_pop, info_message);
    for (i = 0; i < NTable; i++)
        free(n[i]);
    j = (int32)(ch - 'a');
    if (j < 0 || j >= NTable) {
        ggets_err_msg("Not a valid table");
        return -1;
    }
    return j;
}

void
many_pops_get_intern_set(void) {
    char *n[MAX_INTERN_SET], key[MAX_INTERN_SET], ch;
    int32 i;
    int32 j;
    int32 count = Nintern_set;
    Window temp = main_win;
    if (count == 0)
        return;
    for (i = 0; i < Nintern_set; i++) {
        n[i] = xmalloc(256);
        key[i] = 'a' + (char)i;
        sprintf(n[i], "%c: %s", key[i], intern_set[i].name);
    }
    key[count] = 0;
    ch = (char)pop_up_list(&temp, "Param set", n, key, count, 12, 0, 10, 0,
                           no_hint, info_pop, info_message);
    for (i = 0; i < count; i++)
        free(n[i]);
    j = (int32)(ch - 'a');
    if (j < 0 || j >= Nintern_set) {
        ggets_err_msg("Not a valid set");
        return;
    }
    graphics_get_graph();
    load_eqn_extract_internset(j);
    numerics_chk_delay();
    init_conds_redraw_params();
    init_conds_redraw_ics();
    graphics_reset_graph();
}

void
many_pops_make_icon(char *icon, int32 wid, int32 hgt, Window window) {
    Pixmap icon_map;
    XWMHints wm_hints;
    icon_map =
        XCreateBitmapFromData(display, window, icon, (uint)wid, (uint)hgt);
    wm_hints.initial_state = NormalState;
    wm_hints.input = True;
    wm_hints.icon_pixmap = icon_map;
    wm_hints.flags = StateHint | IconPixmapHint | InputHint;

    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";
        XSetWMProperties(display, window, NULL, NULL, NULL, 0, NULL, &wm_hints,
                         &class_hints);
    }
    return;
}

void
many_pops_title_text(char *string) {
    many_pops_gtitle_text(string, draw_win);
    return;
}

void
many_pops_gtitle_text(char *string, Window win) {
    XTextProperty wname;
    XTextProperty iname;
    many_pops_gr_col();
    if (win != graph[0].window) {
        XStringListToTextProperty(&string, 1, &wname);
        XStringListToTextProperty(&string, 1, &iname);
        XSetWMProperties(display, win, &wname, &iname, NULL, 0, NULL, NULL,
                         NULL);
    } else {
        int32 len = (int32)strlen(string);
        int32 x, y;
        uint32 w, h, bw, de;
        int32 xs, ys = 2;
        Window root;
        XGetGeometry(display, win, &root, &x, &y, &w, &h, &bw, &de);
        xs = ((int32)w - len*DCURX) / 2;
        if (xs < 0)
            xs = 0;
        ggets_f_text(xs, ys, string, win);
        color_set(0);
        ggets_xline(0, 18, (int32)w, 18, win);
    }
    many_pops_base_col();
    return;
}

void
many_pops_restore_off(void) {
    MyGraph->Restore = 0;
    /* MyGraph->Nullrestore=0; */
    return;
}

void
many_pops_restore_on(void) {
    MyGraph->Restore = 1;
    /*  MyGraph->Nullrestore=1; */
    return;
}

void
many_pops_add_label(char *s, int32 x, int32 y, int32 size, int32 font) {
    int32 i;
    double xp;
    double yp;
    scale_to_real(x, y, &xp, &yp);
    for (i = 0; i < MAXLAB; i++) {
        if (lb[i].use == 0) {
            lb[i].use = 1;
            lb[i].x = xp;
            lb[i].y = yp;
            lb[i].window = draw_win;
            lb[i].font = font;
            lb[i].size = size;
            strcpy(lb[i].s, s);
            return;
        }
    }
    return;
}

void
draw_marker(double x, double y, double size, int32 type) {
    int32 pen = 0;
    double x1 = x, y1 = y, x2, y2;
    int32 ind = 0;
    int32 offset;

    static int32 sym_dir[] = {
        /*          box              */
        0,
        -6,
        -6,
        1,
        12,
        0,
        1,
        0,
        12,
        1,
        -12,
        0,
        1,
        0,
        -12,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,

        /*          diamond             */
        0,
        8,
        0,
        1,
        -8,
        -8,
        1,
        8,
        -8,
        1,
        8,
        8,
        1,
        -8,
        8,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        /*          triangle         */
        0,
        -6,
        -6,
        1,
        12,
        0,
        1,
        -6,
        12,
        1,
        -6,
        -12,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,

        /*          plus            */
        0,
        -6,
        0,
        1,
        12,
        0,
        0,
        -6,
        -6,
        1,
        0,
        12,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,

        /*          cross            */
        0,
        -6,
        6,
        1,
        12,
        -12,
        0,
        -12,
        0,
        1,
        12,
        12,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,

        /*          circle           */
        0,
        6,
        0,
        1,
        -1,
        3,
        1,
        -2,
        2,
        1,
        -3,
        1,
        1,
        -3,
        -1,
        1,
        -2,
        -2,
        1,
        -1,
        -3,
        1,
        1,
        -3,
        1,
        2,
        -2,
        1,
        3,
        -1,
        1,
        3,
        1,
        1,
        2,
        2,
        1,
        1,
        3,
        3,
        0,
        0,
        3,
        0,
        0,
        3,
        0,
        0,

    };
    double dx = (MyGraph->xhi - MyGraph->xlo)*WDMARK*size;
    double dy = (MyGraph->yhi - MyGraph->ylo)*HTMARK*size;
    while (true) {
        offset = 48*type + 3*ind;
        pen = sym_dir[offset];
        if (pen == 3)
            break;
        x2 = dx*sym_dir[offset + 1] + x1;
        y2 = dy*sym_dir[offset + 2] + y1;
        if (pen == 1)
            graphics_line_abs(x1, y1, x2, y2);
        x1 = x2;
        y1 = y2;
        ind++;
    }
    return;
}

void
draw_grob(int32 i) {
    double xs = grob[i].xs, ys = grob[i].ys, xe = grob[i].xe, ye = grob[i].ye;
    graphics_set_linestyle(grob[i].color);
    if (grob[i].type == POINTER)
        graphics_line_abs(xs, ys, xe, ye);
    if (grob[i].type == ARROW || grob[i].type == POINTER)
        arrow_head(xs, ys, xe, ye, grob[i].size);
    if (grob[i].type >= MARKER)
        draw_marker(xs, ys, grob[i].size, grob[i].type - 2);
    return;
}

void
arrow_head(double xs, double ys, double xe, double ye, double size) {
    double l = xe - xs, h = ye - ys;
    double ar = (MyGraph->xhi - MyGraph->xlo) / (MyGraph->yhi - MyGraph->ylo);
    double x0 = xs + size*l, y0 = ys + size*h;
    /* double tot=(double)sqrt((double)(l*l+h*h)); */

    double xp = x0 + .5*size*h*ar, yp = y0 - .5*size*l / ar;
    double xm = x0 - .5*size*h*ar, ym = y0 + .5*size*l / ar;
    graphics_line_abs(xs, ys, xp, yp);
    graphics_line_abs(xs, ys, xm, ym);
    return;
}

void
destroy_grob(Window window) {
    int32 i;
    for (i = 0; i < MAXGROB; i++) {
        if ((grob[i].use == 1) && (grob[i].window == window)) {
            grob[i].use = 0;
            grob[i].window = (Window)0;
        }
    }
    return;
}

void
destroy_label(Window window) {
    int32 i;
    for (i = 0; i < MAXLAB; i++) {
        if ((lb[i].use == 1) && (lb[i].window == window)) {
            lb[i].use = 0;
            lb[i].window = (Window)0;
        }
    }
    return;
}

void
many_pops_draw_label(Window window) {
    int32 i;
    many_pops_gr_col();
    for (i = 0; i < MAXLAB; i++) {
        if ((lb[i].use == 1) && (lb[i].window == window))
            graphics_fancy_text_abs(lb[i].x, lb[i].y, lb[i].s, lb[i].size);
    }
    for (i = 0; i < MAXGROB; i++) {
        if ((grob[i].use == 1) && (grob[i].window == window))
            draw_grob(i);
    }
    many_pops_base_col();
    return;
}

void
add_grob(double xs, double ys, double xe, double ye, double size, int32 type,
         int32 color) {
    int32 i;
    for (i = 0; i < MAXGROB; i++) {
        if (grob[i].use == 0) {
            grob[i].use = 1;
            grob[i].xs = xs;
            grob[i].xe = xe;
            grob[i].ys = ys;
            grob[i].ye = ye;
            grob[i].window = draw_win;
            grob[i].size = size;
            grob[i].color = color;
            grob[i].type = type;
            /*     main_redraw_all(); */
            return;
        }
    }
    return;
}

int32
select_marker_type(int32 *type) {
    int32 ival = *type - MARKER;
    int32 i;
    char *list[] = {"Box", "Diamond", "Triangle", "Plus", "X", "circle2"};
    static char key[] = "bdtpxc";
    Window temp = main_win;
    char ch;
    ch = (char)pop_up_list(&temp, "Markers", list, key, 6, 9, ival, 10,
                           4*DCURY + 8, no_hint, info_pop, info_message);
    if (ch == 27)
        return 0;
    for (i = 0; i < 6; i++) {
        if (ch == key[i])
            ival = i;
    }
    if (ival < 6)
        *type = MARKER + ival;

    return 1;
}

int32
get_marker_info(void) {
    static char *n[] = {"*5Type", "*4Color", "Size"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    int32 status;
    snprintf(values[0], sizeof(values[0]), "%d", markinfo.type);
    snprintf(values[1], sizeof(values[1]), "%d", markinfo.color);
    snprintf(values[2], sizeof(values[2]), "%g", markinfo.size);
    status = do_string_box(3, 3, 1, "Add Marker", n, values, 25);
    if (status != 0) {
        markinfo.type = atoi(values[0]);
        markinfo.size = atof(values[2]);
        markinfo.color = atoi(values[1]);
        return 1;
    }
    return 0;
}

int32
get_markers_info(void) {
    static char *n[] = {"*5Type", "*4Color", "Size", "Number", "Row1", "Skip"};
    char values[LENGTH(n)][MAX_LEN_SBOX];
    int32 status;
    snprintf(values[0], sizeof(values[0]), "%d", markinfo.type);
    snprintf(values[1], sizeof(values[1]), "%d", markinfo.color);
    snprintf(values[2], sizeof(values[2]), "%g", markinfo.size);
    snprintf(values[3], sizeof(values[3]), "%d", markinfo.number);
    snprintf(values[4], sizeof(values[4]), "%d", markinfo.start);
    snprintf(values[5], sizeof(values[5]), "%d", markinfo.skip);
    status = do_string_box(6, 6, 1, "Add Markers", n, values, 25);
    if (status != 0) {
        markinfo.type = atoi(values[0]);
        markinfo.size = atof(values[2]);
        markinfo.color = atoi(values[1]);
        markinfo.number = atoi(values[3]);
        markinfo.start = atoi(values[4]);
        markinfo.skip = atoi(values[5]);

        return 1;
    }
    return 0;
}

void
add_marker(void) {
    int32 flag, i1, j1, status;
    double xe = 0.0, ye = 0.0, xs, ys;
    status = get_marker_info();
    if (status == 0)
        return;
    menudrive_message_box("Position");
    flag = menudrive_get_mouse_xy(&i1, &j1);
    menudrive_message_box_kill();
    XFlush(display);
    if (flag == 0)
        return;
    scale_to_real(i1, j1, &xs, &ys);
    add_grob(xs, ys, xe, ye, markinfo.size, markinfo.type, markinfo.color);
    main_redraw_all();
}

void
add_markers(void) {
    int32 i;
    double xe = 0.0, ye = 0.0, xs, ys, x, y, z;

    if (get_markers_info() == 0)
        return;
    for (i = 0; i < markinfo.number; i++) {
        browse_get_data_xyz(&x, &y, &z, MyGraph->xv[0], MyGraph->yv[0], MyGraph->zv[0],
                     markinfo.start + i*markinfo.skip);
        if (MyGraph->ThreeDFlag == 0) {
            xs = x;
            ys = y;
        } else {
            graphics_threed_proj(x, y, z, &xs, &ys);
        }
        add_grob(xs, ys, xe, ye, markinfo.size, markinfo.type, markinfo.color);
    }
    main_redraw_all();
}

void
add_pntarr(int32 type) {
    double size = .1;
    int32 i1, j1, i2, j2, color = 0;
    double xe, ye, xs, ys;
    /*Window temp;*/
    int32 flag;
    /*temp=main_win;*/
    if (ggets_new_float("Size: ", &size))
        return;
    if (ggets_new_int("Color: ", &color))
        return;
    /* pop_list_message_box(&temp,0,SCALEY-5*DCURY,"Choose start/end"); */
    menudrive_message_box("Choose start/end");
    flag = rubber(&i1, &j1, &i2, &j2, draw_win, 1);
    /*  XDestroyWindow(display,temp); */
    menudrive_message_box_kill();
    XFlush(display);
    if (flag) {
        scale_to_real(i1, j1, &xs, &ys);
        scale_to_real(i2, j2, &xe, &ye);
        if (i1 == i2 && j1 == j2)
            return;
        add_grob(xs, ys, xe, ye, size, type, color);
        main_redraw_all();
    }
    return;
}

void
many_pops_edit_object_com(int32 com) {
    char ans, str[80];
    int32 i, j, ilab = -1, flag, type;
    double x;
    double y;
    double dist = 1e20, dd;

    menudrive_message_box("Choose Object");
    flag = menudrive_get_mouse_xy(&i, &j);

    menudrive_message_box_kill();
    XFlush(display);
    if (flag) {
        scale_to_real(i, j, &x, &y);
        /* now search all labels to find the best */
        type = 0; /* label =  0, arrows, etc =1 */
        for (i = 0; i < MAXLAB; i++) {
            if (lb[i].use == 1 && lb[i].window == draw_win) {
                dd = (x - lb[i].x)*(x - lb[i].x) +
                     (y - lb[i].y)*(y - lb[i].y);
                if (dd < dist) {
                    ilab = i;
                    dist = dd;
                }
            }
        }
        for (i = 0; i < MAXGROB; i++) {
            if (grob[i].use == 1 && grob[i].window == draw_win) {
                dd = (x - grob[i].xs)*(x - grob[i].xs) +
                     (y - grob[i].ys)*(y - grob[i].ys);
                if (dd < dist) {
                    ilab = i;
                    dist = dd;
                    type = 1;
                }
            }
        }
        if (ilab >= 0 && type == 0) {
            switch (com) {
            case 0:
                snprintf(str, sizeof(str), "Move %s ?", lb[ilab].s);
                ans = (char)menudrive_two_choice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    menudrive_message_box("Click on new position");
                    flag = menudrive_get_mouse_xy(&i, &j);

                    menudrive_message_box_kill();
                    XFlush(display);
                    if (flag) {
                        scale_to_real(i, j, &x, &y);
                        lb[ilab].x = x;
                        lb[ilab].y = y;
                        main_clr_scrn();
                        main_redraw_all();
                    }
                }
                break;
            case 1:
                snprintf(str, sizeof(str), "Change %s ?", lb[ilab].s);
                ans = (char)menudrive_two_choice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    ggets_new_string("Text: ", lb[ilab].s);
                    ggets_new_int("Size 0-4 :", &lb[ilab].size);
                    /* ggets_new_int("Font  0-times/1-symbol :",&lb[ilab].font);
                     */
                    if (lb[ilab].size > 4)
                        lb[ilab].size = 4;
                    if (lb[ilab].size < 0)
                        lb[ilab].size = 0;
                    main_clr_scrn();
                    main_redraw_all();
                }
                break;
            case 2:
                snprintf(str, sizeof(str), "Delete %s ?", lb[ilab].s);
                ans = (char)menudrive_two_choice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    lb[ilab].window = 0;
                    lb[ilab].use = 0;
                    main_clr_scrn();
                    main_redraw_all();
                }
                break;
            default:
                break;
            }
        }
        if (ilab >= 0 && type == 1) {
            switch (com) {
            case 0:
                snprintf(str, sizeof(str), "Move graphic at (%f,%f)",
                         grob[ilab].xs, grob[ilab].ys);
                ans = (char)menudrive_two_choice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    menudrive_message_box("Reposition");
                    flag = menudrive_get_mouse_xy(&i, &j);

                    menudrive_message_box_kill();
                    XFlush(display);
                    if (flag) {
                        scale_to_real(i, j, &x, &y);
                        grob[ilab].xe = grob[ilab].xe - grob[ilab].xs + x;
                        grob[ilab].ye = grob[ilab].ye - grob[ilab].ys + y;
                        grob[ilab].xs = x;
                        grob[ilab].ys = y;
                        main_clr_scrn();
                        main_redraw_all();
                    }
                }
                break;
            case 1:
                snprintf(str, sizeof(str), "Change graphic at (%f,%f)",
                         grob[ilab].xs, grob[ilab].ys);
                ans = (char)menudrive_two_choice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    if (grob[ilab].type >= MARKER)
                        select_marker_type(&grob[ilab].type);
                    ggets_new_float("Size ", &grob[ilab].size);
                    ggets_new_int("Color :", &grob[ilab].color);
                    main_clr_scrn();
                    main_redraw_all();
                }
                break;
            case 2:
                snprintf(str, sizeof(str), "Delete graphic at (%f,%f)",
                         grob[ilab].xs, grob[ilab].ys);
                ans = (char)menudrive_two_choice("Yes", "No", str, "yn");
                if (ans == 'y') {
                    grob[ilab].window = 0;
                    grob[ilab].use = 0;
                    main_clr_scrn();
                    main_redraw_all();
                }
                break;
            default:
                break;
            }
        }
    }
    return;
}

void
many_pops_do_gr_objs_com(int32 com) {
    switch (com) {
    case 0:
        ggets_cput_text();
        break;
    case 1:
        add_pntarr(ARROW);
        break;
    case 2:
        add_pntarr(POINTER);
        break;
    case 3:
        add_marker();
        break;
    case 6:
        add_markers();
        break;
    case 5:
        destroy_label(draw_win);
        destroy_grob(draw_win);
        main_clr_scrn();
        main_redraw_all();
        break;
    default:
        break;
    }
    return;
}

void
many_pops_set_active_windows(void) {
    int32 i, np = 0;
    for (i = 0; i < MAXPOP; i++) {
        if (graph[i].Use == 1) {
            ActiveWinList[np] = i;
            np++;
        }
    }
    num_pops = np;
    return;
}

void
many_pops_do_windows_com(int32 c) {
    switch (c) {
    case 0:
        many_pops_create_a_pop();
        break;
    case 1:
        if (pop_list_yes_no_box())
            kill_all_pops();
        break;
    case 3:
        XLowerWindow(display, draw_win);
        break;
    case 2:
        destroy_a_pop();
        break;
    case 5:
        set_restore(0);
        break;
    case 4:
        set_restore(1);
        break;
    case 6:
        SimulPlotFlag = 1 - SimulPlotFlag;
        break;
    default:
        break;
    }

    many_pops_set_active_windows();
    return;
}

void
set_restore(int32 flag) {
    int32 i;
    for (i = 0; i < MAXPOP; i++) {
        if (graph[i].window == draw_win) {
            graph[i].Restore = flag;
            graph[i].Nullrestore = flag;
            return;
        }
    }
    return;
}

void
destroy_a_pop(void) {
    int32 i;
    if (draw_win == graph[0].window) {
        pop_list_respond_box("Okay", "Can't destroy big window!");
        /*pop_list_respond_box(main_win,0,0,"Okay","Can't destroy big
         * window!");*/
        return;
    }
    for (i = 1; i < MAXPOP; i++)
        if (graph[i].window == draw_win)
            break;
    if (i >= MAXPOP)
        return;
    select_window(graph[0].window);
    graph[i].Use = 0;
    destroy_label(graph[i].window);
    destroy_grob(graph[i].window);
    browse_wait_a_sec(ClickTime);
    XDestroySubwindows(display, graph[i].window);
    XDestroyWindow(display, graph[i].window);
    num_pops--;
    return;
}

void
many_pops_init_grafs(int32 x, int32 y, int32 w, int32 h) {
    int32 i;
    many_pops_gr_col();
    for (i = 0; i < MAXLAB; i++) {
        lb[i].use = 0;
        lb[i].window = (Window)0;
    }
    for (i = 0; i < MAXGROB; i++) {
        grob[i].window = (Window)0;
        grob[i].use = 0;
    }
    graf_par_init_bd();
    for (i = 0; i < MAXFRZ; i++)
        frz[i].use = 0;
    /* for(i=0;i<MAXNCLINE ... */
    for (i = 0; i < MAXPOP; i++)
        graph[i].Use = 0;
    ActiveWinList[0] = 0;
    graphics_init_all();

    graph[0].window = XCreateSimpleWindow(display, main_win, x, y + 4, (uint)w,
                                          (uint)h, 2, GrFore, MyDrawWinColor);
    graph[0].w_info = info_pop;

    info_message = graph[0].gr_info;
    graph[0].Use = 1;
    graph[0].Restore = 1;
    graph[0].Nullrestore = 1;
    graph[0].x0 = x;
    graph[0].y0 = y;
    graph[0].Height = h;
    graph[0].Width = w;
    XSelectInput(display, graph[0].window,
                 KeyPressMask | ButtonPressMask | ExposureMask |
                     ButtonReleaseMask | ButtonMotionMask |
                     StructureNotifyMask);
    num_pops = 1;
    XMapWindow(display, graph[0].window);
    draw_win = graph[0].window;
    current_pop = 0;
    graphics_get_draw_area();
    select_sym(graph[0].window);
    many_pops_base_col();
}

void
many_pops_ps_restore(void) {
    if (Xup) {
        nullcline_redraw_dfield();
        ps_do_color(0);
        if (MyGraph->Nullrestore) {
            restore_nullclines();
            ps_stroke();
        }
    }
    ps_last_pt_off();

    integrate_restore(0, my_browser.maxrow);

    nullcline_do_batch_nclines();
    nullcline_do_batch_dfield();
    axes2_do();

    ps_do_color(0);
    if (Xup) {
        many_pops_draw_label(draw_win);
        graf_par_draw_freeze(draw_win);
    }
    ps_end();
    return;
}

void
many_pops_svg_restore(void) {
    /* integrate_restore(0,my_browser.maxrow);
     */
    /*ps_do_color(0);
    if(MyGraph->Nullrestore){restore_nullclines();ps_stroke();}
     */

    nullcline_redraw_dfield();
    if (MyGraph->Nullrestore)
        restore_nullclines();
    svg_last_pt_off();
    /*ps_do_color(0);*/
    integrate_restore(0, my_browser.maxrow);
    axes2_do();
    if (Xup) {
        many_pops_draw_label(draw_win);
        graf_par_draw_freeze(draw_win);
    }
    nullcline_do_batch_nclines();
    nullcline_do_batch_dfield();
    svg_end();
    return;
}

int32
many_pops_rotate_3dcheck(XEvent event) {
    Window window = event.xbutton.window;
    XEvent z;
    int32 xini, yini, dx, dy;
    double theta;
    double phi;
    if (window == draw_win && MyGraph->ThreeDFlag) {
        xini = event.xbutton.x;
        yini = event.xbutton.y;
        phi = MyGraph->Phi;
        theta = MyGraph->Theta;
        while (true) {
            XNextEvent(display, &z);
            if (z.type == ButtonRelease) {
                axes2_do();
                main_redraw_all();
                many_pops_hi_lite(draw_win);
                return 1;
            }
            if (z.type == MotionNotify) {
                dx = z.xmotion.x - xini;
                dy = z.xmotion.y - yini;
                MyGraph->Phi = phi - (double)dy;
                MyGraph->Theta = theta - (double)dx;
                axes2_redraw_cube_pt(MyGraph->Theta, MyGraph->Phi);
            }
        }
    }
    return 0;
}

void
many_pops_do_motion_events(XEvent event) {
    int32 i = event.xmotion.x;
    int32 j = event.xmotion.y;
    double x;
    double y;
    char buf[256];
    init_conds_slider_motion(event);
#ifdef AUTO
    auto_x11_motion(event);
#endif
    if (event.xmotion.window == draw_win) {
        scale_to_real(i, j, &x, &y);
        snprintf(buf, sizeof(buf), "x=%f y=%f ", x, y);
        many_pops_canvas_xy(buf);
    }
    return;
}

void
many_pops_do_expose(XEvent event) {
    int32 i;
    int32 cp = current_pop;
    Window temp;

    temp = draw_win;
    main_top_button_draw(event.xany.window);
    array_plot_expose(event.xany.window);
    /* txt_redraw_view(ev.xany.window);  */
    ani_expose(event.xany.window);
    expose_my_browser(event);
    /* ggets_draw_info_pop(ev.xany.window); */
    menudrive_message_box_redraw(event.xany.window);
    eig_list_draw_eq_list(event.xany.window);
    eig_list_draw_eq_box(event.xany.window);
    init_conds_do_box_expose(event.xany.window);
    init_conds_expose_slides(event.xany.window);
    menu_expose(event.xany.window);
#ifdef AUTO
    auto_x11_display(event.xany.window);
#endif
    /* if(ev.xexpose.window==menu_pop){
          menu_draw_help();

     return;
      }
      */
    if (manual_expose == 0) {
        many_pops_gr_col();

        for (i = 0; i < MAXPOP; i++) {
            if ((graph[i].Use) && (event.xexpose.window == graph[i].w_info)) {
                XClearWindow(display, graph[i].w_info);
                if (i == 0) {
                    many_pops_base_col();
                    XDrawString(display, graph[i].w_info, gc, 5, CURY_OFF,
                                graph[i].gr_info,
                                (int)strlen(graph[i].gr_info));
                } else {
                    many_pops_small_base();
                    XDrawString(display, graph[i].w_info, small_gc, 0,
                                CURY_OFFs, graph[i].gr_info,
                                (int)strlen(graph[i].gr_info));
                    many_pops_small_gr();
                }
            }
            if ((event.type == Expose) && (graph[i].Use) &&
                (event.xexpose.window == graph[i].window)) {
                /* nullcline_redraw_dfield(); */

                current_pop = i;
                MyGraph = &graph[i];
                draw_win = graph[i].window;
                graphics_get_draw_area();
                axes2_do();
                if (graph[i].Restore)
                    integrate_restore(0, my_browser.maxrow);
                many_pops_draw_label(graph[i].window);
                graf_par_draw_freeze(graph[i].window);
                if (graph[i].Nullrestore)
                    restore_nullclines();
            }
        }
    } /* namual expose */
    draw_win = temp;
    MyGraph = &graph[cp];
    current_pop = cp;
    many_pops_hi_lite(draw_win);
    graphics_get_draw_area();
    many_pops_base_col();
    many_pops_small_base();
    return;
}

void
many_pops_resize_all(int32 wid, int32 hgt) {
    int32 nw = wid - 16 - 16*DCURX + 7,
          nh = hgt - 3*DCURYb - 4*DCURYs - 24;
    nw = 4*((nw / 4));
    nh = 4*((nh / 4));
    XResizeWindow(display, graph[0].window, (uint)nw, (uint)nh);
    graph[0].Width = nw;
    graph[0].Height = nh;
    graphics_get_draw_area();
    return;
}

void
kill_all_pops(void) {
    int32 i;
    select_window(graph[0].window);

    for (i = 1; i < MAXPOP; i++)
        if (graph[i].Use) {
            graph[i].Use = 0;
            destroy_label(graph[i].window);
            destroy_grob(graph[i].window);

            XDestroySubwindows(display, graph[i].window);
            XDestroyWindow(display, graph[i].window);
        }
    num_pops = 1;
    return;
}

void
many_pops_create_a_pop(void) {
    int32 i;
    int32 index;

    for (i = 1; i < MAXPOP; i++)
        if (graph[i].Use == 0)
            break;
    if (i >= MAXPOP) {
        /*pop_list_respond_box(main_win,0,0,"Okay","Too many windows!");*/
        pop_list_respond_box("Okay", "Too many windows!");
        return;
    }
    index = i;

    graph[index].window =
        XCreateSimpleWindow(display, RootWindow(display, screen), 0, 0,
                            (uint)MINI_W, (uint)MINI_H, 2, GrFore, GrBack);
    graph[index].w_info =
        make_window(graph[index].window, 10, 0, 40*DCURXs, DCURYs, 0);
    XSetWindowBackground(display, graph[i].window, MyDrawWinColor);

    copy_graph(index, current_pop);
    graph[index].Width = MINI_W;
    graph[index].Height = MINI_H;
    graph[index].x0 = 0;
    graph[index].y0 = 0;
    num_pops++;
    many_pops_make_icon((char *)graph_bits, graph_width, graph_height,
                        graph[index].window);
    XSelectInput(display, graph[index].window,
                 KeyPressMask | ButtonPressMask | ExposureMask |
                     ButtonReleaseMask | ButtonMotionMask);
    XMapWindow(display, graph[index].window);
    XRaiseWindow(display, graph[index].window);
    XSetWMProtocols(display, graph[index].window, &deleteWindowAtom, 1);
    select_window(graph[index].window);
    /*  select_window(graph[0].window);
        select_window(graph[index].window); */
    XRaiseWindow(display, graph[0].window);
    /*  XDestroyWindow(display,temp); */
    return;
}

void
many_pops_gr_col(void) {
    XSetForeground(display, gc, GrFore);
    XSetBackground(display, gc, GrBack);
    return;
}

void
many_pops_base_col(void) {
    XSetForeground(display, gc, MyForeColor);
    XSetBackground(display, gc, MyBackColor);
    return;
}

void
many_pops_small_gr(void) {
    XSetForeground(display, small_gc, GrFore);
    XSetBackground(display, small_gc, GrBack);
    return;
}

void
many_pops_small_base(void) {
    XSetForeground(display, small_gc, MyForeColor);
    XSetBackground(display, small_gc, MyBackColor);
    return;
}

void
many_pops_change_plot_vars(int32 k) {
    int32 i;
    int32 ip;
    int32 np;
    for (i = 0; i < MAXPOP; i++) {
        if (graph[i].Use) {
            np = graph[i].nvars;
            for (ip = 0; ip < np; ip++) {
                if (graph[i].xv[ip] > k)
                    graph[i].xv[ip] = graph[i].xv[ip] - 1;
                if (graph[i].yv[ip] > k)
                    graph[i].yv[ip] = graph[i].yv[ip] - 1;
                if (graph[i].zv[ip] > k)
                    graph[i].zv[ip] = graph[i].zv[ip] - 1;
            }
        }
    }
    return;
}

int32
many_pops_check_active_plot(int32 k) {
    int32 i;
    int32 ip;
    int32 np;
    for (i = 0; i < MAXPOP; i++) {
        if (graph[i].Use) {
            np = graph[i].nvars;
            for (ip = 0; ip < np; ip++) {
                if (graph[i].xv[ip] == k || graph[i].yv[ip] == k ||
                    graph[i].zv[ip] == k)
                    return 1;
            }
        }
    }
    return 0;
}

void
many_pops_make_active(int32 i, int32 flag) {
    current_pop = i;
    MyGraph = &graph[current_pop];
    draw_win = MyGraph->window;
    graphics_get_draw_area_flag(flag);
    return;
}

void
select_window(Window window) {
    int32 i;

    if (window == draw_win)
        return;
    many_pops_gr_col();
    if (window == graph[0].window)
        current_pop = 0;
    else {
        for (i = 1; i < MAXPOP; i++)
            if ((graph[i].Use) && (window == graph[i].window))
                current_pop = i;
    }
    MyGraph = &graph[current_pop];
    lo_lite(draw_win);
    draw_win = window;
    many_pops_hi_lite(window);
    XRaiseWindow(display, window);
    graphics_get_draw_area();
    many_pops_base_col();
    return;
}

void
set_gr_fore(void) {
    XSetForeground(display, gc, GrFore);
    return;
}

void
set_gr_back(void) {
    XSetForeground(display, gc, GrBack);
    return;
}

void
many_pops_hi_lite(Window wi) {
    set_gr_fore();
    select_sym(wi);
    return;
}

void
lo_lite(Window wi) {
    set_gr_back();
    ggets_bar(0, 0, 5, 5, wi);
    return;
}

void
select_sym(Window window) {
    ggets_bar(0, 0, 5, 5, window);
    return;
}

void
many_pops_canvas_xy(char *buf) {
    XClearWindow(display, MyGraph->w_info);
    strcpy(MyGraph->gr_info, buf);
    if (MyGraph->w_info == info_pop) {
        many_pops_base_col();
        XDrawString(display, MyGraph->w_info, gc, 5, CURY_OFF, buf,
                    (int)strlen(buf));
    } else {
        many_pops_small_base();
        XDrawString(display, MyGraph->w_info, small_gc, 0, CURY_OFFs, buf,
                    (int)strlen(buf));
        /* many_pops_small_gr(); */
    }
    return;
}

void
many_pops_check_draw_button(XEvent event) {
    int32 k;
    char buf[256];

    int32 button;
    int32 i;
    int32 j;
    double x, y;
    int32 flag = 0;
    Window window;
    button = (int32)event.xbutton.button;
    window = event.xbutton.window;
    i = event.xbutton.x;
    j = event.xbutton.y;
    if (button == 1) { /* select window   */

        for (k = 1; k < MAXPOP; k++)
            if ((graph[k].Use) && (window == graph[k].window))
                flag = 1;
        if ((window == graph[1].window) || (flag == 1))
            select_window(window);
    } else /* any other button   */
    {
        if (window != draw_win)
            return;
        scale_to_real(i, j, &x, &y);
        snprintf(buf, sizeof(buf), "x=%f y=%f ", x, y);
        many_pops_canvas_xy(buf);
    }
}
