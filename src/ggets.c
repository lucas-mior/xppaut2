#include "functions.h"
#include "integers.h"
#include <stdbool.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/XKBlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include "mykeydef.h"
#include <stdarg.h>

#define XK_PgUp XK_Prior
#define XK_PgDn XK_Next

#define ESCAPE 27
char *info_message;
int32 ms_style = 0;

static int32 curs_x;
static int32 curs_y;

int32 xor_flag;

static void ggets_put_string_at(Window window, int32 col, char *s, int32 off);
static void ggets_clr_line_at(Window window, int32 col0, int32 pos, int32 n);

void
ggets_ping(void) {
    if (tfBell && !xpp_batch) {
        /*
        XkbBell allows window managers to react
        to bell events using possibly user-specified
        accessiblity options (e.g. visual bell)
        */

        XkbBell(display, command_pop, 100, (Atom)NULL);
    }
    /* Call to XBell seems to be ignored by many window managers
     * where XkbBell is not. if(tfBell&&!xpp_batch) XBell(display,100); */
    return;
}

void
ggets_reset_graphics(void) {
    ggets_blank_screen(draw_win);
    axes_do();
    many_pops_hi_lite(draw_win);
    return;
}

void
ggets_blank_screen(Window window) {
    curs_x = 0;
    curs_y = 0;
    xor_flag = 0;
    XClearWindow(display, window);
    return;
}

void
ggets_set_fore(void) {
    XSetForeground(display, gc, my_fore_color);
    return;
}

void
ggets_set_back(void) {
    XSetForeground(display, gc, my_back_color);
    return;
}

void
ggets_show_char(int32 ch, int32 col, int32 row, Window or) {
    char bob[2];
    bob[0] = (char)ch;
    ggets_chk_xor();
    XDrawString(display, or, gc, col, row + cury_off, bob, 1);
    return;
}

void
ggets_chk_xor(void) {
    if (xor_flag == 1) {
        XSetFunction(display, gc, GXxor);
    } else {
        XSetFunction(display, gc, GXcopy);
    }
    return;
}

void
ggets_clr_command(void) {
    ggets_blank_screen(command_pop);
    return;
}

void
ggets_draw_info_pop(Window window) {
    if (window == info_pop) {
        XClearWindow(display, info_pop);
        many_pops_base_col();
        XDrawString(display, info_pop, gc, 5, cury_off, info_message, (int)strlen(info_message));
    }
    return;
}

void
ggets_bottom_msg(char *msg) {
    XClearWindow(display, info_pop);
    many_pops_base_col();
    strcpy(info_message, msg);
    XDrawString(display, info_pop, gc, 5, cury_off, msg, (int)strlen(msg));
}

void
ggets_err_msg(char *string) {
    if (Xup) {
        pop_list_respond_box("OK", string);
    } else {
        ggets_plintf("%s\n", string);
    }
    return;
}

int32
ggets_plintf(char *fmt, ...) {
    int32 nchar = 0;
    va_list arglist;

    if (!xpp_verbose) {
        return nchar;  // Don't print at all!
    }

    if (logfile == NULL) {
        printf("The log file is NULL!\n");
        logfile = stdout;
    }

    va_start(arglist, fmt);
    nchar = vfprintf(logfile, fmt, arglist);
    va_end(arglist);
    /*Makes sense to flush to the output file to
    prevent loss of log info if program crashes.
    Then maybe user can figure out what happened and when.*/
    fflush(logfile);

    return nchar;
}

int32
ggets_show_position(XEvent event) {
    many_pops_check_draw_button(event);
    return 0;
}

void
ggets_put_command(char *string) {
    ggets_clr_command();
    ggets_f_text(0, 0, string, command_pop);
    curs_x = (int32)strlen(string);
    return;
}

int32
ggets_get_key_press(XEvent *event) {
    int32 maxlen = 64;
    char buf[65];
    XComposeStatus comp;
    KeySym ks;

    XLookupString((XKeyEvent *)event, buf, maxlen, &ks, &comp);

    if (ks == XK_Escape) {
        return KEY_ESC;
    }

    if ((ks == XK_Return) || (ks == XK_KP_Enter) || (ks == XK_Linefeed)) {
        return KEY_FINE;
    } else if (((ks >= XK_KP_Space) && (ks <= XK_KP_9)) ||
               ((ks >= XK_space) && (ks <= XK_asciitilde))) {
        return (int32)buf[0];
    }

    else if (ks == XK_BackSpace) {
        return KEY_BKSP;
    } else if (ks == XK_Delete) {
        return KEY_DEL;
    } else if (ks == XK_Tab) {
        return KEY_TAB;
    } else if (ks == XK_Home) {
        return KEY_HOME;
    } else if (ks == XK_End) {
        return KEY_END;
    } else if (ks == XK_Left) {
        return KEY_LEFT;
    } else if (ks == XK_Right) {
        return KEY_RIGHT;
    } else if (ks == XK_Up) {
        return KEY_UP;
    } else if (ks == XK_Down) {
        return KEY_DOWN;
    } else if (ks == XK_PgUp) {
        return KEY_PGUP;
    } else if (ks == XK_PgDn) {
        return KEY_PGDN;
    } else {
        return KEY_BADKEY;
    }
}

void
ggets_cput_text(void) {
    char string[256];
    char new[256];
    int32 x;
    int32 y;
    int32 size = 2;
    int32 font = 0;
    Window temp;
    temp = main_win;
    strcpy(string, "");
    if (ggets_new_string("Text: ", string) == 0) {
        return;
    }
    if (string[0] == '%') {
        graphics_fillin_text(&string[1], new);
        strcpy(string, new);
    }  // this makes it permanent

    ggets_new_int("Size 0-4 :", &size);
    if (size > 4) {
        size = 4;
    }
    if (size < 0) {
        size = 0;
    }
    pop_list_message_box(&temp, 0, scale_y - 5*dcur_y, "Place text with mouse");
    if (menudrive_get_mouse_xy(&x, &y)) {
        many_pops_gr_col();
        graphics_fillin_text(string, new);
        graphics_special_put_text_x11(x, y, new, size);
        many_pops_add_label(string, x, y, size, font);
        many_pops_base_col();
    }
    browser_wait_a_sec(CLICK_TIME);
    XDestroyWindow(display, temp);
    return;
}

int32
ggets_mouse_xy(int32 *x, int32 *y, Window window) {
    int32 no_but = 1;
    char ch;
    XEvent event;
    *x = 0;
    *y = 0;
    while (no_but) {
        XNextEvent(display, &event);
        switch (event.type) {
        case Expose:
            many_pops_do_expose(event);
            break;
        case KeyPress:
            ch = (char)ggets_get_key_press(&event);
            if (ch == ESCAPE) {
                return 0;
            }
            if (ch == KEY_FINE) {
                return -2;
            }
            if (ch == KEY_TAB) {
                return -3;
            }
            break;
        case ButtonPress:
            if (event.xbutton.window != window) {
                return 0;
            }
            no_but = 0;
            *x = event.xbutton.x;
            *y = event.xbutton.y;
            return 1;
        default:
            break;
        }
    }
    return 0;
}

void
ggets_f_text(int32 x, int32 y, char *string, Window o) {
    ggets_chk_xor();
    XDrawString(display, o, gc, x, y + cury_off, string, (int)strlen(string));
    return;
}

void
ggets_bar(int32 x, int32 y, int32 x2, int32 y2, Window window) {
    XFillRectangle(display, window, gc, x, y, (uint)(x2 - x), (uint)(y2 - y));
    return;
}

void
ggets_rectangle(int32 x, int32 y, int32 x2, int32 y2, Window window) {
    XDrawRectangle(display, window, gc, x, y, (uint)(x2 - x), (uint)(y2 - y));
    return;
}

void
ggets_circle(int32 x, int32 y, int32 radius, Window window) {
    XDrawArc(display, window, gc, x - radius, y - radius, (uint)(2*radius), (uint)(2*radius), 0,
             360*64);
    return;
}

void
ggets_xline(int32 x0, int32 y0, int32 x1, int32 y1, Window window) {
    XDrawLine(display, window, gc_graph, x0, y0, x1, y1);
    return;
}

int32
ggets_new_float(char *name, double *value) {
    int32 done2;
    int32 flag;
    double newz;
    char tvalue[200];
    snprintf(tvalue, sizeof(tvalue), "%.16g", *value);
    done2 = ggets_new_string(name, tvalue);
    if (done2 == 0 || strlen(tvalue) == 0) {
        return -1;
    }

    if (tvalue[0] == '%') {
        flag = calc_do_calc(&tvalue[1], &newz);
        if (flag != -1) {
            *value = newz;
        }
        return 0;
    }
    *value = atof(tvalue);

    return 0;
}

int32
ggets_new_int(char *name, int32 *value) {
    char svalue[200];
    snprintf(svalue, sizeof(svalue), "%d", *value);
    if (ggets_new_string(name, svalue) == 0 || strlen(svalue) == 0) {
        return -1;
    }
    *value = atoi(svalue);
    return 0;
}

void
ggets_display_command(char *name, char *value, int32 pos) {
    int32 l = (int32)strlen(name);
    int32 m = (int32)strlen(value);

    ggets_set_fore();
    ggets_bar(0, 0, l*dcur_x, dcur_y + 4, command_pop);
    ggets_set_back();
    XDrawString(display, command_pop, gc, 0, cury_off, name, l);
    ggets_set_fore();
    if (m > 0) {
        XDrawString(display, command_pop, gc, l*dcur_x, cury_off, value, m);
        ggets_put_cursor_at(command_pop, dcur_x*l, pos);
    }
    return;
}

void
ggets_clr_line_at(Window window, int32 col0, int32 pos, int32 n) {
    XClearArea(display, window, col0 + pos*dcur_x, 0, (uint)((n + 2)*dcur_x), 2*(uint)dcur_y,
               False);
    return;
}

void
ggets_put_cursor_at(Window window, int32 col0, int32 pos) {
    int32 x1 = col0 + pos*dcur_x;
    int32 x2 = x1 + 1;
    int32 y1 = dcur_y - 2, y2 = 2;
    XDrawLine(display, window, gc, x1, y1, x1, y2);
    XDrawLine(display, window, gc, x2, y1, x2, y2);
    return;
}

void
ggets_put_string_at(Window window, int32 col, char *s, int32 off) {
    int32 l = (int32)strlen(s) - off;

    XDrawString(display, window, gc, col, cury_off, s + off, l);
    return;
}

void
ggets_mov_mem(char *s1, char *s2, int32 len) {
    for (int32 i = len - 1; i >= 0; i--) {
        s1[i] = s2[i];
    }
    return;
}

void
ggets_mem_mov(char *s1, char *s2, int32 len) {
    for (int32 i = 0; i < len; i++) {
        s1[i] = s2[i];
    }
    return;
}

void
ggets_edit_window(Window window, int32 *pos, char *value, int32 *col, int32 *done2, int32 ch) {
    int32 col0 = *col - *pos*dcur_x;

    *done2 = 0;
    switch (ch) {
    case KEY_LEFT:
        if (*pos > 0) {
            *pos = *pos - 1;
            *col -= dcur_x;
        } else {
            ggets_ping();
        }
        break;
    case KEY_RIGHT:
        if (*pos < (int32)strlen(value)) {
            *pos = *pos + 1;
            *col += dcur_x;
        } else {
            ggets_ping();
        }
        break;
    case KEY_HOME: {
        *pos = 0;
        *col = col0;
    } break;
    case KEY_END: {
        *pos = (int32)strlen(value);
        *col = *pos*dcur_x + col0;
    } break;
    case KEY_BADKEY:
    case KEY_DOWN:
    case KEY_UP:
    case KEY_PGUP:
    case KEY_PGDN:
        return;  // junk key
    case KEY_ESC:
        *done2 = -1;  // quit without saving
        return;
    case KEY_FINE:
        if (ms_style == 0) {
            *done2 = 1;
        } else {
            *done2 = 2;
        }
        return;  // save this guy
    case KEY_BKSP:
        /*
        *pos=0;
        *col=col0;
        value[0]=0;
        ggets_clr_line_at(w,col0,0,80);
        break; */
    case KEY_DEL:
        if (*pos > 0) {
            ggets_mem_mov(&value[*pos - 1], &value[*pos], (int32)strlen(value) - *pos + 1);
            *pos = *pos - 1;
            *col -= dcur_x;
        } else {
            ggets_ping();
        }
        break;
    case KEY_TAB:
        if (ms_style == 0) {
            *done2 = 2;
        } else {
            *done2 = 1;
        }
        return;
    default:
        if ((ch >= ' ') && (ch <= '~')) {
            ggets_mov_mem(&value[*pos + 1], &value[*pos], (int32)strlen(value) - *pos + 1);
            value[*pos] = (char)ch;
            *pos = *pos + 1;
            *col += dcur_x;
        }
        break;
    }  // end key cases
       // Now redraw everything !!
    ggets_clr_line_at(window, col0, 0, (int32)strlen(value));
    ggets_put_string_at(window, col0, value, 0);
    ggets_put_cursor_at(window, col0, *pos);

    XFlush(display);
    return;
}

void
ggets_edit_command_string(XEvent event, char *name, char *value, int32 *done2, int32 *pos,
                          int32 *col) {
    char ch;
    switch (event.type) {
    case ConfigureNotify:
    case Expose:
    case MapNotify:
        many_pops_do_expose(event);
        if (event.xexpose.window == command_pop) {
            ggets_display_command(name, value, *pos);
        }
        break;
    case ButtonPress:
        if (event.xbutton.window == command_pop) {
            XSetInputFocus(display, command_pop, RevertToParent, CurrentTime);
        }
        break;
    case KeyPress:
        ch = (char)ggets_get_key_press(&event);
        ggets_edit_window(command_pop, pos, value, col, done2, ch);
        break;
    default:
        break;
    }
    return;
}

int32
ggets_new_string(char *name, char *value) {
    char old_value[80];
    int32 done2 = 0;
    int32 pos = (int32)strlen(value);
    int32 col = (pos + (int32)strlen(name))*dcur_x;

    XEvent event;
    strcpy(old_value, value);
    ggets_clr_command();
    ggets_display_command(name, value, pos);
    while (done2 == 0) {
        XNextEvent(display, &event);
        ggets_edit_command_string(event, name, value, &done2, &pos, &col);
    }
    ggets_clr_command();
    if (done2 == 1 || done2 == 2) {
        return done2;
    }
    strcpy(value, old_value);
    return 0;
}
