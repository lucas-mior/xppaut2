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
extern int32 XPPBatch;
extern int32 SCALEX;
extern int32 SCALEY;
int32 MSStyle = 0;
extern int32 Xup;

extern int32 tfBell;
extern int32 screen;
int32 CURS_X;
int32 CURS_Y;
extern int32 DCURY;
extern int32 DCURX;
extern int32 CURY_OFF;
extern Window win;
extern Window command_pop;
extern Window info_pop;
extern Window draw_win;
extern Window main_win;
extern GC gc;
extern GC gc_graph;
extern uint32 MyBackColor;
extern uint32 MyForeColor;
int32 xor_flag;
extern FILE *logfile;
extern int32 XPPVERBOSE;

static void put_string_at(Window window, int32 col, char *s, int32 off);
static void clr_line_at(Window window, int32 col0, int32 pos, int32 n);

void
ggets_ping(void) {
    if (tfBell && !XPPBatch) {
        /*
        XkbBell allows window managers to react
        to bell events using possibly user-specified
        accessiblity options (e.g. visual bell)
        */

        XkbBell(display, command_pop, 100, (Atom)NULL);
    }
    /* Call to XBell seems to be ignored by many window managers
     * where XkbBell is not. if(tfBell&&!XPPBatch) XBell(display,100); */
    return;
}

void
ggets_reset_graphics(void) {
    ggets_blank_screen(draw_win);
    axes2_do();
    hi_lite(draw_win);
    return;
}

void
ggets_blank_screen(Window window)

{
    CURS_X = 0;
    CURS_Y = 0;
    xor_flag = 0;
    XClearWindow(display, window);
    return;
}

void
ggets_set_fore(void) {
    XSetForeground(display, gc, MyForeColor);
    return;
}

void
ggets_set_back(void) {
    XSetForeground(display, gc, MyBackColor);
    return;
}

void
ggets_show_char(int32 ch, int32 col, int32 row, Window or) {
    char bob[2];
    bob[0] = (char)ch;
    ggets_chk_xor();
    XDrawString(display, or, gc, col, row + CURY_OFF, bob, 1);
    return;
}

void
ggets_chk_xor(void) {
    if (xor_flag == 1)
        XSetFunction(display, gc, GXxor);
    else
        XSetFunction(display, gc, GXcopy);
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
        XDrawString(display, info_pop, gc, 5, CURY_OFF, info_message,
                    (int)strlen(info_message));
    }
    return;
}

void
ggets_bottom_msg(char *msg) {
    XClearWindow(display, info_pop);
    many_pops_base_col();
    strcpy(info_message, msg);
    XDrawString(display, info_pop, gc, 5, CURY_OFF, msg, (int)strlen(msg));
}

void
ggets_err_msg(char *string) {
    if (Xup)
        respond_box("OK", string);
    else {
        ggets_plintf("%s\n", string);
    }
    return;
}

int32
ggets_plintf(char *fmt, ...) {
    int32 nchar = 0;
    va_list arglist;

    if (!XPPVERBOSE)
        return nchar; /*Don't print at all!*/

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
    check_draw_button(event);
    return 0;
}

void
ggets_put_command(char *string) {
    ggets_clr_command();
    ggets_f_text(0, 0, string, command_pop);
    CURS_X = (int32)strlen(string);
    return;
}

int32
ggets_get_key_press(XEvent *event) {
    int32 maxlen = 64;
    char buf[65];
    XComposeStatus comp;
    KeySym ks;

    XLookupString((XKeyEvent *)event, buf, maxlen, &ks, &comp);
    /*       printf(" ks=%d buf[0]=%d char=%c \n",ks,(int32)buf[0],buf[0]); */

    if (ks == XK_Escape)
        return KEY_ESC;
    if ((ks == XK_Return) || (ks == XK_KP_Enter) || (ks == XK_Linefeed))
        return KEY_FINE;
    else if (((ks >= XK_KP_Space) && (ks <= XK_KP_9)) ||
             ((ks >= XK_space) && (ks <= XK_asciitilde)))
        return (int32)buf[0];
    /*   else if ((ks>=XK_Shift_L)&&(ks<=XK_Hyper_R)) return 0;
       else if ((ks>=XK_F1)&&(ks<=XK_F35))  return 0; */

    else if (ks == XK_BackSpace)
        return KEY_BKSP;
    else if (ks == XK_Delete)
        return KEY_DEL;
    else if (ks == XK_Tab)
        return KEY_TAB;
    else if (ks == XK_Home)
        return KEY_HOME;
    else if (ks == XK_End)
        return KEY_END;
    else if (ks == XK_Left)
        return KEY_LEFT;
    else if (ks == XK_Right)
        return KEY_RIGHT;
    else if (ks == XK_Up)
        return KEY_UP;
    else if (ks == XK_Down)
        return KEY_DOWN;
    else if (ks == XK_PgUp)
        return KEY_PGUP;
    else if (ks == XK_PgDn)
        return KEY_PGDN;
    else {
        return KEY_BADKEY;
    }
}

void
ggets_cput_text(void) {
    char string[256], new[256];
    int32 x, y, size = 2, font = 0;
    Window temp;
    temp = main_win;
    strcpy(string, "");
    if (ggets_new_string("Text: ", string) == 0)
        return;
    if (string[0] == '%') {
        graphics_fillin_text(&string[1], new);
        strcpy(string, new);
    } /* this makes it permanent */

    ggets_new_int("Size 0-4 :", &size);
    /* ggets_new_int("Font  0-times/1-symbol :",&font); */
    if (size > 4)
        size = 4;
    if (size < 0)
        size = 0;
    message_box(&temp, 0, SCALEY - 5*DCURY, "Place text with mouse");
    if (menudrive_get_mouse_xy(&x, &y)) {
        many_pops_gr_col();
        /* fancy_put_text_x11(x,y,string,size,font); */
        graphics_fillin_text(string, new);
        special_put_text_x11(x, y, new, size);
        add_label(string, x, y, size, font);
        many_pops_base_col();
    }
    browse_wait_a_sec(ClickTime);
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
            do_expose(event);
            break;
        case KeyPress:
            ch = (char)ggets_get_key_press(&event);
            if (ch == ESCAPE)
                return 0;
            if (ch == KEY_FINE)
                return -2;
            if (ch == KEY_TAB)
                return -3;
            break;
        case ButtonPress:
            if (event.xbutton.window != window)
                return 0;
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
    XDrawString(display, o, gc, x, y + CURY_OFF, string, (int)strlen(string));
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
    XDrawArc(display, window, gc, x - radius, y - radius, (uint)(2*radius),
             (uint)(2*radius), 0, 360*64);
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
    if (done2 == 0 || strlen(tvalue) == 0)
        return -1;

    if (tvalue[0] == '%') {
        flag = do_calc(&tvalue[1], &newz);
        if (flag != -1)
            *value = newz;
        return 0;
    }
    *value = atof(tvalue);

    return 0;
}

int32
ggets_new_int(char *name, int32 *value) {
    char svalue[200];
    snprintf(svalue, sizeof(svalue), "%d", *value);
    if (ggets_new_string(name, svalue) == 0 || strlen(svalue) == 0)
        return -1;
    *value = atoi(svalue);
    return 0;
}

void
ggets_display_command(char *name, char *value, int32 pos) {
    int32 l = (int32)strlen(name);
    int32 m = (int32)strlen(value);

    ggets_set_fore();
    ggets_bar(0, 0, l*DCURX, DCURY + 4, command_pop);
    ggets_set_back();
    XDrawString(display, command_pop, gc, 0, CURY_OFF, name, l);
    ggets_set_fore();
    if (m > 0) {
        XDrawString(display, command_pop, gc, l*DCURX, CURY_OFF, value, m);
        /* ggets_show_char('_',DCURX*(l+m),0,command_pop); */
        ggets_put_cursor_at(command_pop, DCURX*l, pos);
    }
    return;
}

void
clr_line_at(Window window, int32 col0, int32 pos, int32 n) {
    XClearArea(display, window, col0 + pos*DCURX, 0, (uint)((n + 2)*DCURX),
               2*(uint)DCURY, False);
    return;
}

void
ggets_put_cursor_at(Window window, int32 col0, int32 pos) {
    int32 x1 = col0 + pos*DCURX;
    int32 x2 = x1 + 1;
    int32 y1 = DCURY - 2, y2 = 2;
    /* XDrawString(display,w,gc,col0+pos*DCURX-1,DCURY,"^",1);*/
    XDrawLine(display, window, gc, x1, y1, x1, y2);
    XDrawLine(display, window, gc, x2, y1, x2, y2);
    return;
}

void
put_string_at(Window window, int32 col, char *s, int32 off) {
    int32 l = (int32)strlen(s) - off;

    XDrawString(display, window, gc, col, CURY_OFF, s + off, l);
    return;
}

void
ggets_mov_mem(char *s1, char *s2, int32 len) {
    int32 i;
    for (i = len - 1; i >= 0; i--)
        s1[i] = s2[i];
    return;
}

void
ggets_mem_mov(char *s1, char *s2, int32 len) {
    int32 i;
    for (i = 0; i < len; i++)
        s1[i] = s2[i];
    return;
}

void
edit_window(Window window, int32 *pos, char *value, int32 *col, int32 *done2,
            int32 ch) {
    int32 col0 = *col - *pos*DCURX;

    *done2 = 0;
    switch (ch) {
    case KEY_LEFT:
        if (*pos > 0) {
            *pos = *pos - 1;
            *col -= DCURX;
        } else
            ggets_ping();
        break;
    case KEY_RIGHT:
        if (*pos < (int32)strlen(value)) {
            *pos = *pos + 1;
            *col += DCURX;
        } else
            ggets_ping();
        break;
    case KEY_HOME: {
        *pos = 0;
        *col = col0;
    } break;
    case KEY_END: {
        *pos = (int32)strlen(value);
        *col = *pos*DCURX + col0;
    } break;
    case KEY_BADKEY:
    case KEY_DOWN:
    case KEY_UP:
    case KEY_PGUP:
    case KEY_PGDN:
        return; /* junk key  */
    case KEY_ESC:
        *done2 = -1; /* quit without saving */
        return;
    case KEY_FINE:
        if (MSStyle == 0)
            *done2 = 1;
        else
            *done2 = 2;
        return; /* save this guy */
    case KEY_BKSP:
        /*
        *pos=0;
        *col=col0;
        value[0]=0;
        clr_line_at(w,col0,0,80);
        break; */
    case KEY_DEL:
        if (*pos > 0) {
            ggets_mem_mov(&value[*pos - 1], &value[*pos],
                   (int32)strlen(value) - *pos + 1);
            *pos = *pos - 1;
            *col -= DCURX;
        } else
            ggets_ping();
        break;
    case KEY_TAB:
        if (MSStyle == 0)
            *done2 = 2;
        else
            *done2 = 1;
        return;
    default:
        if ((ch >= ' ') && (ch <= '~')) {
            ggets_mov_mem(&value[*pos + 1], &value[*pos],
                   (int32)strlen(value) - *pos + 1);
            value[*pos] = (char)ch;
            *pos = *pos + 1;
            *col += DCURX;
        }
        break;
    } /* end key cases */
    /* Now redraw everything !!  */
    clr_line_at(window, col0, 0, (int32)strlen(value));
    put_string_at(window, col0, value, 0);
    ggets_put_cursor_at(window, col0, *pos);

    XFlush(display);
    return;
}

void
edit_command_string(XEvent event, char *name, char *value, int32 *done2,
                    int32 *pos, int32 *col) {
    char ch;
    switch (event.type) {
    case ConfigureNotify:
    case Expose:
    case MapNotify:
        do_expose(event);
        if (event.xexpose.window == command_pop)
            ggets_display_command(name, value, *pos);
        break;
    case ButtonPress:
        if (event.xbutton.window == command_pop)
            XSetInputFocus(display, command_pop, RevertToParent, CurrentTime);
        break;
    case KeyPress:
        ch = (char)ggets_get_key_press(&event);
        /* printf("ch= %ld \n",ch); */
        edit_window(command_pop, pos, value, col, done2, ch);
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
    int32 col = (pos + (int32)strlen(name))*DCURX;

    XEvent event;
    strcpy(old_value, value);
    ggets_clr_command();
    ggets_display_command(name, value, pos);
    while (done2 == 0) {
        XNextEvent(display, &event);
        edit_command_string(event, name, value, &done2, &pos, &col);
    }
    ggets_clr_command();
    if (done2 == 1 || done2 == 2)
        return done2;
    strcpy(value, old_value);
    return 0;
}
