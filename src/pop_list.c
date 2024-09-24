#include <stdlib.h>
#include <string.h>
#include "functions.h"
#include "integers.h"
#include <stdbool.h>

#include "info.bitmap"
#include "alert.bitmap"

typedef struct ScrollBox {
    Window base;
    Window slide;
    Window *w;
    int32 nw, nent, i0;
    int32 len;
    int32 exist;
    char **list;
} ScrollBox;

static void destroy_scroll_box(ScrollBox *sb);
static void create_scroll_box(Window root, int32 x0, int32 y0, int32 nent,
                              int32 nw, char **list, ScrollBox *sb);
static void expose_scroll_box(Window window, ScrollBox sb);
static void redraw_scroll_box(ScrollBox sb);
static void crossing_scroll_box(Window window, int32 c, ScrollBox sb);
static int32 scroll_box_motion(XEvent event, ScrollBox *sb);
static int32 select_scroll_item(Window window, ScrollBox sb);
static void scroll_popup(STRING_BOX *sb, ScrollBox *scrb);
static int32 s_box_event_loop(STRING_BOX *sb, int32 *pos, int32 *col,
                              ScrollBox *scrb);
static void draw_pop_up(POP_UP p, Window window);
static void set_sbox_item(STRING_BOX *sb, int32 item);
static void reset_hot(int32 inew, STRING_BOX *sb);
static void expose_sbox(STRING_BOX sb, Window window, int32 pos);
static int32 get_x_coord_win(Window win);

TEXTWIN mytext;
SCRBOX_LIST scrbox_list[10];
void
pop_list_set_window_title(Window win, char *string) {
    XTextProperty wname;
    XTextProperty iname;
    XStringListToTextProperty(&string, 1, &wname);
    XStringListToTextProperty(&string, 1, &iname);

    {
        XClassHint class_hints;
        class_hints.res_name = "";
        class_hints.res_class = "";

        XSetWMProperties(display, win, &wname, &iname, NULL, 0, NULL, NULL,
                         &class_hints);
    }
    return;
}

/* these are the standard lists that are possible */

void
pop_list_make_scrbox_lists(void) {
    int32 i;
    int32 n;
    static char *method[] = {
        "Discrete", "Euler",    "Mod. Euler", "Runge-Kutta", "Adams",
        "Gear",     "Volterra", "BackEul",    "QualRK",      "Stiff",
        "CVode",    "DoPri5",   "DoPri8(3)",  "Rosenbrock",  "Symplectic"};
    /* plottable list */
    scrbox_list[0].n = NEQ + 1;
    scrbox_list[0].list = xmalloc((usize)(NEQ + 1)*sizeof(char *));
    scrbox_list[0].list[0] = xmalloc(5);
    strcpy(scrbox_list[0].list[0], "T");
    for (i = 0; i < NEQ; i++) {
        scrbox_list[0].list[i + 1] = xmalloc(15);
        strcpy(scrbox_list[0].list[i + 1], uvar_names[i]);
    }
    /* variable list */
    scrbox_list[1].n = NODE + NMarkov;
    scrbox_list[1].list = xmalloc((usize)(NODE + NMarkov)*sizeof(char *));
    for (i = 0; i < NODE + NMarkov; i++) {
        scrbox_list[1].list[i] = xmalloc(15);
        strcpy(scrbox_list[1].list[i], uvar_names[i]);
    }

    /* parameter list */
    scrbox_list[2].n = NUPAR;
    scrbox_list[2].list = xmalloc((usize)NUPAR*sizeof(char *));
    for (i = 0; i < NUPAR; i++) {
        scrbox_list[2].list[i] = xmalloc(15);
        strcpy(scrbox_list[2].list[i], upar_names[i]);
    }

    /* parvar list */
    n = NODE + NMarkov + NUPAR;
    scrbox_list[3].n = n;
    scrbox_list[3].list = xmalloc((usize)n*sizeof(char *));
    for (i = 0; i < NODE + NMarkov; i++) {
        scrbox_list[3].list[i] = xmalloc(15);
        strcpy(scrbox_list[3].list[i], uvar_names[i]);
    }
    for (i = NODE + NMarkov; i < n; i++) {
        scrbox_list[3].list[i] = xmalloc(15);
        strcpy(scrbox_list[3].list[i], upar_names[i - NODE - NMarkov]);
    }
    /* color list */
    scrbox_list[4].n = 11;
    scrbox_list[4].list = xmalloc(11*sizeof(char *));
    for (i = 0; i < 11; i++) {
        scrbox_list[4].list[i] = xmalloc(20);
        sprintf(scrbox_list[4].list[i], "%d %s", i, color_names[i]);
    }
    /* marker list */
    scrbox_list[5].n = 6;
    scrbox_list[5].list = xmalloc(6*sizeof(char *));
    for (i = 0; i < 6; i++)
        scrbox_list[5].list[i] =
            xmalloc(13*sizeof(*(scrbox_list[5].list[i])));
    strcpy(scrbox_list[5].list[0], "2 Box");
    strcpy(scrbox_list[5].list[1], "3 Diamond");
    strcpy(scrbox_list[5].list[2], "4 Triangle");
    strcpy(scrbox_list[5].list[3], "5 Plus");
    strcpy(scrbox_list[5].list[4], "6 X");
    strcpy(scrbox_list[5].list[5], "7 circle2");
    /* method list */
    scrbox_list[6].list = xmalloc(15*sizeof(char *));
    scrbox_list[6].n = 15;
    for (i = 0; i < 15; i++) {
        scrbox_list[6].list[i] =
            xmalloc(22*sizeof(*(scrbox_list[6].list[i])));
        sprintf(scrbox_list[6].list[i], "%d %s", i, method[i]);
    }
    return;
}

int32
get_x_coord_win(Window win) {
    int32 x;
    int32 y;
    uint32 h, w, bw, d;
    Window root;
    XGetGeometry(display, win, &root, &x, &y, &w, &h, &bw, &d);
    return x;
}

void
destroy_scroll_box(ScrollBox *sb) {
    if (sb->exist == 1) {
        sb->exist = 0;
        browse_wait_a_sec(ClickTime);
        XDestroySubwindows(display, sb->base);
        XDestroyWindow(display, sb->base);
    }
    return;
}

void
create_scroll_box(Window root, int32 x0, int32 y0, int32 nent, int32 nw,
                  char **list, ScrollBox *sb) {
    int32 slen = 0;
    int32 hgt;
    int32 wid;
    int32 ww, len;
    int32 hw = DCURYs + 4;
    for (int32 i = 0; i < nent; i++) {
        if (slen < (int32)strlen(list[i]))
            slen = (int32)strlen(list[i]);
    }
    wid = (slen + 2)*(DCURXs);
    ww = slen*DCURXs + DCURXs / 2;
    hgt = hw*(nw + 1);
    len = hgt - 6;
    sb->base = (Window)make_plain_window(root, x0, y0, wid, hgt, 2);
    sb->w = xmalloc((usize)nw*sizeof(*(sb->w)));
    for (int32 i = 0; i < nw; i++)
        sb->w[i] = make_window(sb->base, 1, hw / 2 + i*hw, ww, DCURYs, 0);
    sb->i0 = 0;
    sb->nw = nw;
    sb->nent = nent;
    sb->list = list;
    if (sb->nw < sb->nent)
        sb->slide = make_window(sb->base, ww + DCURXs / 2 + 2, 2,
                                ww + DCURXs / 2 + 6, 2 + len, 1);
    sb->len = len - 4;
    sb->exist = 1;
    return;
}

void
expose_scroll_box(Window window, ScrollBox sb) {
    int32 i;
    /*int32 flag=-1;*/
    for (i = 0; i < sb.nw; i++)
        if (window == sb.w[i]) {
            redraw_scroll_box(sb);
            return;
        }
    if (sb.nw < sb.nent && window == sb.slide)
        redraw_scroll_box(sb);
    return;
}

void
redraw_scroll_box(ScrollBox sb) {
    int32 i;
    int32 p;
    int32 i0 = sb.i0;
    for (i = 0; i < sb.nw; i++) {
        XClearWindow(display, sb.w[i]);
        XDrawString(display, sb.w[i], small_gc, 0, CURY_OFFs, sb.list[i + i0],
                    (int)strlen(sb.list[i + i0]));
    }
    if (sb.nw < sb.nent) {
        XClearWindow(display, sb.slide);
        /* now calculate the slide position */
        p = 2 + (sb.i0*sb.len) / (sb.nent - sb.nw);
        for (i = -2; i <= 2; i++)
            XDrawLine(display, sb.slide, small_gc, 0, p + i, 5, p + i);
    }
    return;
}

void
crossing_scroll_box(Window window, int32 c, ScrollBox sb) {
    int32 i;
    for (i = 0; i < sb.nw; i++) {
        if (window == sb.w[i]) {
            XSetWindowBorderWidth(display, window, (uint)c);
            return;
        }
    }
    return;
}

int32
scroll_box_motion(XEvent event, ScrollBox *sb) {
    int32 x;
    Window window;
    int32 pos;
    int32 len;
    window = event.xmotion.window;
    x = event.xmotion.y;
    if (sb->nw >= sb->nent)
        return 0;
    if (window == sb->slide) {
        len = sb->len;
        if (x < 2)
            x = 2;
        if (x > (len + 2))
            x = len + 2;
        pos = ((x - 2)*(sb->nent - sb->nw)) / len;
        if (pos < 0)
            pos = 0;
        if (pos > (sb->nent - sb->nw))
            pos = sb->nent - sb->nw;
        sb->i0 = pos;
        redraw_scroll_box(*sb);
    }
    return 1;
}

int32
select_scroll_item(Window window, ScrollBox sb) {
    int32 i;
    int32 item = -1;
    for (i = 0; i < sb.nw; i++) {
        if (window == sb.w[i]) {
            item = i + sb.i0;
            return item;
        }
    }
    return -1;
}

void
scroll_popup(STRING_BOX *sb, ScrollBox *scrb) {
    int32 hw = DCURYs + 4;
    int32 ihot = sb->hot;
    int32 id = sb->hh[ihot];
    int32 xx;
    int32 maxhgt = sb->hgt;
    int32 maxw;
    if (id < 0)
        return; /* shouldnt happen */
    maxw = maxhgt / hw - 1;
    if (maxw > scrbox_list[id].n)
        maxw = scrbox_list[id].n;
    xx = get_x_coord_win(sb->win[ihot]);
    create_scroll_box(sb->base, xx, 3, scrbox_list[id].n, maxw,
                      scrbox_list[id].list, scrb);
    return;
}

int32
do_string_box(int32 n, int32 row, int32 col, char *title, char **names,
              char values[][MAX_LEN_SBOX], int32 maxchar) {
    STRING_BOX sb;
    int32 i;
    int32 status;
    int32 colm, pos;
    ScrollBox scrb;
    scrb.exist = 0;

    for (i = 0; i < n; i++) {
        sb.hh[i] = -1;
        if (names[i][0] == '*') {
            sb.hh[i] = atoi(names[i] + 1);
            snprintf(sb.name[i], sizeof(sb.name[i]), "*%s:", names[i] + 2);
        } else
            snprintf(sb.name[i], sizeof(sb.name[i]), "%s:", names[i]);
        strcpy(sb.value[i], values[i]);
    }
    sb.n = n;
    sb.hot = 0;
    make_sbox_windows(&sb, row, col, title, maxchar);
    XSelectInput(display, sb.cancel, BUT_MASK);
    XSelectInput(display, sb.ok, BUT_MASK);
    pos = (int32)strlen(sb.value[0]);
    colm = (pos + (int32)strlen(sb.name[0]))*DCURX;

    while (true) {
        status = s_box_event_loop(&sb, &pos, &colm, &scrb);
        if (status != -1)
            break;
    }
    XSelectInput(display, sb.cancel, EV_MASK);
    XSelectInput(display, sb.ok, EV_MASK);

    browse_wait_a_sec(ClickTime);
    XDestroySubwindows(display, sb.base);
    XDestroyWindow(display, sb.base);

    if (status == FORGET_ALL)
        return status;
    for (i = 0; i < n; i++)
        strcpy(values[i], sb.value[i]);
    return status;
}

void
expose_sbox(STRING_BOX sb, Window window, int32 pos) {
    int32 i;
    int32 flag;

    if (window == sb.ok) {
        XDrawString(display, window, gc, 5, CURY_OFF, "Ok", 2);
        return;
    }
    if (window == sb.cancel) {
        XDrawString(display, window, gc, 5, CURY_OFF, "Cancel", 6);
        return;
    }
    for (i = 0; i < sb.n; i++) {
        if (window != sb.win[i])
            continue;
        flag = 0;
        if (i == sb.hot)
            flag = 1;
        do_hilite_text(sb.name[i], sb.value[i], flag, window, pos);
    }
    return;
}

void
do_hilite_text(char *name, char *value, int32 flag, Window window, int32 pos) {
    int32 l = (int32)strlen(name);
    int32 m = (int32)strlen(value);
    if (flag) {
        ggets_set_fore();
        ggets_bar(0, 0, l*DCURX, DCURY + 4, window);
        ggets_set_back();
    }
    XDrawString(display, window, gc, 0, CURY_OFF, name, l);
    ggets_set_fore();
    if (m > 0)
        XDrawString(display, window, gc, l*DCURX, CURY_OFF, value, m);
    /* if(flag) ggets_show_char('_',DCURX*(l+m),0,w); */
    if (flag)
        ggets_put_cursor_at(window, DCURX*l, pos);
    return;
}

void
reset_hot(int32 inew, STRING_BOX *sb) {
    int32 i = sb->hot;
    sb->hot = inew;
    XClearWindow(display, sb->win[inew]);
    do_hilite_text(sb->name[inew], sb->value[inew], 1, sb->win[inew],
                   (int)strlen(sb->value[inew]));
    XClearWindow(display, sb->win[i]);
    do_hilite_text(sb->name[i], sb->value[i], 0, sb->win[i],
                   (int)strlen(sb->value[i]));
    return;
}

void
new_editable(STRING_BOX *sb, int32 inew, int32 *pos, int32 *col, int32 *done,
             Window *w) {
    reset_hot(inew, sb);
    *pos = (int32)strlen(sb->value[inew]);
    *col = (*pos + (int32)strlen(sb->name[inew]))*DCURX;
    *done = 0;
    *w = sb->win[inew];
    return;
}

void
set_sbox_item(STRING_BOX *sb, int32 item) {
    int32 i = sb->hot;
    int32 id = sb->hh[i];
    if (id < 0)
        return;
    strcpy(sb->value[i], scrbox_list[id].list[item]);
    return;
}

int32
s_box_event_loop(STRING_BOX *sb, int32 *pos, int32 *col, ScrollBox *scrb) {
    XEvent event;
    int32 status = -1, inew;
    int32 nn = sb->n;
    int32 done = 0, j;
    int32 item;
    char ch;
    int32 ihot = sb->hot;
    Window wt;
    Window window = sb->win[ihot]; /* active window   */
    char *s;
    s = sb->value[ihot];

    XNextEvent(display, &event);
    switch (event.type) {
    case ConfigureNotify:
    case Expose:
    case MapNotify:
        many_pops_do_expose(event); /*  menus and graphs etc  */
        expose_sbox(*sb, event.xany.window, *pos);
        if (scrb->exist)
            expose_scroll_box(event.xany.window, *scrb);
        break;
    case MotionNotify:
        if (scrb->exist)
            scroll_box_motion(event, scrb);
        break;
    case ButtonPress:
        if (scrb->exist) {
            item = select_scroll_item(event.xbutton.window, *scrb);
            if (item >= 0) {
                set_sbox_item(sb, item);
                new_editable(sb, sb->hot, pos, col, &done, &window);
                destroy_scroll_box(scrb);
            }
        }
        if (event.xbutton.window == sb->ok) {
            destroy_scroll_box(scrb);
            status = DONE_ALL;
            break;
        }
        if (event.xbutton.window == sb->cancel) {
            status = FORGET_ALL;
            break;
        }
        for (int32 i = 0; i < nn; i++) {
            if (event.xbutton.window == sb->win[i]) {
                XSetInputFocus(display, sb->win[i], RevertToParent,
                               CurrentTime);
                if (i != sb->hot) {
                    destroy_scroll_box(scrb);
                    new_editable(sb, i, pos, col, &done, &window);
                } else { /* i==sb->hot */
                    if (event.xbutton.x < DCURX) {
                        j = sb->hot;
                        if (sb->hh[j] >= 0) {
                            scroll_popup(sb, scrb);
                        }
                    }
                }

                break;
            }
        }
        break;

    case EnterNotify:
        wt = event.xcrossing.window;
        if (scrb->exist)
            crossing_scroll_box(wt, 1, *scrb);
        if (wt == sb->ok || wt == sb->cancel)
            XSetWindowBorderWidth(display, wt, 2);
        break;

    case LeaveNotify:
        wt = event.xcrossing.window;
        if (scrb->exist)
            crossing_scroll_box(wt, 0, *scrb);
        if (wt == sb->ok || wt == sb->cancel)
            XSetWindowBorderWidth(display, wt, 1);
        break;

    case KeyPress:
        ch = (char)ggets_get_key_press(&event);
        edit_window(window, pos, s, col, &done, ch);
        if (done != 0) {
            if (done == DONE_ALL) {
                status = DONE_ALL;
                break;
            }
            inew = (sb->hot + 1) % nn;
            new_editable(sb, inew, pos, col, &done, &window);
        }
        break;
    default:
        break;
    }
    return status;
}

void
make_sbox_windows(STRING_BOX *sb, int32 row, int32 col, char *title,
                  int32 maxchar) {
    int32 width;
    int32 height;
    int32 i;
    int32 xpos, ypos, n = sb->n;
    int32 xstart;
    int32 ystart;

    XTextProperty winname;
    XSizeHints size_hints;
    Window base;
    width = (maxchar + 4)*col*DCURX;
    height = (row + 4)*(DCURY + 16);
    base =
        make_plain_window(DefaultRootWindow(display), 0, 0, width, height, 4);
    XStringListToTextProperty(&title, 1, &winname);
    size_hints.flags = PPosition | PSize | PMinSize | PMaxSize;
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

        many_pops_make_icon((char *)info_bits, info_width, info_height, base);
        XSetWMProperties(display, base, &winname, NULL, NULL, 0, &size_hints,
                         NULL, &class_hints);
    }
    sb->base = base;
    sb->hgt = height;
    sb->wid = width;
    ystart = DCURY;
    xstart = DCURX;
    for (i = 0; i < n; i++) {
        xpos = xstart + (maxchar + 4)*DCURX*(i / row);
        ypos = ystart + (i % row)*(DCURY + 10);
        sb->win[i] = make_window(base, xpos, ypos, maxchar*DCURX, DCURY, 1);
    }

    ypos = height - 2*DCURY;
    xpos = (width - 16*DCURX) / 2;
    (sb->ok) = make_window(base, xpos, ypos, 8*DCURX, DCURY, 1);
    (sb->cancel) =
        make_window(base, xpos + 8*DCURX + 4, ypos, 8*DCURX, DCURY, 1);
    XRaiseWindow(display, base);
    return;
}

Window
make_fancy_window(Window root, int32 x, int32 y, int32 width, int32 height,
                  int32 bw) {
    Window win;
    win = XCreateSimpleWindow(display, root, x, y, (uint)width, (uint)height,
                              (uint)bw, MyForeColor, MyBackColor);

    if (UserGradients == 1) {
        int32 xx, yy;
        double cosine;
        XColor bcolour, col2, diffcol;
        Colormap cmap = DefaultColormap(display, DefaultScreen(display));
        Pixmap pmap =
            XCreatePixmap(display, root, (uint)width, (uint)height,
                          (uint)DefaultDepth(display, DefaultScreen(display)));

        /*double l2rads;*/
        xx = 0;

        XParseColor(display, cmap, UserWhite, &bcolour);
        XParseColor(display, cmap, UserBlack, &diffcol);

        /*l2rads = 3.1415926535897932384/(1.0*height);
         */

        /*win=XCreateSimpleWindow(display,root,x,y,width,height,
               bw,diffcol.pixel,bcolour.pixel);
        */

        for (yy = 0; yy < height; yy += 1) {
            if (yy < 1.0) {
                col2.red = 65535;
                col2.green = 65355;
                col2.blue = 65355;
            } else {
                if (yy < (height / 2.0)) {
                    cosine = 1.0;
                } else if ((height - yy) <= 1.0) {
                    cosine = 0.1;
                } else {
                    cosine = 0.93;
                }
                col2.red = (ushort)(bcolour.red*cosine);
                col2.green = (ushort)(bcolour.green*cosine);
                col2.blue = (ushort)(bcolour.blue*cosine);
            }

            XAllocColor(display, cmap, &col2);
            XSetForeground(display, gc, col2.pixel);

            for (xx = 1; xx < width - 1; xx += 1) {
                XDrawPoint(display, pmap, gc, xx, yy);
            }

            /*Now do xx=0 and xx=width-1*/
            xx = 0;
            col2.red = 65535;
            col2.green = 65355;
            col2.blue = 65355;
            XAllocColor(display, cmap, &col2);
            XSetForeground(display, gc, col2.pixel);
            XDrawPoint(display, pmap, gc, xx, yy);
            xx = width - 1;
            cosine = 0.1;
            col2.red = (ushort)(bcolour.red*cosine);
            col2.green = (ushort)(bcolour.green*cosine);
            col2.blue = (ushort)(bcolour.blue*cosine);

            XAllocColor(display, cmap, &col2);
            XSetForeground(display, gc, col2.pixel);
            XDrawPoint(display, pmap, gc, xx, yy);
        }

        XSetWindowBackgroundPixmap(display, win, pmap);
        XFreePixmap(display, pmap);
    }

    XSelectInput(display, win,
                 ExposureMask | KeyPressMask | ButtonPressMask |
                     StructureNotifyMask | ButtonReleaseMask |
                     ButtonMotionMask | LeaveWindowMask | EnterWindowMask);
    XMapWindow(display, win);

    return win;
}

Window
make_unmapped_window(Window root, int32 x, int32 y, int32 width, int32 height,
                     int32 bw) {
    Window win;
    win = XCreateSimpleWindow(display, root, x, y, (uint)width, (uint)height,
                              (uint)bw, MyForeColor, MyBackColor);

    /*Gradient stuff*/

    if (UserGradients == 1) {
        int32 xx, yy;
        double cosine;
        Pixmap pmap =
            XCreatePixmap(display, root, (uint)width, (uint)height,
                          (uint)DefaultDepth(display, DefaultScreen(display)));
        XColor bcolour, col2, diffcol;
        Colormap cmap = DefaultColormap(display, DefaultScreen(display));

        /*double l2rads;
         */
        xx = 0;

        XParseColor(display, cmap, UserWhite, &bcolour);
        XParseColor(display, cmap, UserBlack, &diffcol);

        /*l2rads = 3.1415926535897932384/(1.0*height);*/

        /* win=XCreateSimpleWindow(display,root,x,y,width,height,
                bw,diffcol.pixel,bcolour.pixel);
         */

        for (yy = 0; yy < height; yy += 1) {
            if (yy < 1.0) {
                col2.red = 65535;
                col2.green = 65355;
                col2.blue = 65355;
            } else {
                if (yy < (height / 2.0)) {
                    cosine = 1.0;
                } else if ((height - yy) <= 1.0) {
                    cosine = 0.1;
                } else {
                    cosine = 0.93;
                }
                col2.red = (ushort)(bcolour.red*cosine);
                col2.green = (ushort)(bcolour.green*cosine);
                col2.blue = (ushort)(bcolour.blue*cosine);
            }

            XAllocColor(display, cmap, &col2);
            XSetForeground(display, gc, col2.pixel);

            for (xx = 1; xx < width - 1; xx += 1) {
                XDrawPoint(display, pmap, gc, xx, yy);
            }

            /*Now do xx=0 and xx=width-1*/
            xx = 0;
            col2.red = 65535;
            col2.green = 65355;
            col2.blue = 65355;
            XAllocColor(display, cmap, &col2);
            XSetForeground(display, gc, col2.pixel);
            XDrawPoint(display, pmap, gc, xx, yy);
            xx = width - 1;
            cosine = 0.1;
            col2.red = (ushort)(bcolour.red*cosine);
            col2.green = (ushort)(bcolour.green*cosine);
            col2.blue = (ushort)(bcolour.blue*cosine);

            XAllocColor(display, cmap, &col2);
            XSetForeground(display, gc, col2.pixel);
            XDrawPoint(display, pmap, gc, xx, yy);
        }
        XSetWindowBackgroundPixmap(display, win, pmap);
        XFreePixmap(display, pmap);
    }

    XSelectInput(display, win,
                 ExposureMask | KeyPressMask | ButtonPressMask |
                     StructureNotifyMask | ButtonReleaseMask |
                     ButtonMotionMask | LeaveWindowMask | EnterWindowMask);

    return win;
}

static void
bin_prnt_byte(int32 x, int32 *arr) {
    int32 n = 0;
    for (n = 7; n >= 0; n--) {
        if ((x & 0x80) != 0) {
            arr[n] = 1;
        } else {
            arr[n] = 0;
        }

        x = x << 1;
    }

    return;
}

Window
make_unmapped_icon_window(Window root, int32 x, int32 y, int32 width,
                          int32 height, int32 bw, uchar *icdata) {
    Window win =
        XCreateSimpleWindow(display, root, x, y, (uint)width, (uint)height,
                            (uint)bw, MyForeColor, MyBackColor);

    /*Gradient stuff*/

    Pixmap pmap =
        XCreatePixmap(display, root, (uint)width, (uint)height,
                      (uint)DefaultDepth(display, DefaultScreen(display)));
    int32 xx;
    int32 yy;
    int32 row = 0, col = 0;

    XColor bcolour, col2, diffcol;
    Colormap cmap = DefaultColormap(display, DefaultScreen(display));
    XParseColor(display, cmap, UserWhite, &bcolour);
    XParseColor(display, cmap, UserBlack, &diffcol);

    if (UserGradients == 1) {
        double cosine;

        /*double l2rads;
         */
        xx = 0;

        /*l2rads = 3.1415926535897932384/(1.0*height);*/

        /* win=XCreateSimpleWindow(display,root,x,y,width,height,
                bw,diffcol.pixel,bcolour.pixel);
         */

        for (yy = 0; yy < height; yy += 1) {
            if (yy < 1.0) {
                col2.red = 65535;
                col2.green = 65355;
                col2.blue = 65355;
            } else {
                if (yy < (height / 2.0)) {
                    cosine = 1.0;
                } else if ((height - yy) <= 1.0) {
                    cosine = 0.1;
                } else {
                    cosine = 0.93;
                }
                col2.red = (ushort)(bcolour.red*cosine);
                col2.green = (ushort)(bcolour.green*cosine);
                col2.blue = (ushort)(bcolour.blue*cosine);
            }

            XAllocColor(display, cmap, &col2);
            XSetForeground(display, gc, col2.pixel);

            for (xx = 1; xx < width - 1; xx += 1) {
                XDrawPoint(display, pmap, gc, xx, yy);
            }

            /*Now do xx=0 and xx=width-1*/
            xx = 0;
            col2.red = 65535;
            col2.green = 65355;
            col2.blue = 65355;
            XAllocColor(display, cmap, &col2);
            XSetForeground(display, gc, col2.pixel);
            XDrawPoint(display, pmap, gc, xx, yy);
            xx = width - 1;
            cosine = 0.1;
            col2.red = (ushort)(bcolour.red*cosine);
            col2.green = (ushort)(bcolour.green*cosine);
            col2.blue = (ushort)(bcolour.blue*cosine);

            XAllocColor(display, cmap, &col2);
            XSetForeground(display, gc, col2.pixel);
            XDrawPoint(display, pmap, gc, xx, yy);
        }
    } else {
        col2.red = bcolour.red;
        col2.green = bcolour.green;
        col2.blue = bcolour.blue;
        XAllocColor(display, cmap, &col2);
        XSetForeground(display, gc, col2.pixel);

        for (yy = 0; yy < height; yy += 1) {
            for (xx = 0; xx < width; xx += 1) {
                XDrawPoint(display, pmap, gc, xx, yy);
            }
        }
    }

    if (icdata == NULL) {
        /*Don't do anything...*/

    } else {
        uchar *ps = icdata;
        int32 intstack[8];
        col2.red = diffcol.red;
        col2.green = diffcol.green;
        col2.blue = diffcol.blue;
        XAllocColor(display, cmap, &col2);
        XSetForeground(display, gc, col2.pixel);

        col = 0;
        row = -1;
        while (row < height) {
            col = 0;
            row++;
            while (true) {
                int32 q = 0;
                bin_prnt_byte(*ps, intstack);
                ps++;

                for (q = 0; q < 8; q++) /*8 bits per byte*/
                {
                    if (col >= width) {
                        /*myint[z] = 0;*/
                    } else {
                        if (intstack[q] == 1) {
                            XDrawPoint(display, pmap, gc, col, row);
                        }
                    }
                    col++;
                }

                if (col >= width) {
                    break;
                }
            }
        }
    }

    XSetWindowBackgroundPixmap(display, win, pmap);
    XFreePixmap(display, pmap);

    XSelectInput(display, win,
                 ExposureMask | KeyPressMask | ButtonPressMask |
                     StructureNotifyMask | ButtonReleaseMask |
                     ButtonMotionMask | LeaveWindowMask | EnterWindowMask);

    return win;
}

Window
make_plain_unmapped_window(Window root, int32 x, int32 y, int32 width,
                           int32 height, int32 bw) {
    Window win;
    win = XCreateSimpleWindow(display, root, x, y, (uint)width, (uint)height,
                              (uint)bw, MyForeColor, MyBackColor);

    XSelectInput(display, win,
                 ExposureMask | KeyPressMask | ButtonPressMask |
                     StructureNotifyMask | ButtonReleaseMask |
                     ButtonMotionMask | LeaveWindowMask | EnterWindowMask);

    return win;
}

Window
make_icon_window(Window root, int32 x, int32 y, int32 width, int32 height,
                 int32 bw, uchar *icdata) {
    Window win;
    win = make_unmapped_icon_window(root, x, y, width, height, bw, icdata);
    if (root == RootWindow(display, screen))
        XSetWMProtocols(display, win, &deleteWindowAtom, 1);
    XMapWindow(display, win);
    return win;
}

Window
make_window(Window root, int32 x, int32 y, int32 width, int32 height,
            int32 bw) {
    Window win;
    win = make_unmapped_window(root, x, y, width, height, bw);
    if (root == RootWindow(display, screen))
        XSetWMProtocols(display, win, &deleteWindowAtom, 1);
    XMapWindow(display, win);
    return win;
}

Window
make_plain_window(Window root, int32 x, int32 y, int32 width, int32 height,
                  int32 bw) {
    Window win;
    win = make_plain_unmapped_window(root, x, y, width, height, bw);
    if (root == RootWindow(display, screen))
        XSetWMProtocols(display, win, &deleteWindowAtom, 1);
    XMapWindow(display, win);
    return win;
}

void
expose_resp_box(char *button, char *message, Window wb, Window wm,
                Window window) {
    if (window == wb)
        ggets_f_text(0, 0, button, wb);
    if (window == wm)
        ggets_f_text(0, 0, message, wm);
    return;
}

void
pop_list_respond_box(char *button, char *message) {
    int32 l1 = (int32)strlen(message);
    int32 l2 = (int32)strlen(button);
    int32 width;
    int32 height;
    int32 done = 0;
    XEvent event;
    Window wmain, wb, wm;
    width = l1;
    if (l1 < l2)
        width = l2;
    width = width + 4;
    height = 5*DCURY;
    wmain = make_plain_window(RootWindow(display, screen), DisplayWidth / 2,
                              DisplayHeight / 2, width*DCURX, height, 4);
    many_pops_make_icon((char *)alert_bits, alert_width, alert_height, wmain);
    wm = make_plain_window(wmain, ((width - l1)*DCURX) / 2, DCURY / 2,
                           l1*DCURX, DCURY, 0);
    wb = make_window(wmain, ((width - l2)*DCURX) / 2, 2*DCURY, l2*DCURX,
                     DCURY, 1);

    ggets_ping();
    pop_list_set_window_title(wmain, "!!");
    XSelectInput(display, wb, BUT_MASK);
    while (!done) {
        XNextEvent(display, &event);
        switch (event.type) {
        case Expose:
        case MapNotify:
            many_pops_do_expose(event);
            expose_resp_box(button, message, wb, wm, event.xexpose.window);
            break;
        case KeyPress:
            done = 1;
            break;
        case ButtonPress:
            if (event.xbutton.window == wb)
                done = 1;
            break;
        case EnterNotify:
            if (event.xcrossing.window == wb)
                XSetWindowBorderWidth(display, wb, 2);
            break;
        case LeaveNotify:
            if (event.xcrossing.window == wb)
                XSetWindowBorderWidth(display, wb, 1);
            break;
        default:
            break;
        }
    }

    XSelectInput(display, wb, EV_MASK);
    browse_wait_a_sec(ClickTime);
    XDestroySubwindows(display, wmain);
    XDestroyWindow(display, wmain);
    return;
}

void
pop_list_message_box(Window *w, int32 x, int32 y, char *message) {
    int32 wid = (int32)strlen(message)*DCURX;
    int32 hgt = 4*DCURY;
    Window z;
    z = make_plain_window(*w, x, y, wid + 50, hgt, 4);
    XSelectInput(display, z, 0);
    ggets_f_text(25, 2*DCURY, message, z);
    ggets_ping();
    *w = z;
    return;
}

void
expose_choice(char *choice1, char *choice2, char *msg, Window c1, Window c2,
              Window wm, Window window) {
    if (window == wm)
        ggets_f_text(0, 0, msg, wm);
    if (window == c1)
        ggets_f_text(0, 0, choice1, c1);
    if (window == c2)
        ggets_f_text(0, 0, choice2, c2);
    return;
}

int32
two_choice(char *choice1, char *choice2, char *string, char *key, int32 x,
           int32 y, Window window, char *title) {
    Window base, c1, c2, wm;
    XEvent event;
    int32 not_done = 1;
    int32 value = 0;
    int32 l1 = (int32)strlen(choice1)*DCURX;
    int32 l2 = (int32)strlen(choice2)*DCURX;
    int32 lm = (int32)strlen(string)*DCURX;
    int32 tot = lm, xm, x1, x2;

    if (lm < (l1 + l2 + 4*DCURX))
        tot = (l1 + l2 + 4*DCURX);
    tot = tot + 6*DCURX;
    xm = (tot - lm) / 2;
    x1 = (tot - l1 - l2 - 4*DCURX) / 2;
    x2 = x1 + l1 + 4*DCURX;
    base = make_plain_window(window, x, y, tot, 5*DCURY, 4);

    many_pops_make_icon((char *)alert_bits, alert_width, alert_height, base);

    c1 = make_window(base, x1, 3*DCURY, l1 + DCURX, DCURY + 4, 1);
    c2 = make_window(base, x2, 3*DCURY, l2 + DCURX, DCURY + 4, 1);
    XSelectInput(display, c1, BUT_MASK);
    XSelectInput(display, c2, BUT_MASK);

    wm = make_window(base, xm, DCURY / 2, lm + 2, DCURY, 0);

    ggets_ping();
    if (window == RootWindow(display, screen)) {
        if (title == NULL) {
            pop_list_set_window_title(base, "!!!!");
        } else {
            pop_list_set_window_title(base, title);
        }
    }

    while (not_done) {
        XNextEvent(display, &event);
        switch (event.type) {
        case Expose:
        case MapNotify:
            many_pops_do_expose(event);
            expose_choice(choice1, choice2, string, c1, c2, wm,
                          event.xexpose.window);
            break;

        case ButtonPress:
            if (event.xbutton.window == c1) {
                value = (int32)key[0];
                not_done = 0;
            }
            if (event.xbutton.window == c2) {
                value = (int32)key[1];
                not_done = 0;
            }
            break;
        case KeyPress:
            value = ggets_get_key_press(&event);
            not_done = 0;
            break;
        case EnterNotify:
            if (event.xcrossing.window == c1 || event.xcrossing.window == c2)
                XSetWindowBorderWidth(display, event.xcrossing.window, 2);
            XFlush(display);
            break;
        case LeaveNotify:
            if (event.xcrossing.window == c1 || event.xcrossing.window == c2)
                XSetWindowBorderWidth(display, event.xcrossing.window, 1);
            XFlush(display);
            break;
        default:
            break;
        }
    }
    browse_wait_a_sec(2*ClickTime);
    XFlush(display);
    XSelectInput(display, c1, EV_MASK);
    XSelectInput(display, c2, EV_MASK);
    XFlush(display);
    XDestroySubwindows(display, base);
    XDestroyWindow(display, base);
    return value;
}

int32
pop_list_yes_no_box(void) {
    char ans;
    ans = (char)menudrive_two_choice("YES", "NO", "Are you sure?", "yn");
    if (ans == 'y')
        return 1;
    return 0;
}

/*  new pop_up_list   */

int32
pop_up_list(Window *root, char *title, char **list, char *key, int32 n,
            int32 max, int32 def, int32 x, int32 y, char **hints, Window hwin,
            char *httxt)

{
    POP_UP p;
    XEvent event;
    Window window;
    Cursor txt;
    int32 i, done = 0, value;
    int32 width = DCURX*(max + 5);
    int32 length = (DCURY + 6)*(n + 2);
    window = make_plain_window(*root, x, y, width, length, 2);
    txt = XCreateFontCursor(display, XC_hand2);
    XDefineCursor(display, window, txt);
    p.base = window;
    p.entries = list;
    p.title = title;
    p.n = n;
    p.hints = hints;
    p.max = max;
    p.key = key;
    p.hot = def;
    value = (int32)key[def];
    p.w = xmalloc((usize)n*sizeof(*(p.w)));
    p.tit = make_window(window, 0, 0, width, DCURY + 7, 0);
    for (i = 0; i < n; i++) {
        p.w[i] = make_window(window, DCURX, DCURY + 10 + i*(DCURY + 6),
                             DCURX*(max + 3), DCURY + 3, 0);
        XSelectInput(display, p.w[i], BUT_MASK);
    }

    while (!done) {
        XNextEvent(display, &event);
        switch (event.type) {
        case Expose:
        case MapNotify:
            many_pops_do_expose(event);
            draw_pop_up(p, event.xexpose.window);
            break;
        case KeyPress:
            value = ggets_get_key_press(&event);
            done = 1;
            break;
        case ButtonPress:
            for (i = 0; i < n; i++) {
                if (event.xbutton.window == p.w[i]) {
                    value = (int32)p.key[i];
                    done = 1;
                }
            }

            break;
        case EnterNotify:
            for (i = 0; i < p.n; i++)
                if (event.xcrossing.window == p.w[i]) {
                    XSetWindowBorderWidth(display, p.w[i], 1);
                    if (TipsFlag) {
                        strcpy(httxt, hints[i]);
                        XClearWindow(display, hwin);
                        XDrawString(display, hwin, gc, 5, CURY_OFF, hints[i],
                                    (int)strlen(hints[i]));
                    }
                }

            break;
        case LeaveNotify:
            for (i = 0; i < p.n; i++)
                if (event.xcrossing.window == p.w[i])
                    XSetWindowBorderWidth(display, p.w[i], 0);
            break;
        default:
            break;
        }
    }

    for (i = 0; i < n; i++)
        XSelectInput(display, p.w[i], EV_MASK);
    /*browse_wait_a_sec(ClickTime); Not here. Don't want to delay short cuts*/
    XDestroySubwindows(display, p.base);
    XDestroyWindow(display, p.base);
    XFlush(display);
    if (value == 13)
        value = (int32)key[def];
    return value;
}

void
draw_pop_up(POP_UP p, Window window) {
    int32 i;

    if (window == p.tit) {
        ggets_set_fore();
        ggets_bar(0, 0, DCURX*(p.max + 5), (DCURY + 7), window);
        ggets_set_back();
        ggets_f_text(DCURX*2, 4, p.title, window);
        ggets_set_fore();
        return;
    }
    for (i = 0; i < p.n; i++) {
        if (window == p.w[i]) {
            ggets_f_text(DCURX / 2, 3, p.entries[i], window);
            if (i == p.hot)
                ggets_f_text(DCURX*(p.max + 1), 4, "X", window);
            return;
        }
    }
    return;
}
