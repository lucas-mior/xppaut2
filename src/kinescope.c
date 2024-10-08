#include "functions.h"
#include "integers.h"
#include <stdbool.h>

/* Kinescope for X windows */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <stdio.h>
#include <sys/time.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include "mykeydef.h"

#define MAXFILM 250

static struct Movie {
    int32 h;
    int32 w;
    Pixmap xi;
} movie[MAXFILM];

static void kinescope_too_small(void);
static void kinescope_auto_play(void);
static void kinescope_save_movie(char *basename, int32 fmat);
static void kinescope_make_anigif(void);
static void kinescope_save_kine(void);
static void kinescope_play_back(void);
static int32 kinescope_show_frame(int32 i, int32 h, int32 w);

void
kinescope_do_movie_com(int32 c) {
    switch (c) {
    case 0:
        if (kinescope_film_clip() == 0) {
            pop_list_respond_box("Okay", "Out of film!");
        }
        break;
    case 1:
        kinescope_reset_film();
        break;
    case 2:
        kinescope_play_back();
        break;
    case 3:
        kinescope_auto_play();
        break;
    case 4:
        kinescope_save_kine();
        break;
    case 5:
        kinescope_make_anigif();
        break;
    default:
        break;
    }
    return;
}

void
kinescope_reset_film(void) {
    if (mov_ind == 0) {
        return;
    }
    for (int32 i = 0; i < mov_ind; i++) {
        XFreePixmap(display, movie[i].xi);
    }
    mov_ind = 0;
    return;
}

int32
kinescope_film_clip(void) {
    int32 x;
    int32 y;
    uint32 h;
    uint32 w;
    uint32 bw;
    uint32 d;
    Window root;

    if (mov_ind >= MAXFILM) {
        return 0;
    }
    XGetGeometry(display, draw_win, &root, &x, &y, &w, &h, &bw, &d);
    movie[mov_ind].h = (int32)h;
    movie[mov_ind].w = (int32)w;
    movie[mov_ind].xi = XCreatePixmap(display, RootWindow(display, screen), w, h,
                                      (uint)DefaultDepth(display, screen));
    XCopyArea(display, draw_win, movie[mov_ind].xi, gc_graph, 0, 0, w, h, 0, 0);
    mov_ind++;
    return 1;
}

int32
kinescope_show_frame(int32 i, int32 h, int32 w) {
    if (h < movie[i].h || w < movie[i].w) {
        kinescope_too_small();
        return 1;
    }
    XCopyArea(display, movie[i].xi, draw_win, gc_graph, 0, 0, (uint)w, (uint)h, 0, 0);
    XFlush(display);

    return 0;
}

void
kinescope_play_back(void) {
    int32 x;
    int32 y;
    int32 h;
    int32 w;
    int32 bw;
    int32 d;

    Window root;
    XEvent event;
    int32 i = 0;
    XGetGeometry(display, draw_win, &root, &x, &y, (uint32 *)&w, (uint32 *)&h, (uint32 *)&bw,
                 (uint32 *)&d);
    if (mov_ind == 0) {
        return;
    }
    if (h < movie[i].h || w < movie[i].w) {
        kinescope_too_small();
        return;
    }

    XCopyArea(display, movie[i].xi, draw_win, gc_graph, 0, 0, (uint)w, (uint)h, 0, 0);
    XFlush(display);
    while (true) {
        XNextEvent(display, &event);
        switch (event.type) {
        case ButtonPress:
            i++;
            if (i >= mov_ind) {
                i = 0;
            }
            if (kinescope_show_frame(i, h, w)) {
                return;
            }
            break;
        case KeyPress:
            switch (ggets_get_key_press(&event)) {
            case KEY_ESC:
                return;
            case KEY_RIGHT:
                i++;
                if (i >= mov_ind) {
                    i = 0;
                }
                if (kinescope_show_frame(i, h, w)) {
                    return;
                }
                break;
            case KEY_LEFT:
                i--;
                if (i < 0) {
                    i = mov_ind - 1;
                }
                if (kinescope_show_frame(i, h, w)) {
                    return;
                }
                break;
            case KEY_HOME:
                i = 0;
                if (kinescope_show_frame(i, h, w)) {
                    return;
                }
                break;
            case KEY_END:
                i = mov_ind - 1;
                if (kinescope_show_frame(i, h, w)) {
                    return;
                }
                break;
            default:
                break;
            }
            break;
        default:
            break;
        }
    }
}

void
kinescope_save_kine(void) {
    char base[128];
    int32 fmat = 2;
    sprintf(base, "frame");
    ggets_new_string("Base file name", base);
    if (strlen(base) > 0) {
        kinescope_save_movie(base, fmat);
    }
    return;
}

void
kinescope_make_anigif(void) {
    int32 i = 0;
    int32 x;
    int32 y;
    FILE *fp;
    Window root;
    int32 h;
    int32 w;
    int32 bw;
    int32 d;
    XGetGeometry(display, draw_win, &root, &x, &y, (uint32 *)&w, (uint32 *)&h, (uint32 *)&bw,
                 (uint32 *)&d);
    if (mov_ind == 0) {
        return;
    }
    if (h < movie[i].h || w < movie[i].w) {
        kinescope_too_small();
        return;
    }
    h = movie[0].h;
    w = movie[0].w;
    for (i = 0; i < mov_ind; i++) {
        if ((movie[i].h != h) || (movie[i].w != w)) {
            ggets_err_msg("All clips must be same size");
            return;
        }
    }
    fp = fopen("anim.gif", "wb");
    scrngif_set_global_map(1);
    for (i = 0; i < mov_ind; i++) {
        XCopyArea(display, movie[i].xi, draw_win, gc_graph, 0, 0, (uint)w, (uint)h, 0, 0);
        XFlush(display);
        scrngif_add_ani_gif(movie[i].xi, fp, i);
    }

    scrngif_end_ani_gif(fp);
    fclose(fp);
    scrngif_set_global_map(0);
    return;
}

void
kinescope_save_movie(char *basename, int32 fmat) {
    char file[XPP_MAX_NAME];
    int32 i = 0;
    int32 x;
    int32 y;
    FILE *fp;
    Window root;
    Pixmap xi;
    int32 h;
    int32 w;
    int32 bw;
    int32 d;
    XGetGeometry(display, draw_win, &root, &x, &y, (uint32 *)&w, (uint32 *)&h, (uint32 *)&bw,
                 (uint32 *)&d);
    if (mov_ind == 0) {
        return;
    }
    if (h < movie[i].h || w < movie[i].w) {
        kinescope_too_small();
        return;
    }

    for (i = 0; i < mov_ind; i++) {
        if (fmat == 1) {
            sprintf(file, "%s_%d.ppm", basename, i);
        } else {
            sprintf(file, "%s_%d.gif", basename, i);
        }
        XCopyArea(display, movie[i].xi, draw_win, gc_graph, 0, 0, (uint)w, (uint)h, 0, 0);
        XFlush(display);
        if (fmat == 1) {
            ani_write_frame(file, draw_win, w, h);
        }
#ifndef NOGIF
        else {
            XGetGeometry(display, draw_win, &root, &x, &y, (uint32 *)&w, (uint32 *)&h,
                         (uint32 *)&bw, (uint32 *)&d);
            xi = XCreatePixmap(display, RootWindow(display, screen), (uint)w, (uint)h,
                               (uint)DefaultDepth(display, screen));
            XCopyArea(display, draw_win, xi, gc_graph, 0, 0, (uint)w, (uint)h, 0, 0);

            fp = fopen(file, "wb");
            scrngif_screen_to_gif(xi, fp);
            fclose(fp);
        }
#endif
    }
    return;
}

void
kinescope_auto_play(void) {
    static int32 ks_ncycle = 1;
    static int32 ks_speed = 50;
    int32 x;
    int32 y;
    int32 h;
    int32 w;
    int32 bw;
    int32 d;
    int32 key;
    Window root;

    int32 dt = 20;
    int32 smax = 500;
    XEvent event;
    int32 i = 0;
    int32 cycle = 0;

    ggets_new_int("Number of cycles", &ks_ncycle);
    ggets_new_int("Msec between frames", &ks_speed);
    if (ks_speed < 0) {
        ks_speed = 0;
    }
    if (ks_ncycle <= 0) {
        return;
    }
    XGetGeometry(display, draw_win, &root, &x, &y, (uint32 *)&w, (uint32 *)&h, (uint32 *)&bw,
                 (uint32 *)&d);
    if (mov_ind == 0) {
        return;
    }
    if (h < movie[i].h || w < movie[i].w) {
        kinescope_too_small();
        return;
    }

    XCopyArea(display, movie[i].xi, draw_win, gc_graph, 0, 0, (uint)w, (uint)h, 0, 0);
    XFlush(display);

    while (true) {
        // check for events
        if (XPending(display) > 0) {
            XNextEvent(display, &event);
            switch (event.type) {
            case ButtonPress:
                return;
            case KeyPress:
                key = ggets_get_key_press(&event);
                if (key == PAUSE_NUMBER) {
                    return;
                }
                if (key == ',') {
                    ks_speed -= dt;
                    if (ks_speed < dt) {
                        ks_speed = dt;
                    }
                }
                if (key == '.') {
                    ks_speed += dt;
                    if (ks_speed > smax) {
                        ks_speed = smax;
                    }
                }

                break;
            default:
                break;
            }
        }  // done checking  now increment pix

        browser_wait_a_sec(ks_speed);
        i++;
        if (i >= mov_ind) {
            cycle++;
            i = 0;
        }
        if (h < movie[i].h || w < movie[i].w) {
            kinescope_too_small();
            return;
        }
        XCopyArea(display, movie[i].xi, draw_win, gc_graph, 0, 0, (uint)w, (uint)h, 0, 0);
        XFlush(display);
        if (cycle >= ks_ncycle) {
            return;
        }
    }
}

void
kinescope_too_small(void) {
    pop_list_respond_box("Okay", "Window too small for film!");
    return;
}
