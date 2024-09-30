#include "functions.h"
#include "parserslow.h"
#include "integers.h"
#include <stdbool.h>

#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include <X11/cursorfont.h>
#ifndef WCTYPE
#include <ctype.h>
#else
#include <wctype.h>
#endif

#define NEQMAXFOREDIT 20
#define MAX_N_EBOX MAX_ODE
#define MAX_LEN_EBOX 86

typedef struct EditBox {
    Window base;
    Window ok;
    Window cancel;
    Window reset;
    Window win[MAX_N_EBOX];
    char name[MAX_N_EBOX][MAX_LEN_EBOX];
    char value[MAX_N_EBOX][MAX_LEN_EBOX];
    char rval[MAX_N_EBOX][MAX_LEN_EBOX];
    int32 n;
    int32 hot;
} EditBox;

static void edit_rhs_reset_box(EditBox *sb, int32 *pos, int32 *col);
static void edit_rhs_expose_box(EditBox *sb, Window window, int32 pos);
static void edit_rhs_reset_hot(int32 inew, EditBox *sb);
static void edit_rhs_new_editable(EditBox *sb, int32 inew, int32 *pos,
                                  int32 *col, int32 *done, Window *w);
static int32 edit_rhs_box_event_loop(EditBox *sb, int32 *pos, int32 *col);
static void edit_rhs_make_box_windows(EditBox *sb, char *title);

static int32 edit_rhs_box(int32 n, char *title, char **names, char **values);

void
edit_rhs_reset_box(EditBox *sb, int32 *pos, int32 *col) {
    int32 n = sb->n;
    int32 l;
    Window window;
    for (int32 i = 0; i < n; i++) {
        strcpy(sb->value[i], sb->rval[i]);
        window = sb->win[i];
        l = (int32)strlen(sb->name[i]);
        XClearWindow(display, window);
        XDrawString(display, window, gc, 0, cury_off, sb->name[i], l);
        XDrawString(display, window, gc, l*dcur_x, cury_off, sb->value[i],
                    (int32)strlen(sb->value[i]));
    }
    XFlush(display);
    sb->hot = 0;
    *pos = (int32)strlen(sb->value[0]);
    *col = (*pos + (int32)strlen(sb->name[0]))*dcur_x;
    ggets_put_cursor_at(sb->win[0], dcur_x*(int32)strlen(sb->name[0]), *pos);
    return;
}

int32
edit_rhs_box(int32 n, char *title, char **names, char **values) {
    EditBox sb;
    int32 status;
    int32 colm;
    int32 pos;

    for (int32 i = 0; i < n; i++) {
        sprintf(sb.name[i], "%s=", names[i]);
        strcpy(sb.value[i], values[i]);
        strcpy(sb.rval[i], values[i]);
    }
    sb.n = n;
    sb.hot = 0;
    edit_rhs_make_box_windows(&sb, title);
    XSelectInput(display, sb.cancel, MASK_BUTTON);
    XSelectInput(display, sb.ok, MASK_BUTTON);
    XSelectInput(display, sb.reset, MASK_BUTTON);
    pos = (int32)strlen(sb.value[0]);
    colm = (pos + (int32)strlen(sb.name[0]))*dcur_x;

    while (true) {
        status = edit_rhs_box_event_loop(&sb, &pos, &colm);
        if (status != -1) {
            break;
        }
    }
    XSelectInput(display, sb.cancel, MASK_EVENT);
    XSelectInput(display, sb.ok, MASK_EVENT);
    XSelectInput(display, sb.reset, MASK_EVENT);

    browser_wait_a_sec(CLICK_TIME);
    XDestroySubwindows(display, sb.base);
    XDestroyWindow(display, sb.base);

    if (status == ALL_FORGET) {
        return status;
    }
    for (int32 i = 0; i < n; i++) {
        strcpy(values[i], sb.value[i]);
    }
    return status;
}

void
edit_rhs_expose_box(EditBox *sb, Window window, int32 pos) {
    int32 flag;

    if (window == sb->ok) {
        XDrawString(display, window, gc, 0, cury_off, "Ok", 2);
        return;
    }
    if (window == sb->cancel) {
        XDrawString(display, window, gc, 0, cury_off, "Cancel", 6);
        return;
    }
    if (window == sb->reset) {
        XDrawString(display, window, gc, 0, cury_off, "Reset", 5);
        return;
    }
    for (int32 i = 0; i < sb->n; i++) {
        if (window != sb->win[i]) {
            continue;
        }
        flag = 0;
        if (i == sb->hot) {
            flag = 1;
        }
        pop_list_do_hilite_text(sb->name[i], sb->value[i], flag, window, pos);
    }
    return;
}

void
edit_rhs_reset_hot(int32 inew, EditBox *sb) {
    int32 i = sb->hot;
    sb->hot = inew;
    XClearWindow(display, sb->win[inew]);
    pop_list_do_hilite_text(sb->name[inew], sb->value[inew], 1, sb->win[inew],
                            (int32)strlen(sb->value[inew]));
    XClearWindow(display, sb->win[i]);
    pop_list_do_hilite_text(sb->name[i], sb->value[i], 0, sb->win[i],
                            (int32)strlen(sb->value[i]));
    return;
}

void
edit_rhs_new_editable(EditBox *sb, int32 inew, int32 *pos, int32 *col,
                      int32 *done, Window *w) {
    edit_rhs_reset_hot(inew, sb);
    *pos = (int32)strlen(sb->value[inew]);
    *col = (*pos + (int32)strlen(sb->name[inew]))*dcur_x;
    *done = 0;
    *w = sb->win[inew];
    return;
}

int32
edit_rhs_box_event_loop(EditBox *sb, int32 *pos, int32 *col) {
    XEvent event;
    int32 status = -1;
    int32 inew;
    int32 nn = sb->n;
    int32 done = 0;
    char ch;
    int32 ihot = sb->hot;
    Window wt;
    Window window = sb->win[ihot];  // active window
    char *s;
    s = sb->value[ihot];

    XNextEvent(display, &event);
    switch (event.type) {
    case ConfigureNotify:
    case Expose:
    case MapNotify:
        many_pops_do_expose(event);  //  menus and graphs etc
        edit_rhs_expose_box(sb, event.xany.window, *pos);
        break;

    case ButtonPress:
        if (event.xbutton.window == sb->ok) {
            status = ALL_DONE;
            break;
        }
        if (event.xbutton.window == sb->cancel) {
            status = ALL_FORGET;
            break;
        }
        if (event.xbutton.window == sb->reset) {
            edit_rhs_reset_box(sb, pos, col);
            break;
        }
        for (int32 i = 0; i < nn; i++) {
            if (event.xbutton.window == sb->win[i]) {
                XSetInputFocus(display, sb->win[i], RevertToParent,
                               CurrentTime);
                if (i != sb->hot) {
                    edit_rhs_new_editable(sb, i, pos, col, &done, &window);
                }
                break;
            }
        }
        break;

    case EnterNotify:
        wt = event.xcrossing.window;
        if (wt == sb->ok || wt == sb->cancel || wt == sb->reset) {
            XSetWindowBorderWidth(display, wt, 2);
        }
        break;

    case LeaveNotify:
        wt = event.xcrossing.window;
        if (wt == sb->ok || wt == sb->cancel || wt == sb->reset) {
            XSetWindowBorderWidth(display, wt, 1);
        }
        break;

    case KeyPress:
        ch = (char)ggets_get_key_press(&event);
        ggets_edit_window(window, pos, s, col, &done, ch);
        if (done != 0) {
            if (done == ALL_DONE) {
                status = ALL_DONE;
                break;
            }
            inew = (sb->hot + 1) % nn;
            edit_rhs_new_editable(sb, inew, pos, col, &done, &window);
        }
        break;
    default:
        break;
    }
    return status;
}

void
edit_rhs_make_box_windows(EditBox *sb, char *title) {
    int32 width;
    int32 height;
    int32 xpos;
    int32 ypos;
    int32 n = sb->n;
    int32 xstart;
    int32 ystart;

    XTextProperty winname;
    XSizeHints size_hints;
    Window base;
    width = (MAX_LEN_EBOX + 4)*dcur_x;
    height = (n + 4)*(dcur_y + 16);
    base = pop_list_make_plain_window(DefaultRootWindow(display), 0, 0, width,
                                      height, 4);
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
    XSetWMProperties(display, base, &winname, NULL, NULL, 0, &size_hints, NULL,
                     NULL);
    sb->base = base;

    ystart = dcur_y;
    xstart = dcur_x;
    for (int32 i = 0; i < n; i++) {
        xpos = xstart;
        ypos = ystart + i*(dcur_y + 10);
        sb->win[i] = pop_list_make_window(base, xpos, ypos,
                                          MAX_LEN_EBOX*dcur_x, dcur_y, 1);
    }

    ypos = height - 2*dcur_y;
    xpos = (width - 19*dcur_x) / 2;
    (sb->ok) = pop_list_make_window(base, xpos, ypos, 2*dcur_x, dcur_y, 1);
    (sb->cancel) = pop_list_make_window(base, xpos + 4*dcur_x, ypos,
                                        6*dcur_x, dcur_y, 1);
    (sb->reset) = pop_list_make_window(base, xpos + 12*dcur_x, ypos,
                                       5*dcur_x, dcur_y, 1);
    XRaiseWindow(display, base);
    return;
}

void
edit_rhs_menu(void) {
    Window temp = main_win;
    static char *n[] = {"RHS's", "Functions", "Save as", "Load DLL"};
    static char key[] = "rfsl";
    char ch;
    int32 edtype = 0;
    ch = (char)pop_list_popup_list_new(&temp, "Edit Stuff", n, key, 4, 11,
                                       edtype, 10, 13*dcur_y + 8, edrh_hint,
                                       info_pop, info_message);
    edtype = -1;
    for (int32 i = 0; i < 4; i++) {
        if (ch == key[i]) {
            edtype = i;
        }
    }
    switch (edtype) {
    case 0:
        edit_rhs();
        break;
    case 1:
        edit_rhs_functions();
        break;
    case 2:
        edit_rhs_save_as();
        break;
    case 3:
        extra_load_new_dll();
        break;
    default:
        fprintf(stderr, "Unexpected switch case in %s.\n", __func__);
        exit(EXIT_FAILURE);
    }
    return;
}

void
edit_rhs(void) {
    char **names, **values;
    int32 **command;
    int32 status;
    int32 err;
    int32 len;
    int32 i0;
    int32 n = NEQ;
    char fstr[20];
    char msg[200];
    if (NEQ > NEQMAXFOREDIT) {
        return;
    }
    names = xmalloc((usize)n*sizeof(char *));
    values = xmalloc((usize)n*sizeof(char *));
    command = xmalloc((usize)n*sizeof(int32 *));
    for (int32 i = 0; i < n; i++) {
        values[i] = xmalloc(MAX_LEN_EBOX*sizeof(*(values[i])));
        names[i] = xmalloc(MAX_LEN_EBOX*sizeof(*(names[i])));
        command[i] = xmalloc(200*sizeof(*(command[i])));
        if (i < NODE && METHOD > 0) {
            strcpy(fstr, "d%s/dT");
        }
        if (i < NODE && METHOD == 0) {
            strcpy(fstr, "%s(n+1)");
        }
        if (i < NODE && eq_type[i] == 1) {
            strcpy(fstr, "%s(T)");
        }
        if (i >= NODE) {
            strcpy(fstr, "%s");
        }
        sprintf(names[i], fstr, uvar_names[i]);
        strcpy(values[i], ode_names[i]);
    }
    status = edit_rhs_box(n, "Right Hand Sides", names, values);
    if (status != 0) {
        for (int32 i = 0; i < n; i++) {
            if (i < NODE || (i >= (NODE + NMarkov))) {
                err = parserslow_add_expr(values[i], command[i], &len);
                if (err == 1) {
                    sprintf(msg, "Bad rhs:%s=%s", names[i], values[i]);
                    ggets_err_msg(msg);
                } else {
                    free(ode_names[i]);
                    ode_names[i] = xmalloc(strlen(values[i]) + 5);
                    strcpy(ode_names[i], values[i]);
                    i0 = i;
                    if (i >= NODE) {
                        i0 = i0 + fix_var - NMarkov;
                    }

                    for (int32 j = 0; j < len; j++) {
                        my_ode[i0][j] = command[i][j];
                    }
                }
            }
        }
    }

    for (int32 i = 0; i < n; i++) {
        free(values[i]);
        free(names[i]);
        free(command[i]);
    }
    free(values);
    free(names);
    free(command);
    return;
}

void
edit_rhs_user_fun_info(FILE *fp) {
    char fundef[256];
    for (int32 j = 0; j < NFUN; j++) {
        sprintf(fundef, "%s(", ufun_names[j]);
        for (int32 i = 0; i < narg_fun[j]; i++) {
            strcat(fundef, ufun_arg[j].args[i]);
            if (i < narg_fun[j] - 1) {
                strcat(fundef, ",");
            }
        }
        strcat(fundef, ") = ");
        strcat(fundef, ufun_def[j]);
        fprintf(fp, "%s\n", fundef);
    }
    return;
}

void
edit_rhs_functions(void) {
    char **names, **values;
    int32 **command;
    int32 status;
    int32 err;
    int32 len;
    int32 n = NFUN;
    char msg[200];
    if (n == 0 || n > NEQMAXFOREDIT) {
        return;
    }
    names = xmalloc((usize)n*sizeof(char *));
    values = xmalloc((usize)n*sizeof(char *));
    command = xmalloc((usize)n*sizeof(int32 *));
    for (int32 i = 0; i < n; i++) {
        values[i] = xmalloc(MAX_LEN_EBOX*sizeof(*(values[i])));
        names[i] = xmalloc(MAX_LEN_EBOX*sizeof(*(names[i])));
        command[i] = xmalloc(200*sizeof(*(command[i])));
        sprintf(values[i], "%s", ufun_def[i]);

        if (narg_fun[i] == 0) {
            sprintf(names[i], "%s()", ufun_names[i]);
        }
        if (narg_fun[i] == 1) {
            sprintf(names[i], "%s(%s)", ufun_names[i], ufun_arg[i].args[0]);
        }
        if (narg_fun[i] > 1) {
            sprintf(names[i], "%s(%s,...,%s)", ufun_names[i],
                    ufun_arg[i].args[0], ufun_arg[i].args[narg_fun[i] - 1]);
        }
    }

    status = edit_rhs_box(n, "Functions", names, values);
    if (status != 0) {
        for (int32 i = 0; i < n; i++) {
            set_new_arg_names(narg_fun[i], ufun_arg[i].args);
            err = parserslow_add_expr(values[i], command[i], &len);
            set_old_arg_names(narg_fun[i]);
            if (err == 1) {
                sprintf(msg, "Bad func.:%s=%s", names[i], values[i]);
                ggets_err_msg(msg);
            } else {
                strcpy(ufun_def[i], values[i]);
                for (int32 j = 0; j <= len; j++) {
                    ufun[i][j] = command[i][j];
                }
                fixup_endfun(ufun[i], len, narg_fun[i]);
            }
        }
    }

    for (int32 i = 0; i < n; i++) {
        free(values[i]);
        free(names[i]);
        free(command[i]);
    }
    free(values);
    free(names);
    free(command);
    return;
}

int32
edit_rhs_save_as(void) {
    int32 ok;
    FILE *fp;
    double z;
    char filename[sizeof(this_file)];
    strncpy(filename, this_file, sizeof(filename));
    ggets_ping();
    if (!init_conds_file_selector("Save As", filename, "*.ode")) {
        return -1;
    }
    browser_open_write_file(&fp, filename, &ok);
    if (!ok) {
        return -1;
    }
    fp = fopen(filename, "w");
    if (fp == NULL) {
        return 0;
    }
    fprintf(fp, "%d", NEQ);
    for (int32 i = 0; i < NODE; i++) {
        if (i % 5 == 0) {
            fprintf(fp, "\nvariable ");
        }
        fprintf(fp, " %s=%.16g ", uvar_names[i], last_ic[i]);
    }
    fprintf(fp, "\n");
    for (int32 i = NODE; i < NEQ; i++) {
        if ((i - NODE) % 5 == 0) {
            fprintf(fp, "\naux ");
        }
        fprintf(fp, " %s ", uvar_names[i]);
    }
    fprintf(fp, "\n");
    for (int32 i = 0; i < NUPAR; i++) {
        if (i % 5 == 0) {
            fprintf(fp, "\nparam  ");
        }
        get_val(upar_names[i], &z);
        fprintf(fp, " %s=%.16g   ", upar_names[i], z);
    }
    fprintf(fp, "\n");
    for (int32 i = 0; i < NFUN; i++) {
        fprintf(fp, "user %s %d %s\n", ufun_names[i], narg_fun[i], ufun_def[i]);
    }
    for (int32 i = 0; i < NODE; i++) {
        if (eq_type[i] == 1) {
            fprintf(fp, "i ");
        } else {
            fprintf(fp, "o ");
        }
        fprintf(fp, "%s\n", ode_names[i]);
    }
    for (int32 i = NODE; i < NEQ; i++) {
        fprintf(fp, "o %s\n", ode_names[i]);
    }
    for (int32 i = 0; i < NODE; i++) {
        fprintf(fp, "b %s \n", my_bc[i].string);
    }
    fprintf(fp, "done\n");
    fclose(fp);

    return 1;
}
