#ifndef _pop_list_h
#define _pop_list_h
#include "integers.h"

#include "phsplan.h"
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <stdio.h>
#include "xpplim.h"
#include "math.h"

#define MAX_N_SBOX 22
#include "max_len_sbox.h"

#define FORGET_ALL 0
#define DONE_ALL 2
#define FORGET_THIS 3
#define DONE_THIS 1

#define EV_MASK                                                                \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask)

#define BUT_MASK                                                               \
    (ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |     \
     EnterWindowMask | LeaveWindowMask)

extern Display *display;
extern int32 DisplayWidth, DisplayHeight;
extern int32 screen;
extern Atom deleteWindowAtom;
extern Window main_win, info_pop, draw_win;
extern int32 DCURY, DCURX, CURY_OFF, DCURXs, DCURYs, CURY_OFFs, xor_flag;
extern GC gc, small_gc;
extern uint32 MyBackColor, MyForeColor;

extern int32 TipsFlag;
char *get_next(), *get_first();
extern char UserBlack[8];
extern char UserWhite[8];
extern int32 UserGradients;

Window make_window();
Window make_plain_window();

/*  This is a string box widget which handles a list of
        editable strings
 */

typedef struct {
    Window base, ok, cancel;
    Window win[MAX_N_SBOX];
    char name[MAX_N_SBOX][MAX_LEN_SBOX], value[MAX_N_SBOX][MAX_LEN_SBOX];
    int32 n, hot;
    int32 hgt, wid;
    int32 hh[MAX_N_SBOX];
} STRING_BOX;

typedef struct {
    char **list;
    int32 n;
} SCRBOX_LIST;

extern int32 NUPAR, NEQ, NODE, NMarkov;
extern char upar_names[MAXPAR][14], uvar_names[MAXODE][12];
extern char *color_names[];
extern SCRBOX_LIST scrbox_list[10];

/*  This is a new improved pop_up widget */
typedef struct {
    Window base, tit;
    Window *w;
    char *title;
    char **entries;
    char **hints;
    int32 n, max;
    char *key;
    int32 hot;
} POP_UP;

typedef struct {
    Window base, slide, close, text;
    int32 i0;
    int32 exist, len, nlines;
    char **list;
} TEXTWIN;

typedef struct {
    Window base, slide;
    Window *w;
    int32 nw, nent, i0;
    int32 len, exist;
    char **list;
} SCROLLBOX;

extern TEXTWIN mytext;
#define SB_PLOTTABLE 0
#define SB_VARIABLE 1
#define SB_PARAMETER 2
#define SB_PARVAR 3
#define SB_COLOR 4
#define SB_MARKER 5
#define SB_METHOD 6

void set_window_title(Window win, char *string);
void make_scrbox_lists(void);
int32 get_x_coord_win(Window win);
void destroy_scroll_box(SCROLLBOX *sb);
void create_scroll_box(Window root, int32 x0, int32 y0, int32 nent, int32 nw,
                       char **list, SCROLLBOX *sb);
void expose_scroll_box(Window w, SCROLLBOX sb);
void redraw_scroll_box(SCROLLBOX sb);
void crossing_scroll_box(Window w, int32 c, SCROLLBOX sb);
int32 scroll_box_motion(XEvent ev, SCROLLBOX *sb);
int32 select_scroll_item(Window w, SCROLLBOX sb);
void scroll_popup(STRING_BOX *sb, SCROLLBOX *scrb);
int32 do_string_box(int32 n, int32 row, int32 col, char *title, char **names,
                    char values[][MAX_LEN_SBOX], int32 maxchar);
void expose_sbox(STRING_BOX sb, Window w, int32 pos, int32 col);
void do_hilite_text(char *name, char *value, int32 flag, Window w, int32 pos,
                    int32 col);
void reset_hot(int32 inew, STRING_BOX *sb);
void new_editable(STRING_BOX *sb, int32 inew, int32 *pos, int32 *col,
                  int32 *done, Window *w);
void set_sbox_item(STRING_BOX *sb, int32 item);
int32 s_box_event_loop(STRING_BOX *sb, int32 *pos, int32 *col, SCROLLBOX *scrb);
void make_sbox_windows(STRING_BOX *sb, int32 row, int32 col, char *title,
                       int32 maxchar);
Window make_fancy_window(Window root, int32 x, int32 y, int32 width,
                         int32 height, int32 bw, int32 fc, int32 bc);
Window make_unmapped_window(Window root, int32 x, int32 y, int32 width,
                            int32 height, int32 bw);
Window make_plain_unmapped_window(Window root, int32 x, int32 y, int32 width,
                                  int32 height, int32 bw);
Window make_window(Window root, int32 x, int32 y, int32 width, int32 height,
                   int32 bw);
Window make_plain_window(Window root, int32 x, int32 y, int32 width,
                         int32 height, int32 bw);
void expose_resp_box(char *button, char *message, Window wb, Window wm,
                     Window w);
void respond_box(char *button, char *message);
void message_box(Window *w, int32 x, int32 y, char *message);
void expose_choice(char *choice1, char *choice2, char *msg, Window c1,
                   Window c2, Window wm, Window w);
int32 two_choice(char *choice1, char *choice2, char *string, char *key, int32 x,
                 int32 y, Window w, char *title);
int32 yes_no_box(void);
int32 pop_up_list(Window *root, char *title, char **list, char *key, int32 n,
                  int32 max, int32 def, int32 x, int32 y, char **hints,
                  Window hwin, char *httxt);
void draw_pop_up(POP_UP p, Window w);
Window make_unmapped_icon_window(Window root, int32 x, int32 y, int32 width,
                                 int32 height, int32 bw, int32 icx, int32 icy,
                                 unsigned char *icdata);
Window make_icon_window(Window root, int32 x, int32 y, int32 width,
                        int32 height, int32 bw, int32 icx, int32 icy,
                        unsigned char *icdata);

#endif
