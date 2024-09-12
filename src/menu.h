#ifndef _xppmenu_h_
#define _xppmenu_h_
#include "integers.h"

#include <X11/Xlib.h>

void flash(int32 num);
void add_menu(Window base, int32 j, int32 n, char **names, char *key,
              char **hint);
void create_the_menus(Window base);
void show_menu(int32 j);
void unshow_menu(int32 j);
void help(void);
void help_num(void);
void help_file(void);
void menu_crossing(Window win, int32 yn);
void menu_expose(Window win);
void menu_button(Window win);
void draw_help(void);

#endif
