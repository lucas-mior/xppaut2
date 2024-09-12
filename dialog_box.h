
#ifndef _dialog_box_h
#define _dialog_box_h
#include "integers.h"

#include <X11/Xlib.h>
#include "struct.h"

int32 get_dialog(char *wname, char *name, char *value, char *ok, char *cancel,
               int32 max);
int32 dialog_event_loop(DIALOG *d, int32 max, int32 *pos, int32 *col);
void display_dialog(Window w, DIALOG d, int32 pos, int32 col);

#endif
