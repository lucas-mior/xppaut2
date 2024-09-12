#ifndef _kinescope_h_
#define _kinescope_h_
#include "integers.h"

void do_movie_com(int32 c);
void reset_film(void);
int32 film_clip(void);
int32 show_frame(int32 i, int32 h, int32 w);
void play_back(void);
void save_kine(void);
void make_anigif(void);
void save_movie(char *basename, int32 fmat);
void auto_play(void);
void too_small(void);

#endif
