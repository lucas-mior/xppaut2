#include "functions.h"
#include "integers.h"
#include <stdbool.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define MAXPSLINE 100

#define PS_XOFF 50
#define PS_YOFF 50
#define PS_XMAX 7200
#define PS_YMAX 5040

#define PS_VTIC (PS_YMAX / 80)
#define PS_HTIC (PS_YMAX / 80)

#define PS_SC (10)  // scale is 1pt = 10 units
#define PS_VCHAR (ps_fontsize*PS_SC)

#define LEFT 0
#define RIGHT 2
#define CENTER 1
#define POINT_TYPES 8
int32 last_pt_line;
int32 no_break_line = 0;
int32 ps_fontsize = 14;
double ps_lw = 5;
char ps_font[200] = "Times-Roman";
FILE *psfile;
/*Default is now with color*/
int32 plt_fmt_flag;
int32 ps_color_flag = 1;
int32 ps_lines;
int32 last_psx;
int32 last_psy;

/* this header stuff was stolen from GNUPLOT I have added  filled circles
 * and open circles for bifurcation diagrams I also use Times Roman
 * since Courier is an ugly font!!  */

static char *PS_header[] = {
    "/vpt2 vpt 2 mul def\n", "/hpt2 hpt 2 mul def\n",
    "/Romfnt {/Times-Roman findfont exch scalefont setfont} def ",
    "/Symfnt {/Symbol findfont exch scalefont setfont} def ",
    // flush left show
    "/Lshow { currentpoint stroke moveto\n", "  0 vshift rmoveto show } def\n",
    // flush right show
    "/Rshow { currentpoint stroke moveto\n",
    "  dup stringwidth pop neg vshift rmoveto show } def\n",
    // centred show
    "/Cshow { currentpoint stroke moveto\n",
    "  dup stringwidth pop -2 div vshift rmoveto show } def\n",
    // Dash or Color Line
    "/DL { Color {setrgbcolor [] 0 setdash pop}\n", " {pop pop pop 0 setdash} ifelse } def\n",
    // Border Lines
    "/BL { stroke xpplinewidth 2 mul setlinewidth } def\n",
    // Axes Lines
    "/AL { stroke xpplinewidth 2 div setlinewidth } def\n",
    // Plot Lines
    "/PL { stroke xpplinewidth setlinewidth } def\n",
    // Line Types
    "/LTb { BL [] 0 0 0 DL } def\n",                             // border
    "/LTa { AL [1 dl 2 dl] 0 setdash 0 0 0 setrgbcolor } def\n", /* axes
                                                                  */
    "/LT0 { PL [] 0 0 0 DL } def\n", "/LT1 { PL [4 dl 2 dl] 1 0 0 DL } def\n",
    "/LT2 { PL [2 dl 3 dl] .95 .4 0 DL } def\n", "/LT3 { PL [1 dl 1.5 dl] 1 .65 0 DL } def\n",
    "/LT4 { PL [5 dl 2 dl 1 dl 2 dl] 1 .8 0 DL } def\n",
    "/LT5 { PL [4 dl 3 dl 1 dl 3 dl] .85 .85 0 DL } def\n",
    "/LT6 { PL [2 dl 2 dl 2 dl 4 dl]  .6 .8 .2 DL } def\n",
    "/LT7 { PL [2 dl 2 dl 2 dl 2 dl 2 dl 4 dl] 0 .9 0 DL } def\n",
    "/LT8 { stroke 16. setlinewidth [] 0 .85 .85 DL } def\n", /* really
                                                                 fat line
                                                               */
    "/LT9 { stroke 16. setlinewidth [4 dl 2 dl] 0 0 1 DL } def\n",
    "/LTc { stroke 16. setlinewidth [2 dl 3 dl] .62 .125 .93 DL } def\n", "/M {moveto} def\n",
    "/L {lineto} def\n", "/R {rlineto} def\n",
    "/P { stroke [] 0 setdash\n",  // Point
    "  currentlinewidth 2 div sub moveto\n", "  0 currentlinewidth rlineto  stroke } def\n",
    "/D { stroke [] 0 setdash  2 copy  vpt add moveto\n",  // Diamond
    "  hpt neg vpt neg rlineto  hpt vpt neg rlineto\n",
    "  hpt vpt rlineto  hpt neg vpt rlineto  closepath  stroke\n", "  P  } def\n",
    "/A { stroke [] 0 setdash  vpt sub moveto  0 vpt2 rlineto\n", /* Plus
                                                                     (Add)
                                                                   */
    "  currentpoint stroke moveto\n", "  hpt neg vpt neg rmoveto  hpt2 0 rlineto stroke\n",
    "  } def\n",
    "/B { stroke [] 0 setdash  2 copy  exch hpt sub exch vpt add "
    "moveto\n",  // Box
    "  0 vpt2 neg rlineto  hpt2 0 rlineto  0 vpt2 rlineto\n",
    "  hpt2 neg 0 rlineto  closepath  stroke\n", "  P  } def\n",
    "/C { stroke [] 0 setdash  exch hpt sub exch vpt add moveto\n", /* Cross
                                                                     */
    "  hpt2 vpt2 neg rlineto  currentpoint  stroke  moveto\n",
    "  hpt2 neg 0 rmoveto  hpt2 vpt2 rlineto stroke  } def\n",
    "/T { stroke [] 0 setdash  2 copy  vpt 1.12 mul add moveto\n", /* Triangle
                                                                    */
    "  hpt neg vpt -1.62 mul rlineto\n", "  hpt 2 mul 0 rlineto\n",
    "  hpt neg vpt 1.62 mul rlineto  closepath  stroke\n", "  P  } def\n",
    "/S { 2 copy A C} def\n",                                    // Star
    "/K { stroke [] 0 setdash vpt 0 360 arc stroke} def ",       // circle2
    "/F { stroke [] 0 setdash vpt 0 360 arc fill stroke } def ", /* Filled
                                                                    circle
                                                                  */
    NULL};

static void ps_write(char *str);
static void ps_check_lines(void);
static void ps_rel(int32 x, int32 y);

int32
ps_init(char *filename, int32 color) {
    if ((psfile = fopen(filename, "w")) == NULL) {
        ggets_err_msg("Cannot open file ");
        return 0;
    }
    graphics_init_ps();
    plt_fmt_flag = 1;
    ps_lines = 0;
    last_psx = -10000;
    last_psy = -10000;
    fprintf(psfile, "%%!PS-Adobe-2.0\n");
    fprintf(psfile, "%%Creator: xppaut\n");
    fprintf(psfile, "%%%%BoundingBox: %d %d %d %d\n", PS_XOFF, PS_YOFF,
            (int32)(PS_YMAX / PS_SC + .5 + PS_YOFF + 0.1*PS_VCHAR),
            (int32)(PS_XMAX / PS_SC + .5 + PS_XOFF + 0.1*PS_VCHAR));
    fprintf(psfile, "/xppdict 40 dict def\nxppdict begin\n");
    if (color == 0) {
        fprintf(psfile, "/Color false def \n");
        ps_color_flag = 0;
    } else {
        fprintf(psfile, "/Color true def \n");
        fprintf(psfile, "/RGB {setrgbcolor currentpoint stroke moveto} def\n");
        fprintf(psfile, "/RGb {setrgbcolor } def\n");
        ps_color_flag = 1;
    }
    fprintf(psfile, "/xpplinewidth %.3f def\n", ps_lw);
    fprintf(psfile, "/vshift %d def\n", (int32)(PS_VCHAR) / (-3));
    fprintf(psfile, "/dl {%d mul} def\n", PS_SC);  // dash length
    fprintf(psfile, "/hpt %.1f def\n", PS_HTIC / 2.0);
    fprintf(psfile, "/vpt %.1f def\n", PS_VTIC / 2.0);
    for (int32 i = 0; PS_header[i] != NULL; i++) {
        fprintf(psfile, "%s", PS_header[i]);
    }
    fprintf(psfile, "end\n");
    fprintf(psfile, "%%%%EndProlog\n");
    fprintf(psfile, "xppdict begin\n");
    fprintf(psfile, "gsave\n");
    fprintf(psfile, "%d %d translate\n", PS_XOFF, PS_YOFF);
    fprintf(psfile, "%.3f %.3f scale\n", 1. / PS_SC, 1. / PS_SC);
    if (!ps_port) {
        fprintf(psfile, "90 rotate\n0 %d translate\n", -PS_YMAX);
    }
    fprintf(psfile, "/%s findfont %d ", ps_font, ps_fontsize*PS_SC);
    fprintf(psfile, "scalefont setfont\n");
    fprintf(psfile, "newpath\n");
    return 1;
}

void
ps_stroke(void) {
    fprintf(psfile, "stroke\n");
    return;
}

void
ps_do_color(int32 color) {
    double r;
    double g;
    double b;
    // this doesn work very well
    if (plt_fmt_flag == 0) {
        return;
    }
    if (ps_color_flag == 0) {
        return;
    }
    color_get_ps(color, &r, &g, &b);
    fprintf(psfile, "%f %f %f RGb\n", r, g, b);
    return;
}

void
ps_end(void) {
    ps_write("stroke");
    ps_write("grestore");
    ps_write("end");
    ps_write("showpage");
    lunch_ps_write_pars(psfile);
    fclose(psfile);
    plt_fmt_flag = 0;
    if (Xup) {
        graphics_init_x11();
    }
    return;
}

void
ps_frect(int32 x, int32 y, int32 w, int32 h) {
    fprintf(psfile, " newpath %d %d M %d %d R %d %d R %d %d R closepath fill\n", x, y, 0, -h, w, 0,
            0, h);
    return;
}

void
ps_last_pt_off(void) {
    last_pt_line = 0;
    return;
}

void
ps_line(int32 xp1, int32 yp1, int32 xp2, int32 yp2) {
    last_pt_line = 1;
    if (no_break_line == 1) {
        fprintf(psfile, "%d %d M\n%d %d L\n", xp1, yp1, xp2, yp2);
        last_psx = xp2;
        last_psy = yp2;
        ps_check_lines();
        return;
    }
    if (xp1 == last_psx && yp1 == last_psy) {
        last_psx = xp2;
        last_psy = yp2;
        fprintf(psfile, "%d %d L\n", xp2, yp2);
        ps_check_lines();
        return;
    }
    if (xp2 == last_psx && yp2 == last_psy) {
        last_psx = xp1;
        last_psy = yp1;
        fprintf(psfile, "%d %d L\n", xp1, yp1);
        ps_check_lines();
        return;
    }
    fprintf(psfile, "%d %d M\n%d %d L\n", xp1, yp1, xp2, yp2);
    last_psx = xp2;
    last_psy = yp2;
    ps_check_lines();
    return;
}

void
ps_check_lines(void) {
    ps_lines++;
    if (ps_lines >= MAXPSLINE) {
        fprintf(psfile, "currentpoint stroke moveto\n");
        ps_lines = 0;
    }
    return;
}

void
ps_linetype(int32 linetype) {
    char *line = "ba0123456789c";

    fprintf(psfile, "LT%c\n", line[(linetype % 11) + 2]);
    ps_lines = 0;
    last_psx = -100000000;
    last_psy = -100000000;
    return;
}

void
ps_point(int32 x, int32 y) {
    int32 number = point_type;
    char *point = "PDABCTSKF";
    number %= POINT_TYPES;
    if (number < -1) {
        number = -1;
    }
    if (point_radius > 0) {
        number = 7;
    }
    fprintf(psfile, "%d %d %c\n", x, y, point[number + 1]);
    ps_lines = 0;
    last_pt_line = 0;
    return;
}

void
ps_write(char *str) {
    fprintf(psfile, "%s\n", str);
    return;
}

void
ps_fnt(int32 cf, int32 scale) {
    if (cf == 0) {
        fprintf(psfile, "/%s findfont %d scalefont setfont \n", ps_font, scale);
    } else {
        fprintf(psfile, "%d Symfnt\n", scale);
    }
    return;
}

void
ps_show(char *str, int32 type) {
    char ch;
    putc('(', psfile);
    ch = *str++;
    while (ch != '\0') {
        if ((ch == '(') || (ch == ')') || (ch == '\\')) {
            putc('\\', psfile);
        }
        putc(ch, psfile);
        ch = *str++;
    }
    if (type == 1) {
        fprintf(psfile, ") Lshow\n");
    } else {
        fprintf(psfile, ") show\n");
    }
    ps_lines = 0;
    return;
}

void
ps_abs(int32 x, int32 y) {
    fprintf(psfile, "%d %d moveto \n", x, y);
    return;
}

void
ps_rel(int32 x, int32 y) {
    fprintf(psfile, "%d %d rmoveto \n", x, y);
    return;
}

void
ps_special_put_text(int32 x, int32 y, char *str, int32 size) {
    int32 i = 0;
    int32 j = 0;
    int32 type = 1;
    int32 cf = 0;
    int32 n = (int32)strlen(str);
    int32 cy = 0;
    char tmp[256];
    char c;
    int32 sub;
    int32 sup;
    int32 pssz;
    static int32 sz[] = {8, 10, 14, 18, 24};
    fprintf(psfile, "0 0 0 setrgbcolor \n");
    ps_abs(x, y);
    pssz = sz[size]*PS_SC;
    sub = (int32)(.3*pssz);
    sup = (int32)(.6*pssz);
    // set the size here!
    ps_fnt(cf, pssz);
    while (i < n) {
        c = str[i];
        if (c == '\\') {
            i++;
            c = str[i];
            tmp[j] = 0;  // end the current buffer
            if (strlen(tmp) > 0) {
                ps_show(tmp, type);
                type = 0;
            }

            j = 0;
            if (c == '0') {
                cf = 0;
                ps_fnt(cf, pssz);
            }
            if (c == 'n') {
                ps_rel(0, -cy);
                cy = 0;
                pssz = PS_SC*sz[size];
                ps_fnt(cf, pssz);
            }
            if (c == 's') {
                cy = cy - sub;
                ps_rel(0, -sub);
                pssz = 3*PS_SC*sz[size] / 5;
                ps_fnt(cf, pssz);
            }
            if (c == 'S') {
                pssz = 3*PS_SC*sz[size] / 5;
                cy = cy + sup;
                ps_rel(0, sup);
                ps_fnt(cf, pssz);
            }
            if (c == '1') {
                cf = 1;
                ps_fnt(cf, pssz);
            }

            i++;
        } else {
            tmp[j] = c;
            j++;
            i++;
        }
    }
    tmp[j] = 0;
    if (strlen(tmp) > 0) {
        ps_show(tmp, type);
    }
    return;
}

void
ps_text(int32 x, int32 y, char *str) {
    char ch;
    fprintf(psfile, "0 0 0 setrgbcolor \n");
    fprintf(psfile, "/%s findfont %d ", ps_font, ps_fontsize*PS_SC);
    fprintf(psfile, "scalefont setfont\n");
    fprintf(psfile, "%d %d moveto\n", x, y);
    if (text_angle != 0) {
        fprintf(psfile, "currentpoint gsave translate %d rotate 0 0 moveto\n", text_angle*90);
    }
    putc('(', psfile);
    ch = *str++;
    while (ch != '\0') {
        if ((ch == '(') || (ch == ')') || (ch == '\\')) {
            putc('\\', psfile);
        }
        putc(ch, psfile);
        ch = *str++;
    }
    switch (text_justify) {
    case LEFT:
        fprintf(psfile, ") Lshow\n");
        break;
    case CENTER:
        fprintf(psfile, ") Cshow\n");
        break;
    case RIGHT:
        fprintf(psfile, ") Rshow\n");
        break;
    default:
        break;
    }
    if (text_angle != 0) {
        fprintf(psfile, "grestore\n");
    }
    ps_lines = 0;
}
