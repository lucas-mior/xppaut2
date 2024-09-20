#include "functions.h"
#include "integers.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define LEFT 0
#define RIGHT 2
#define CENTER 1
#define POINT_TYPES 8

char SVGLINETYPE;

extern int32 TextJustify;
extern int32 TextAngle;

extern int32 PointType;
extern int32 PointRadius;
FILE *svgfile;
extern int32 PltFmtFlag, PSColorFlag;
extern int32 PSLines;
extern int32 LastPtLine;
int32 cur_RGB[3];
extern int32 LastPSX, LastPSY;
extern int32 NoBreakLine;
int32 DOING_SVG_COLOR = 0;

int32 DO_MARKER = 0;
extern int32 DOING_DFIELD;

extern int32 Xup;

int32
svg_init(char *filename) {
    FILE *fp;
    char css[256];

    init_svg();

    LastPSX = -10000;
    LastPSY = -10000;

    if ((svgfile = fopen(filename, "w")) == NULL) {
        err_msg("Cannot open file ");
        return 0;
    }
    PltFmtFlag = SVGFMT;

    fprintf(svgfile, "<!-- Uncomment following when using your own custom "
                     "external stylesheet.-->\n");
    fprintf(svgfile, "<!--\n");
    fprintf(svgfile, "<?xml-stylesheet type=\"text/css\" "
                     "href=\"xppaut-stylesheet.css\" ?>\n");
    fprintf(svgfile, "-->\n");
    fprintf(svgfile, "<svg  xmlns=\"http://www.w3.org/2000/svg\"\n");
    fprintf(svgfile, "      xmlns:xlink=\"http://www.w3.org/1999/xlink\" "
                     "font-size=\"12pt\" width=\"640\" height=\"400\">\n");
    fprintf(svgfile, "\n\n      <defs>\n");
    fprintf(svgfile,
            "          <circle class=\"xpppointP\" id = \"xpppointP\"  "
            "r = \"1\"  stroke-width = \"1\"/>\n");
    fprintf(svgfile,
            "          <circle class=\"xppbead\" id = \"xppbead\"  r = "
            "\"1\"  stroke-width = \"1\"/>\n");
    fprintf(svgfile,
            "          <circle class=\"xpppointD\" id = \"xpppointD\"  "
            "r = \"1\"  stroke-width = \"1\"/>\n");
    fprintf(svgfile,
            "          <circle class=\"xpppointA\" id = \"xpppointA\"  "
            "r = \"1\"  stroke-width = \"1\"/>\n");
    fprintf(svgfile,
            "          <circle class=\"xpppointB\" id = \"xpppointB\"  "
            "r = \"1\"  stroke-width = \"1\"/>\n");
    fprintf(svgfile,
            "          <circle class=\"xpppointC\" id = \"xpppointC\"  "
            "r = \"1\"  stroke-width = \"1\"/>\n");
    fprintf(svgfile,
            "          <circle class=\"xpppointT\" id = \"xpppointT\"  "
            "r = \"1\"  stroke-width = \"1\"/>\n");
    fprintf(svgfile,
            "          <circle class=\"xpppointS\" id = \"xpppointS\"  "
            "r = \"1\"  stroke-width = \"1\"/>\n");
    fprintf(svgfile,
            "          <circle class=\"xpppointK\" id = \"xpppointK\"  "
            "r = \"3\"  stroke-width = \"0.75\"/>\n");
    fprintf(svgfile,
            "          <circle class=\"xpppointF\" id = \"xpppointF\"  "
            "r = \"2\"  stroke-width = \"0\"/>\n");
    fprintf(svgfile, "      </defs>\n\n");
    fprintf(svgfile, "\n\n");
    fprintf(svgfile, "      <!-- Comment out the following style block when "
                     "using your own custom external stylesheet.-->\n");
    fprintf(svgfile, "      <!-- As a starting point for your custom external "
                     "stylesheet, consider copying the style \n");
    fprintf(svgfile, "           information (between but not including CDATA "
                     "tags) to a file named xppaut-stylesheet.css \n");
    fprintf(svgfile, "       -->\n");
    fprintf(svgfile, "      <style type=\"text/css\">\n");
    fprintf(svgfile, "           <![CDATA[\n");
    fprintf(svgfile, "      \n");

    fprintf(svgfile, "                 circle.xpppointP {\n");
    fprintf(svgfile, "                    stroke-width: 1.0;\n");
    fprintf(svgfile, "                 }\n\n");
    fprintf(svgfile, "                 circle.xpppointD {\n");
    fprintf(svgfile, "                    stroke-width: 1.0;\n");
    fprintf(svgfile, "                 }\n\n");
    fprintf(svgfile, "                 circle.xpppointA {\n");
    fprintf(svgfile, "                    stroke-width: 1.0;\n");
    fprintf(svgfile, "                 }\n\n");
    fprintf(svgfile, "                 circle.xpppointB {\n");
    fprintf(svgfile, "                    stroke-width: 1.0;\n");
    fprintf(svgfile, "                 }\n\n");
    fprintf(svgfile, "                 circle.xpppointC {\n");
    fprintf(svgfile, "                    stroke-width: 1.0;\n");
    fprintf(svgfile, "                 }\n\n");
    fprintf(svgfile, "                 circle.xpppointT {\n");
    fprintf(svgfile, "                    stroke-width: 1.0;\n");
    fprintf(svgfile, "                 }\n\n");
    fprintf(svgfile, "                 circle.xpppointS {\n");
    fprintf(svgfile, "                    stroke-width: 1.0;\n");
    fprintf(svgfile, "                 }\n\n");
    fprintf(svgfile, "                 circle.xpppointK {\n");
    fprintf(svgfile, "                    stroke-width: 1.0;\n");
    fprintf(svgfile, "                 }\n\n");
    fprintf(svgfile, "                 circle.xpppointF {\n");
    fprintf(svgfile, "                    stroke-width: 1.0;\n");
    fprintf(svgfile, "                 }\n\n");
    fprintf(svgfile, "                 line.xppaxes {\n");
    fprintf(svgfile, "                    stroke: #000000;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppboxaxes {\n");
    fprintf(svgfile, "                    stroke: #000000;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppdfield {\n");
    fprintf(svgfile, "                    stroke: #000000;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xpplineb {\n");
    fprintf(svgfile, "                    stroke: #000000;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xpplinea {\n");
    fprintf(svgfile, "                    stroke-dasharray: 2,8;\n");
    fprintf(svgfile, "                    stroke-width: 2;\n");
    fprintf(svgfile, "                    stroke: #000000;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppline0 {\n");
    fprintf(svgfile, "                    stroke: #000000;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppline1 {\n");
    fprintf(svgfile, "                    stroke-width: 1;\n");
    fprintf(svgfile, "                    stroke: #FF0000;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppline2 {\n");
    fprintf(svgfile, "                    stroke: #F06400;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppline3 {\n");
    fprintf(svgfile, "                    stroke: #FFA500;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppline4 {\n");
    fprintf(svgfile, "                    stroke: #FFCD00;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppline5 {\n");
    fprintf(svgfile, "                    stroke: #C8C800;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppline6 {\n");
    fprintf(svgfile, "                    stroke: #00FF00;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppline7 {\n");
    fprintf(svgfile, "                    stroke: #32CD32;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppline8 {\n");
    fprintf(svgfile, "                    stroke: #00C8C8;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xppline9 {\n");
    fprintf(svgfile, "                    stroke: #0000FF;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 line.xpplinec {\n");
    fprintf(svgfile, "                    stroke: #000000;\n");
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 \n");
    fprintf(svgfile, "                 text.xpptext {\n");
    fprintf(svgfile, "                    font-family: sans-serif;\n");
    fprintf(svgfile, "                    font-size  : 1em;\n");
    fprintf(svgfile, "                    stroke	: #000000;\n");
    fprintf(svgfile,
            "                    fill	: #000000;\n"); /*Need to support more
                                                           than 1 vertical
                                                           centering tactic!*/
    fprintf(
        svgfile,
        "                    baseline-shift:-33%%;\n"); /*Supported in Inkscape
                                                           v0.48.2, but not in
                                                           Firefox v13*/
    fprintf(
        svgfile,
        "                    dominant-baseline: central;\n"); /*Supported in
                                                                 Firefox v13,
                                                                 but not in
                                                                 Inkscape
                                                                 v0.48.2*/
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 text.xppyaxislabelh {\n");
    fprintf(svgfile, "                    font-family: sans-serif;\n");
    fprintf(svgfile, "                    font-size  : 1em;\n");
    fprintf(svgfile, "                    stroke	: none;\n");
    fprintf(svgfile,
            "                    fill	: none;\n"); /*Need to support more than
                                                        1 vertical centering
                                                        tactic!*/
    fprintf(
        svgfile,
        "                    baseline-shift:-33%%;\n"); /*Supported in Inkscape
                                                           v0.48.2, but not in
                                                           Firefox v13*/
    fprintf(
        svgfile,
        "                    dominant-baseline: central;\n"); /*Supported in
                                                                 Firefox v13,
                                                                 but not in
                                                                 Inkscape
                                                                 v0.48.2*/
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 text.xppyaxislabelv {\n");
    fprintf(svgfile, "                    font-family: sans-serif;\n");
    fprintf(svgfile, "                    font-size  : 1em;\n");
    fprintf(svgfile, "                    stroke	: #000000;\n");
    fprintf(svgfile,
            "                    fill	: #000000;\n"); /*Need to support more
                                                           than 1 vertical
                                                           centering tactic!*/
    fprintf(
        svgfile,
        "                    baseline-shift:-33%%;\n"); /*Supported in Inkscape
                                                           v0.48.2, but not in
                                                           Firefox v13*/
    fprintf(
        svgfile,
        "                    dominant-baseline: central;\n"); /*Supported in
                                                                 Firefox v13,
                                                                 but not in
                                                                 Inkscape
                                                                 v0.48.2*/
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 text.xppaxestext {\n");
    fprintf(svgfile, "                    font-family: sans-serif;\n");
    fprintf(svgfile, "                    font-size  : 1em;\n");
    fprintf(svgfile, "                    stroke	: #000000;\n");
    fprintf(svgfile,
            "                    fill	: #000000;\n"); /*Need to support more
                                                           than 1 vertical
                                                           centering tactic!*/
    fprintf(
        svgfile,
        "                    baseline-shift:-33%%;\n"); /*Supported in Inkscape
                                                           v0.48.2, but not in
                                                           Firefox v13*/
    fprintf(
        svgfile,
        "                    dominant-baseline: central;\n"); /*Supported in
                                                                 Firefox v13,
                                                                 but not in
                                                                 Inkscape
                                                                 v0.48.2*/
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 \n");
    fprintf(svgfile, "                 text.xpptext0 {\n");
    fprintf(svgfile, "                    font-family: sans-serif;\n");
    fprintf(svgfile, "                    font-size  : 0.5em;\n");
    fprintf(svgfile, "                    stroke	: #000000;\n");
    fprintf(svgfile,
            "                    fill	: #000000;\n"); /*Need to support more
                                                           than 1 vertical
                                                           centering tactic!*/
    fprintf(
        svgfile,
        "                    baseline-shift:-33%%;\n"); /*Supported in Inkscape
                                                           v0.48.2, but not in
                                                           Firefox v13*/
    fprintf(
        svgfile,
        "                    dominant-baseline: central;\n"); /*Supported in
                                                                 Firefox v13,
                                                                 but not in
                                                                 Inkscape
                                                                 v0.48.2*/
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 \n");
    fprintf(svgfile, "                 text.xpptext1 {\n");
    fprintf(svgfile, "                    font-family: sans-serif;\n");
    fprintf(svgfile, "                    font-size  : 0.75em;\n");
    fprintf(svgfile, "                    stroke	: #000000;\n");
    fprintf(svgfile,
            "                    fill	: #000000;\n"); /*Need to support more
                                                           than 1 vertical
                                                           centering tactic!*/
    fprintf(
        svgfile,
        "                    baseline-shift:-33%%;\n"); /*Supported in Inkscape
                                                           v0.48.2, but not in
                                                           Firefox v13*/
    fprintf(
        svgfile,
        "                    dominant-baseline: central;\n"); /*Supported in
                                                                 Firefox v13,
                                                                 but not in
                                                                 Inkscape
                                                                 v0.48.2*/
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 \n");
    fprintf(svgfile, "                 text.xpptext2 {\n");
    fprintf(svgfile, "                    font-family: sans-serif;\n");
    fprintf(svgfile, "                    font-size  : 1em;\n");
    fprintf(svgfile, "                    stroke	: #000000;\n");
    fprintf(svgfile,
            "                    fill	: #000000;\n"); /*Need to support more
                                                           than 1 vertical
                                                           centering tactic!*/
    fprintf(
        svgfile,
        "                    baseline-shift:-33%%;\n"); /*Supported in Inkscape
                                                           v0.48.2, but not in
                                                           Firefox v13*/
    fprintf(
        svgfile,
        "                    dominant-baseline: central;\n"); /*Supported in
                                                                 Firefox v13,
                                                                 but not in
                                                                 Inkscape
                                                                 v0.48.2*/
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 \n");
    fprintf(svgfile, "                 text.xpptext3 {\n");
    fprintf(svgfile, "                    font-family: sans-serif;\n");
    fprintf(svgfile, "                    font-size  : 1.25em;\n");
    fprintf(svgfile, "                    stroke	: #000000;\n");
    fprintf(svgfile,
            "                    fill	: #000000;\n"); /*Need to support more
                                                           than 1 vertical
                                                           centering tactic!*/
    fprintf(
        svgfile,
        "                    baseline-shift:-33%%;\n"); /*Supported in Inkscape
                                                           v0.48.2, but not in
                                                           Firefox v13*/
    fprintf(
        svgfile,
        "                    dominant-baseline: central;\n"); /*Supported in
                                                                 Firefox v13,
                                                                 but not in
                                                                 Inkscape
                                                                 v0.48.2*/
    fprintf(svgfile, "                 }\n");
    fprintf(svgfile, "                 \n");
    fprintf(svgfile, "                 text.xpptext4 {\n");
    fprintf(svgfile, "                    font-family: sans-serif;\n");
    fprintf(svgfile, "                    font-size  : 1.5em;\n");
    fprintf(svgfile, "                    stroke	: #000000;\n");
    fprintf(svgfile,
            "                    fill	: #000000;\n"); /*Need to support more
                                                           than 1 vertical
                                                           centering tactic!*/
    fprintf(
        svgfile,
        "                    baseline-shift:-33%%;\n"); /*Supported in Inkscape
                                                           v0.48.2, but not in
                                                           Firefox v13*/
    fprintf(
        svgfile,
        "                    dominant-baseline: central;\n"); /*Supported in
                                                                 Firefox v13,
                                                                 but not in
                                                                 Inkscape
                                                                 v0.48.2*/
    fprintf(svgfile, "                 }\n");

    snprintf(css, sizeof(css), "%s/xppaut-stylesheet.css", getenv("HOME"));
    fp = fopen(css, "r");
    if (fp != NULL) {
        char bob[256];
        plintf("Styling svg image according to %s\n", css);
        while (!feof(fp)) {
            bob[0] = '\0';
            fgets(bob, 255, fp);
            fprintf(svgfile, "%s", bob);
        }

        fclose(fp);
    }

    fprintf(svgfile, "           ]]>\n");
    fprintf(svgfile, "      </style>\n\n");

    return 1;
}

void
svg_write(char *str) {
    fprintf(svgfile, "%s\n", str);
    return;
}

void
svg_stroke(void) {
}

void
svg_do_color(int32 color) {
    int32 r, g, b;

    if (PltFmtFlag == SCRNFMT)
        return;
    if (PltFmtFlag == PSFMT)
        return;
    if (PSColorFlag == 0)
        return;
    get_svg_color(color, &r, &g, &b);
    cur_RGB[0] = r;
    cur_RGB[1] = g;
    cur_RGB[2] = b;

    DOING_SVG_COLOR = 1;
    return;
}

void
svg_end(void) {
    svg_write("</svg>");
    fclose(svgfile);
    PltFmtFlag = SCRNFMT;
    DOING_SVG_COLOR = 0;
    if (Xup)
        init_x11();
    return;
}

void
svg_bead(void) {
    DO_MARKER = 1;
    return;
}

void
svg_frect(int32 x, int32 y, int32 w, int32 h) {
    double gray = 0;
    if (DOING_SVG_COLOR) {
        fprintf(svgfile,
                "      <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" "
                "style=\"stroke:rgb(%d,%d,%d);fill:rgb(%d,%d,%d);\"/>",
                x, y, w, h, cur_RGB[0], cur_RGB[1], cur_RGB[2], cur_RGB[0],
                cur_RGB[1], cur_RGB[2]);
    } else {
        gray = (0.299*cur_RGB[0] + 0.587*cur_RGB[1] + 0.114*cur_RGB[2]);
        fprintf(svgfile,
                "      <rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" "
                "style=\"stroke:rgb(%d,%d,%d);fill:rgb(%d,%d,%d);\"/>",
                x, y, w, h, (int32)gray, (int32)gray, (int32)gray, (int32)gray,
                (int32)gray, (int32)gray);
    }
    return;
}

void
svg_last_pt_off(void) {
    LastPtLine = 0;
    return;
}

void
svg_line(int32 xp1, int32 yp1, int32 xp2, int32 yp2) {
    if (DOING_SVG_COLOR) {
        if (DOING_AXES) {
            if (DOING_BOX_AXES) {
                fprintf(
                    svgfile,
                    "      <line class=\"xppboxaxes\"  x1=\"%d\"  y1=\"%d\" "
                    "x2=\"%d\"   y2=\"%d\" style=\"stroke:rgb(%d,%d,%d);\"/>\n",
                    xp1, yp1, xp2, yp2, cur_RGB[0], cur_RGB[1], cur_RGB[2]);
            } else {
                fprintf(svgfile,
                        "      <line class=\"xppaxes\"  x1=\"%d\"  y1=\"%d\" "
                        "x2=\"%d\" "
                        "  y2=\"%d\" style=\"stroke:rgb(%d,%d,%d);\"/>\n",
                        xp1, yp1, xp2, yp2, cur_RGB[0], cur_RGB[1], cur_RGB[2]);
            }
        } else {
            if (DOING_DFIELD) {
                if (DO_MARKER) {
                    fprintf(svgfile, "<g>\n");
                }

                fprintf(
                    svgfile,
                    "      <line class=\"xppdfield\"  x1=\"%d\"  y1=\"%d\" "
                    "x2=\"%d\"   y2=\"%d\" style=\"stroke:rgb(%d,%d,%d);\"/>\n",
                    xp1, yp1, xp2, yp2, cur_RGB[0], cur_RGB[1], cur_RGB[2]);
                if (DO_MARKER) {
                    fprintf(svgfile,
                            "      <use xlink:href = \"#xppbead\" x=\"%d\" "
                            "y=\"%d\" "
                            "style=\"stroke:rgb(%d,%d,%d); "
                            "fill:rgb(%d,%d,%d)\"/>\n",
                            xp2, yp2, cur_RGB[0], cur_RGB[1], cur_RGB[2],
                            cur_RGB[0], cur_RGB[1], cur_RGB[2]);
                    fprintf(svgfile, "</g>\n");
                }
            } else {
                fprintf(
                    svgfile,
                    "      <line class=\"xppline%c\"  x1=\"%d\"  y1=\"%d\" "
                    "x2=\"%d\"   y2=\"%d\" style=\"stroke:rgb(%d,%d,%d);\"/>\n",
                    SVGLINETYPE, xp1, yp1, xp2, yp2, cur_RGB[0], cur_RGB[1],
                    cur_RGB[2]);
            }
        }

    } else {
        if (DOING_AXES) {
            if (DOING_BOX_AXES) {
                fprintf(
                    svgfile,
                    "      <line class=\"xppboxaxes\"  x1=\"%d\"  y1=\"%d\" "
                    "x2=\"%d\"   y2=\"%d\" />\n",
                    xp1, yp1, xp2, yp2);
            } else {
                fprintf(svgfile,
                        "      <line class=\"xppaxes\"  x1=\"%d\"  y1=\"%d\" "
                        "x2=\"%d\" "
                        "  y2=\"%d\" />\n",
                        xp1, yp1, xp2, yp2);
            }
        } else {

            if (DOING_DFIELD) {
                if (DO_MARKER) {
                    fprintf(svgfile, "<g>\n");
                }
                fprintf(svgfile,
                        "      <line class=\"xppdfield\"  x1=\"%d\"  y1=\"%d\" "
                        "x2=\"%d\"   y2=\"%d\" />\n",
                        xp1, yp1, xp2, yp2);
                if (DO_MARKER) {
                    fprintf(svgfile,
                            "      <use xlink:href = \"#xppbead\" x=\"%d\" "
                            "y=\"%d\" />\n",
                            xp2, yp2);
                    fprintf(svgfile, "</g>\n");
                }

            } else {
                fprintf(svgfile,
                        "      <line class=\"xppline%c\"  x1=\"%d\"  y1=\"%d\" "
                        "x2=\"%d\"   y2=\"%d\" />\n",
                        SVGLINETYPE, xp1, yp1, xp2, yp2);
            }
        }
    }

    LastPSX = xp2;
    LastPSY = yp2;

    DOING_SVG_COLOR = 0;
    DO_MARKER = 0;
    return;
}

void
chk_svg_lines(void) {
    /*PSLines++;
    if(PSLines>=MAXPSLINE){
      fprintf(psfile,"currentpoint stroke moveto\n");
      PSLines=0;
    } */
    return;
}

void
svg_linetype(int32 linetype) {
    char *line = "ba0123456789c";

    SVGLINETYPE = line[(linetype % 11) + 2];

    PSLines = 0;
    /* LastPSX=-100000000;
     LastPSY=-100000000;
     */
    return;
}

void
svg_point(int32 x, int32 y) {
    char svgcol[8];
    char svgfill[8];
    int32 number = PointType;
    char *point = "PDABCTSKF";

    snprintf(svgfill, sizeof(svgfill), "none");
    svgcol[0] = '\0';

    number %= POINT_TYPES;
    if (number < -1)
        number = -1;
    if (PointRadius > 0)
        number = 7;

    if (number == 7) {
        snprintf(svgcol, sizeof(svgcol), "00FF00");
        snprintf(svgfill, sizeof(svgfill), "#00FF00");
    } else if (number == 6) {
        snprintf(svgcol, sizeof(svgcol), "0000FF");
    } else {
        snprintf(svgcol, sizeof(svgcol), "000000");
        snprintf(svgfill, sizeof(svgfill), "#000000");
    }

    if (DOING_SVG_COLOR) {
        fprintf(svgfile,
                "      <use xlink:href = \"#xpppoint%c\" x=\"%d\" y=\"%d\" "
                "style=\"stroke:rgb(%d,%d,%d); fill:rgb(%d,%d,%d)\"/>\n",
                point[number + 1], x, y, cur_RGB[0], cur_RGB[1], cur_RGB[2],
                cur_RGB[0], cur_RGB[1], cur_RGB[2]);
    } else {
        fprintf(svgfile,
                "      <use xlink:href = \"#xpppoint%c\" x=\"%d\" y=\"%d\" "
                "style=\"stroke:#%s; fill:%s\"/>\n",
                point[number + 1], x, y, svgcol, svgfill);
    }

    PSLines = 0;
    LastPtLine = 0;
    DOING_SVG_COLOR = 0;
    return;
}

void
special_put_text_svg(int32 x, int32 y, char *str, int32 size) {
    /*int32 i=0,j=0,type=1;
    int32 cf=0;

    int32 n=strlen(str);
    int32 cy=0;
    char tmp[256],c;
    int32 sub,sup,pssz;
    static int32 sz[]={8,10,14,18,24};
    fprintf(psfile, "0 0 0 setrgbcolor \n");
    ps_abs(x,y);
    pssz=sz[size]*PS_SC;
    sub=.3*pssz;
    sup=.6*pssz;
    */
    /* set the size here! */

    /*ps_fnt(cf,pssz);
    while(i<n){
      c=str[i];
      if(c=='\\'){
        i++;
        c=str[i];
        tmp[j]=0;*/ /* end the current buffer */
    /*if(strlen(tmp)>0){
      ps_show(tmp,type);
      type=0;
    }

    j=0;
    if(c=='0'){
      cf=0;
      ps_fnt(cf,pssz);
    }
    if(c=='n'){

      ps_rel(0,-cy);
      cy=0;
      pssz=PS_SC*sz[size];
      ps_fnt(cf,pssz);
    }
    if(c=='s'){

      cy=cy-sub;
      ps_rel(0,-sub);
      pssz=3*PS_SC*sz[size]/5;
      ps_fnt(cf,pssz);
    }
    if(c=='S'){
      pssz=3*PS_SC*sz[size]/5;
      cy=cy+sup;
      ps_rel(0,sup);
      ps_fnt(cf,pssz);
    }
    if(c=='1'){

      cf=1;
      ps_fnt(cf,pssz);
    }

    i++;
  }
  else {
    tmp[j]=c;
    j++;
    i++;
  }
  }
  tmp[j]=0;
  if(strlen(tmp)>0)
  ps_show(tmp,type);
  */

    char anchor[7];

    switch (TextJustify) {
    case LEFT:
        snprintf(anchor, sizeof(anchor), "start");
        break;
    case CENTER:
        snprintf(anchor, sizeof(anchor), "middle");
        break;
    case RIGHT:
        snprintf(anchor, sizeof(anchor), "end");
        break;
    default:
        snprintf(anchor, sizeof(anchor), "start");
        break;
    }

    fprintf(svgfile,
            "\n      <text class=\"xpptext%d\" text-anchor=\"%s\" x=\"%d\"  "
            "y=\"%d\"\n",
            size, anchor, x, y);

    fprintf(svgfile, "      >%s</text>\n", str);
    return;
}

void
svg_text(int32 x, int32 y, char *str) {
    char anchor[7];

    switch (TextJustify) {
    case LEFT:
        snprintf(anchor, sizeof(anchor), "start");
        break;
    case CENTER:
        snprintf(anchor, sizeof(anchor), "middle");
        break;
    case RIGHT:
        snprintf(anchor, sizeof(anchor), "end");
        break;
    default:
        snprintf(anchor, sizeof(anchor), "start");
        break;
    }

    if (DOING_AXES) {
        fprintf(
            svgfile,
            "\n      <text class=\"xppaxestext\" text-anchor=\"%s\" x=\"%d\"  "
            "y=\"%d\"\n",
            anchor, x, y);

    } else {
        fprintf(svgfile,
                "\n      <text class=\"xpptext\" text-anchor=\"%s\" x=\"%d\"  "
                "y=\"%d\"\n",
                anchor, x, y);
    }

    fprintf(svgfile, "      >%s</text>\n", str);

    /* char ch;
     fprintf(psfile, "0 0 0 setrgbcolor \n");
     fprintf(psfile,"/%s findfont %d ",PS_FONT,PS_FONTSIZE*PS_SC);
    fprintf(psfile,"scalefont setfont\n");
    fprintf(psfile,"%d %d moveto\n",x,y);
    if (TextAngle != 0)
      fprintf(psfile,"currentpoint gsave translate %d rotate 0 0 moveto\n"
              ,TextAngle*90);
    putc('(',psfile);
    ch = *str++;
    while(ch!='\0') {
      if ( (ch=='(') || (ch==')') || (ch=='\\') )
        putc('\\',psfile);
      putc(ch,psfile);
      ch = *str++;
    }
    switch(TextJustify) {
    case LEFT : fprintf(psfile,") Lshow\n");
      break;
    case CENTER : fprintf(psfile,") Cshow\n");
      break;
    case RIGHT : fprintf(psfile,") Rshow\n");
      break;
    }
    if (TextAngle != 0)
      fprintf(psfile,"grestore\n");
    PSLines=0;
    */
}
