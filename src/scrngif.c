#include "functions.h"
#include "integers.h"
#include "xmalloc.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xproto.h>

#define MAKE_ONE_GIF 2
#define GET_GLOBAL_CMAP 1
#define FIRST_ANI_GIF 3
#define NEXT_ANI_GIF 4
#define BLOKLEN 255
#define BUFLEN 1000
#define TERMIN 'T'
#define LOOKUP 'L'
#define SEARCH 'S'
#define NUMBER_OF_ARRAYS 20
/* defines the amount of memory set aside in the encoding for the
 * LOOKUP type nodes; for a 256 color GIF, the number of LOOKUP
 * nodes will be <= NUMBER_OF_ARRAYS, for a 128 color GIF the number of
 * LOOKUP nodes will be <= 2*NUMBER_OF_ARRAYS, etc.  */
#define GifPutShort(i, fout)                                                                       \
    do {                                                                                           \
        fputc(i & 0xff, fout);                                                                     \
        fputc(i >> 8, fout);                                                                       \
    } while (0)

static uint32 debugFlag;
static int32 UseGlobalMap = 0;
static int32 gif_frame_delay = 5;
static int32 gif_frame_loop = 1000;
static int32 chainlen = 0;
static int32 maxchainlen = 0;
static int32 nodecount = 0;
static int32 lookuptypes = 0;
static int16 need = 8;
static GifTree *empty[256];
static GifTree gif_root = {LOOKUP, 0, 0, empty, NULL, NULL};
static GifTree *topNode;
static GifTree *baseNode;
static GifTree **node_array;
static GifTree **last_array;

static void scrngif_clear_tree(int32 cc, GifTree *root);
static int32 scrngif_encode(FILE *fout, uchar *pixels, int32 depth, int32 siz);
static void scrngif_make_gif(uchar *pixels, int32 cols, int32 rows, FILE *dst);
static void scrngif_loop(FILE *fout, uint32 repeats);
static void scrngif_write_global_header(int32 cols, int32 rows, FILE *dst);
static void scrngif_stuff(Window win, FILE *fp, int32 task);
static int32 scrngif_make_local_map(uchar *pixels, uchar *ppm, int32 h, int32 w);
static int32 scrngif_use_global_map(uchar *pixels, uchar *ppm, int32 h, int32 w);
static void scrngif_local_to_global(void);
static int32 scrngif_ppm_to_pix(uchar r, uchar g, uchar b, int32 *n);
static uchar *scrngif_add_code_to_buffer(int32, int16, uchar *);
static void scrngif_write_local_header(int32 cols, int32 rows, FILE *fout, int32 colflag,
                                       int32 delay2);

typedef struct GifCol {
    uchar r;
    uchar g;
    uchar b;
} GifCol;

static GifCol gifcol[256];
static GifCol gifGcol[256];
static int32 nglobal_colors = 0;

void
scrngif_set_global_map(int32 flag) {
    if (nglobal_colors == 0) {  // Cant use it if it aint there
        UseGlobalMap = 0;
        return;
    }
    UseGlobalMap = flag;
    return;
}

int32
scrngif_ppm_to_pix(uchar r, uchar g, uchar b, int32 *n) {
    int32 nc = *n;

    if (UseGlobalMap == 1) {
        for (int32 i = 0; i < nglobal_colors; i++) {
            if (r == gifGcol[i].r && g == gifGcol[i].g && b == gifGcol[i].b) {
                return i;
            }
        }

        return -1;
    }
    for (int32 i = 0; i < nc; i++) {
        if (r == gifcol[i].r && g == gifcol[i].g && b == gifcol[i].b) {
            return i;
        }
    }
    if (nc > 255) {
        ggets_plintf("Too many colors \n");
        return -1;
    }
    gifcol[nc].r = r;
    gifcol[nc].g = g;
    gifcol[nc].b = b;
    nc++;
    *n = nc;
    return nc - 1;
}

void
scrngif_end_ani_gif(FILE *fp) {
    fputc(';', fp);
    return;
}

void
scrngif_add_ani_gif(Window win, FILE *fp, int32 count) {
    ggets_plintf("Frame %d \n", count);
    if (count == 0) {
        scrngif_stuff(win, fp, FIRST_ANI_GIF);
    } else {
        scrngif_stuff(win, fp, NEXT_ANI_GIF);
    }
    return;
}

void
scrngif_screen_to_gif(Window win, FILE *fp) {
    scrngif_stuff(win, fp, MAKE_ONE_GIF);
    return;
}

void
scrngif_get_global_colormap(Window win) {
    FILE *junk = NULL;
    scrngif_stuff(win, junk, GET_GLOBAL_CMAP);
    return;
}

void
scrngif_local_to_global(void) {
    for (int32 i = 0; i < 256; i++) {
        gifcol[i].r = gifGcol[i].r;
        gifcol[i].g = gifGcol[i].g;
        gifcol[i].b = gifGcol[i].b;
    }
    return;
}

int32
scrngif_use_global_map(uchar *pixels, uchar *ppm, int32 h, int32 w) {
    uchar r;
    uchar g;
    uchar b;
    int32 k = 0;
    int32 l = 0;
    int32 pix;
    int32 nc;
    for (int32 i = 0; i < h; i++) {
        for (int32 j = 0; j < w; j++) {
            r = ppm[k];
            g = ppm[k + 1];
            b = ppm[k + 2];
            pix = scrngif_ppm_to_pix(r, g, b, &nc);
            if (pix < 0) {
                return 0;
            }
            pixels[l] = (uchar)pix;
            k += 3;
            l++;
        }
    }
    return 1;
}

int32
scrngif_make_local_map(uchar *pixels, uchar *ppm, int32 h, int32 w) {
    uchar r;
    uchar g;
    uchar b;
    int32 k = 0;
    int32 l = 0;
    int32 pix;
    int32 ncol = 0;
    for (int32 i = 0; i < h; i++) {
        for (int32 j = 0; j < w; j++) {
            r = ppm[k];
            g = ppm[k + 1];
            b = ppm[k + 2];
            k += 3;
            pix = scrngif_ppm_to_pix(r, g, b, &ncol);
            if (pix < 0) {
                pix = 255;
            }
            pixels[l] = (uchar)pix;
            l++;
        }
    }
    ggets_plintf("Got %d colors\n", ncol);
    for (int32 i = ncol; i < 256; i++) {
        gifcol[i].r = 255;
        gifcol[i].g = 255;
        gifcol[i].b = 255;
    }
    return ncol;
}

void
scrngif_stuff(Window win, FILE *fp, int32 task) {
    Window root;
    uint32 h;
    uint32 w;
    uint32 bw;
    uint32 d;
    int32 x0;
    int32 y0;
    uchar *ppm;

    uchar *pixels;
    int32 ncol = 0;

    int32 ok;

    XGetGeometry(display, win, &root, &x0, &y0, &w, &h, &bw, &d);
    ppm = xmalloc(w*h*3);
    pixels = xmalloc(h*w);

    ani_get_ppm_bits(win, (int32 *)&w, (int32 *)&h, ppm);
    switch (task) {
    case GET_GLOBAL_CMAP:
        ncol = scrngif_make_local_map(pixels, ppm, (int32)h, (int32)w);
        for (int32 i = 0; i < 256; i++) {
            gifGcol[i].r = gifcol[i].r;
            gifGcol[i].g = gifcol[i].g;
            gifGcol[i].b = gifcol[i].b;
        }
        nglobal_colors = ncol;

        break;
    case MAKE_ONE_GIF:  // don't need global map!
        ncol = scrngif_make_local_map(pixels, ppm, (int32)h, (int32)w);
        scrngif_make_gif(pixels, (int32)w, (int32)h, fp);
        break;
    case FIRST_ANI_GIF:
        if (UseGlobalMap) {
            ok = scrngif_use_global_map(pixels, ppm, (int32)h, (int32)w);
            if (ok == 1) {
                scrngif_local_to_global();
                scrngif_write_global_header((int32)w, (int32)h, fp);
                scrngif_write_local_header((int32)w, (int32)h, fp, 0, gif_frame_delay);
                scrngif_encode(fp, pixels, 8, (int32)(w*h));
            } else  // first map cant be encoded
            {
                UseGlobalMap = 0;
                scrngif_local_to_global();
                scrngif_write_global_header((int32)w, (int32)h,
                                            fp);  // write global header
                scrngif_make_local_map(pixels, ppm, (int32)h, (int32)w);
                scrngif_write_local_header((int32)w, (int32)h, fp, 1, gif_frame_delay);
                scrngif_encode(fp, pixels, 8, (int32)(w*h));
                UseGlobalMap = 1;
            }
        } else {
            scrngif_make_local_map(pixels, ppm, (int32)h, (int32)w);
            scrngif_write_global_header((int32)w, (int32)h, fp);
            scrngif_write_local_header((int32)w, (int32)h, fp, 0, gif_frame_delay);
            scrngif_encode(fp, pixels, 8, (int32)(w*h));
        }
        break;
    case NEXT_ANI_GIF:
        if (UseGlobalMap) {
            ok = scrngif_use_global_map(pixels, ppm, (int32)h, (int32)w);
            if (ok == 1) {
                scrngif_write_local_header((int32)w, (int32)h, fp, 0, gif_frame_delay);
                scrngif_encode(fp, pixels, 8, (int32)(w*h));
            } else {
                UseGlobalMap = 0;
                scrngif_make_local_map(pixels, ppm, (int32)h, (int32)w);
                scrngif_write_local_header((int32)w, (int32)h, fp, 1, gif_frame_delay);
                scrngif_encode(fp, pixels, 8, (int32)(w*h));
                UseGlobalMap = 1;
            }
        } else {
            scrngif_make_local_map(pixels, ppm, (int32)h, (int32)w);
            scrngif_write_local_header((int32)w, (int32)h, fp, 1, gif_frame_delay);
            scrngif_encode(fp, pixels, 8, (int32)(w*h));
        }
        break;
    default:
        break;
    }
    free(pixels);
    free(ppm);
    return;
}

void
scrngif_write_global_header(int32 cols, int32 rows, FILE *dst) {
    uchar *pos;
    uchar *buffer;

    buffer = xmalloc((BUFLEN + 1)*sizeof(*buffer));
    buffer += 1;

    pos = buffer;

    *pos++ = 'G';
    *pos++ = 'I';
    *pos++ = 'F';
    *pos++ = '8';
    *pos++ = '9';
    *pos++ = 'a';

    *pos++ = 0xff & cols;
    *pos++ = (0xff00 & cols) / 0x100;
    *pos++ = 0xff & rows;
    *pos++ = (0xff00 & rows) / 0x100;
    *pos++ = 0x87;
    *pos++ = 0xff;
    *pos++ = 0x0;

    for (int32 i = 0; i < 256; i++) {
        *pos++ = 0xff & gifcol[i].r;
        *pos++ = 0xff & gifcol[i].g;
        *pos++ = 0xff & gifcol[i].b;
    }
    fwrite(buffer, (ulong)(pos - buffer), 1, dst);
    free(buffer - 1);
    scrngif_loop(dst, (uint32)gif_frame_loop);
    return;
}

void
scrngif_loop(FILE *fout, uint32 repeats) {
    fputc(0x21, fout);
    fputc(0xFF, fout);
    fputc(0x0B, fout);
    fputs("NETSCAPE2.0", fout);

    fputc(0x03, fout);
    fputc(0x01, fout);
    GifPutShort(repeats, fout);  // repeat count

    fputc(0x00, fout);  // terminator
    return;
}

void
scrngif_write_local_header(int32 cols, int32 rows, FILE *fout, int32 colflag, int32 delay2) {
    fputc(0x21, fout);
    fputc(0xF9, fout);
    fputc(0x04, fout);
    fputc(0x80, fout);  // flag ???
    GifPutShort(delay2, fout);
    fputc(0x00, fout);
    fputc(0x00, fout);
    fputc(',', fout);  // image separator
    GifPutShort(0, fout);
    GifPutShort(0, fout);
    GifPutShort(cols, fout);
    GifPutShort(rows, fout);
    if (colflag) {
        fputc(0x87, fout);
    } else {
        fputc(0x07, fout);
    }
    if (colflag) {
        for (int32 i = 0; i < 256; i++) {
            fputc(0xff & gifcol[i].r, fout);
            fputc(0xff & gifcol[i].g, fout);
            fputc(0xff & gifcol[i].b, fout);
        }
    }
    return;
}

void
scrngif_make_gif(uchar *pixels, int32 cols, int32 rows, FILE *dst) {
    int32 depth = 8;

    uchar *pos;
    uchar *buffer;

    buffer = xmalloc((BUFLEN + 1)*sizeof(*buffer));
    buffer += 1;

    pos = buffer;

    *pos++ = 'G';
    *pos++ = 'I';
    *pos++ = 'F';
    *pos++ = '8';
    *pos++ = '7';
    *pos++ = 'a';

    *pos++ = 0xff & cols;
    *pos++ = (0xff00 & cols) / 0x100;
    *pos++ = 0xff & rows;
    *pos++ = (0xff00 & rows) / 0x100;
    *pos++ = 0xf0 | (0x7 & (depth - 1));
    *pos++ = 0xff;
    *pos++ = 0x0;

    for (int32 i = 0; i < 256; i++) {
        *pos++ = 0xff & gifcol[i].r;
        *pos++ = 0xff & gifcol[i].g;
        *pos++ = 0xff & gifcol[i].b;
    }
    *pos++ = 0x2c;
    *pos++ = 0x00;
    *pos++ = 0x00;
    *pos++ = 0x00;
    *pos++ = 0x00;
    *pos++ = 0xff & cols;
    *pos++ = (0xff00 & cols) / 0x100;
    *pos++ = 0xff & rows;
    *pos++ = (0xff00 & rows) / 0x100;
    *pos++ = 0x7 & (depth - 1);

    fwrite(buffer, (usize)(pos - buffer), 1, dst);

    // header info done

    scrngif_encode(dst, pixels, depth, rows*cols);
    fputc(';', dst);
    free(buffer - 1);
    return;
}

int32
scrngif_encode(FILE *fout, uchar *pixels, int32 depth, int32 siz) {
    GifTree *first = &gif_root, *newNode, *curNode;
    uchar *end;
    int32 cc;
    int32 eoi;
    int32 next;
    int32 tel = 0;
    int16 cLength;

    uchar *pos;
    uchar *buffer;

    empty[0] = NULL;
    need = 8;

    node_array = empty;
    memmove(++node_array, empty, 255*sizeof(GifTree **));
    if ((buffer = xmalloc((BUFLEN + 1)*sizeof(uchar))) == NULL) {
        return 0;
    }
    buffer++;

    pos = buffer;
    buffer[0] = 0x0;

    cc = (depth == 1) ? 0x4 : 1 << depth;
    fputc((depth == 1) ? 2 : depth, fout);
    eoi = cc + 1;
    next = cc + 2;

    cLength = (int16)((depth == 1) ? 3 : depth + 1);

    if ((topNode = baseNode = xmalloc(sizeof(GifTree)*4094)) == NULL) {
        return 0;
    }
    if ((node_array = first->node = xmalloc(256*sizeof(GifTree *)*NUMBER_OF_ARRAYS)) == NULL) {
        return 0;
    }
    last_array = node_array + (256*NUMBER_OF_ARRAYS - cc);
    scrngif_clear_tree(cc, first);

    pos = scrngif_add_code_to_buffer(cc, cLength, pos);

    end = pixels + siz;
    curNode = first;
    while (pixels < end) {
        if (curNode->node[*pixels] != NULL) {
            curNode = curNode->node[*pixels];
            tel++;
            pixels++;
            chainlen++;
            continue;
        } else if (curNode->typ == SEARCH) {
            newNode = curNode->nxt;
            while (newNode->alt != NULL) {
                if (newNode->ix == *pixels) {
                    break;
                }
                newNode = newNode->alt;
            }
            if (newNode->ix == *pixels) {
                tel++;
                pixels++;
                chainlen++;
                curNode = newNode;
                continue;
            }
        }

        /* ******************************************************
         * If there is no more thread to follow, we create a new node.  If the
         * current node is terminating, it will become a SEARCH node.  If it is
         * a SEARCH node, and if we still have room, it will be converted to a
         * LOOKUP node.
         */
        newNode = ++topNode;
        switch (curNode->typ) {
        case LOOKUP:
            newNode->nxt = NULL;
            newNode->alt = NULL;
            curNode->node[*pixels] = newNode;
            break;
        case SEARCH:
            if (node_array != last_array) {
                node_array += cc;
                curNode->node = node_array;
                curNode->typ = LOOKUP;
                curNode->node[*pixels] = newNode;
                curNode->node[(curNode->nxt)->ix] = curNode->nxt;
                lookuptypes++;
                newNode->nxt = NULL;
                newNode->alt = NULL;
                curNode->nxt = NULL;
                break;
            }
            //   otherwise do as we do with a TERMIN node
            __attribute__((fallthrough));
        case TERMIN:
            newNode->alt = curNode->nxt;
            newNode->nxt = NULL;
            curNode->nxt = newNode;
            curNode->typ = SEARCH;
            break;
        default:
            fprintf(stderr, "Silly node type: %d\n", curNode->typ);
        }
        newNode->code = next;
        newNode->ix = *pixels;
        newNode->typ = TERMIN;
        newNode->node = empty;
        nodecount++;
        /*
         * End of node creation
         * ******************************************************
         */
        if (debugFlag) {
            if (curNode == newNode) {
                fprintf(stderr, "Wrong choice of node\n");
            }
            if (curNode->typ == LOOKUP && curNode->node[*pixels] != newNode) {
                fprintf(stderr, "Wrong pixel coding\n");
            }
            if (curNode->typ == TERMIN) {
                fprintf(stderr, "Wrong Type coding; pixel# = %d; nodecount = %d\n", tel, nodecount);
            }
        }
        pos = scrngif_add_code_to_buffer(curNode->code, cLength, pos);
        if (chainlen > maxchainlen) {
            maxchainlen = chainlen;
        }
        chainlen = 0;
        if (pos - buffer > BLOKLEN) {
            buffer[-1] = BLOKLEN;
            fwrite(buffer - 1, 1, BLOKLEN + 1, fout);
            buffer[0] = buffer[BLOKLEN];
            buffer[1] = buffer[BLOKLEN + 1];
            buffer[2] = buffer[BLOKLEN + 2];
            buffer[3] = buffer[BLOKLEN + 3];
            pos -= BLOKLEN;
        }
        curNode = first;

        if (next == (1 << cLength)) {
            cLength++;
        }
        next++;

        if (next == 0xfff) {
            scrngif_clear_tree(cc, first);
            pos = scrngif_add_code_to_buffer(cc, cLength, pos);
            if (pos - buffer > BLOKLEN) {
                buffer[-1] = BLOKLEN;
                fwrite(buffer - 1, 1, BLOKLEN + 1, fout);
                buffer[0] = buffer[BLOKLEN];
                buffer[1] = buffer[BLOKLEN + 1];
                buffer[2] = buffer[BLOKLEN + 2];
                buffer[3] = buffer[BLOKLEN + 3];
                pos -= BLOKLEN;
            }
            next = cc + 2;
            cLength = (int16)((depth == 1) ? 3 : depth + 1);
        }
    }

    pos = scrngif_add_code_to_buffer(curNode->code, cLength, pos);
    if (pos - buffer > BLOKLEN - 3) {
        buffer[-1] = BLOKLEN - 3;
        fwrite(buffer - 1, 1, BLOKLEN - 2, fout);
        buffer[0] = buffer[BLOKLEN - 3];
        buffer[1] = buffer[BLOKLEN - 2];
        buffer[2] = buffer[BLOKLEN - 1];
        buffer[3] = buffer[BLOKLEN];
        buffer[4] = buffer[BLOKLEN + 1];
        pos -= BLOKLEN - 3;
    }
    pos = scrngif_add_code_to_buffer(eoi, cLength, pos);
    pos = scrngif_add_code_to_buffer(0x0, -1, pos);
    buffer[-1] = (uchar)(pos - buffer);
    pos = scrngif_add_code_to_buffer(0x0, 8, pos);

    fwrite(buffer - 1, (usize)(pos - buffer + 1), 1, fout);
    free(buffer - 1);
    free(first->node);
    free(baseNode);
    if (debugFlag) {
        fprintf(stderr, "pixel count = %d; nodeCount = %d lookup nodes = %d\n", tel, nodecount,
                lookuptypes);
    }
    return 1;
}

void
scrngif_clear_tree(int32 cc, GifTree *root) {
    GifTree *newNode, **xx;

    if (debugFlag > 1) {
        fprintf(stderr, "Clear Tree  cc= %d\n", cc);
    }
    if (debugFlag > 1) {
        fprintf(stderr, "nodeCount = %d lookup nodes = %d\n", nodecount, lookuptypes);
    }
    maxchainlen = 0;
    lookuptypes = 1;
    nodecount = 0;
    node_array = root->node;
    xx = node_array;
    for (int32 i = 0; i < NUMBER_OF_ARRAYS; i++) {
        memmove(xx, empty, 256*sizeof(GifTree **));
        xx += 256;
    }
    topNode = baseNode;
    for (int32 i = 0; i < cc; i++) {
        root->node[i] = newNode = ++topNode;
        newNode->nxt = NULL;
        newNode->alt = NULL;
        newNode->code = i;
        newNode->ix = (uchar)i;
        newNode->typ = TERMIN;
        newNode->node = empty;
        nodecount++;
    }
    return;
}

uchar *
scrngif_add_code_to_buffer(int32 code, int16 n, uchar *buf) {
    int32 mask;

    if (n < 0) {
        if (need < 8) {
            buf++;
            *buf = 0x0;
        }
        need = 8;
        return buf;
    }

    while (n >= need) {
        mask = (1 << need) - 1;
        *buf += (mask & code) << (8 - need);
        buf++;
        *buf = 0x0;
        code = code >> need;
        n -= need;
        need = 8;
    }
    if (n) {
        mask = (1 << n) - 1;
        *buf += (mask & code) << (8 - need);
        need -= n;
    }
    return buf;
}
