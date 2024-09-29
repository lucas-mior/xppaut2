#include <stdlib.h>
#include "integers.h"
#include <stdbool.h>

#include <sys/wait.h>
#include <unistd.h>
/* the menu for XPP commands
 * this calls any command
 * it also has lots of the direct X Gui stuff

*/

#include "functions.h"
#include "parserslow.h"
#include "alert.bitmap"

#include <stdio.h>
#include <string.h>
#include <X11/Xlib.h>

/*
When adding items try to keep to the format of using a strong _verb_ to
start the help topic (e.g. "Press" or "Add") and keep it brief.
*/
#define N_TUTORIAL 12

static char *tutorial[N_TUTORIAL] = {
    "use captital letters on buttons as keyboard shortcuts. Press <n> or <d> "
    "now!",
    "use the <Esc> key to close any unwanted menus.",
    "press the <Tab> key to match/limit directory listings in the file "
    "browser.",
    "zoom in/out using single clicks and/or dragging regions.",
    "set your preferences in your .xpprc file (shortcut='fx').",
    "add buttons for your favorite shortcuts in your .xpprc (e.g. @ "
    "BUT=Fit:wf)",
    "get continuous updating using the left mouse button on a Par/Var slider.",
    "use the mouse, arrows, <PgUp>, <PgDn> etc keys to quickly move around in "
    "AUTO.",
    "link to the HTML manual by setting the environment variable XPPHELP on "
    "your computer.",
    "view documentation using your favorite web browser by setting the "
    "environment variable XPPBROWSER on your computer.",
    "edit files using your favorite text editor by setting the environment "
    "variable XPPEDITOR on your computer.",
    "have XPP open to a default starting directory by setting the environment "
    "variable XPPSTART on your computer."};

static int32 status;

static struct MsgBoxStruct {
    Window window;
    char text[256];
    int32 here;
} MsgBox;

void
menudrive_do_tutorial(void) {
    int32 tut = 0;
    printf("Running tutorial!\n");
    while (true) {
        char ans = (char)pop_list_two_choice(
            "Next", "Done", tutorial[tut], "nd", display_width / 2,
            display_height / 2, RootWindow(display, screen),
            "Did you know you can...");

        if (ans == 'd') {
            break;
        }
        tut++;
        tut = tut % N_TUTORIAL;
    }
    return;
}

void
menudrive_edit_xpprc(void) {
    pid_t child_pid;

    char rc[256];
    char editor[256];
    int32 child_status;

    char *ed = getenv("XPPEDITOR");

    if ((ed == NULL) || (strlen(ed) == 0)) {
        ggets_err_msg("Environment variable XPPEDITOR needs to be set.");

        return;
    } else {
        strncpy(editor, ed, sizeof(editor));
    }

    child_pid = fork();

    if (child_pid == 0) {
        char *args[] = {editor, NULL, NULL};
        snprintf(rc, sizeof(rc), "%s/.xpprc", getenv("HOME"));
        args[1] = rc;

        execvp(editor, args);
        wait(&child_status);
        return;
    } else {
        if (child_pid == -1) {
            ggets_err_msg("Unable to fork process for editor.");
        }

        return;
    }
}

void
menudrive_xpp_hlp(void) {
    char cmd[256];

    if (getenv("XPPHELP") == NULL) {
        ggets_err_msg("Environment variable XPPHELP undefined.");
        return;
    }

    if (getenv("XPPBROWSER") == NULL) {
        ggets_err_msg("Environment variable XPPBROWSER undefined.");
        return;
    }

    snprintf(cmd, sizeof(cmd), "file:///%s", getenv("XPPHELP"));

    if (fork() == 0) {
        execlp(getenv("XPPBROWSER"), getenv("XPPHELP"), cmd, (char *)0);
        perror("Unable to open browser. Check your XPPBROWSER and XPPHELP "
               "environement variables.");
        exit(1);
    } else {
        wait(&status);
    }
    return;
}

void
menudrive_message_box(char *m) {
    int32 wid = (int32)strlen(m)*DCURX + 20;
    int32 hgt = 4*DCURY;
    MsgBox.window = pop_list_make_plain_window(RootWindow(display, screen),
                                               display_width / 2,
                                               display_height / 2, wid, hgt, 4);

    many_pops_make_icon((char *)alert_bits, alert_width, alert_height,
                        MsgBox.window);
    MsgBox.here = 1;
    pop_list_set_window_title(MsgBox.window, "Yo!");
    strcpy(MsgBox.text, m);
    ggets_ping();
    return;
}

void
menudrive_message_box_redraw(Window window) {
    if (window == MsgBox.window) {
        ggets_f_text(10, 2*DCURY, MsgBox.text, MsgBox.window);
    }
    return;
}

void
menudrive_message_box_kill(void) {
    if (MsgBox.here == 0) {
        return;
    }
    MsgBox.here = 0;
    browse_wait_a_sec(ClickTime);
    XDestroyWindow(display, MsgBox.window);
    return;
}

int32
menudrive_two_choice(char *c1, char *c2, char *q, char *key) {
    int32 choice = pop_list_two_choice(c1, c2, q, key, display_width / 2,
                                       display_height / 2,
                                       RootWindow(display, screen), NULL);
    return choice;
}

int32
menudrive_get_mouse_xy(int32 *x, int32 *y) {
    return ggets_mouse_xy(x, y, draw_win);
}

void
menudrive_flush_display(void) {
    XFlush(display);
    return;
}

void
menudrive_clear_draw_window(void) {
    main_clr_scrn();
    many_pops_hi_lite(draw_win);
    return;
}

void
menudrive_drw_all_scrns(void) {
    int32 me = manual_expose;
    int32 ic = current_pop;
    manual_expose = 0;
    if (SimulPlotFlag == 0) {
        main_redraw_all();
        manual_expose = me;
        return;
    }

    for (int32 i = 0; i < num_pops; i++) {
        many_pops_make_active(ActiveWinList[i], 1);
        main_redraw_all();
    }

    many_pops_make_active(ic, 1);
    many_pops_hi_lite(draw_win);
    manual_expose = me;
    return;
}

void
menudrive_clr_all_scrns(void) {
    int32 ic = current_pop;

    if (SimulPlotFlag == 0) {
        main_clr_scrn();
        many_pops_hi_lite(draw_win);
        return;
    }

    for (int32 i = 0; i < num_pops; i++) {
        many_pops_make_active(ActiveWinList[i], 1);
        main_clr_scrn();
    }

    many_pops_make_active(ic, 1);
    many_pops_hi_lite(draw_win);
    return;
}

void
menudrive_run_the_commands(int32 com) {
    if (com < 0) {
        return;
    }
    if (com <= MAX_M_I) {
        integrate_do_init_data(com);
        return;
    }
    if (com == M_C) {
        integrate_cont_integ();
        return;
    }

    if (com >= M_SG && com <= M_SC) {
        integrate_find_equilib_com(com - M_SG);
        return;
    }
    if (com >= M_NFF && com <= M_NFA) {
        nullcline_froz_cline_stuff_com(com - M_NFF);
        return;
    }

    if (com >= M_NN && com <= M_NS) {
        nullcline_new_clines_com(com - M_NN);
        return;
    }

    if (com >= M_DD && com <= M_DS) {
        nullcline_direct_field_com(com - M_DD);
        if ((com - M_DD) == 1) {
            return;
        }
        nullcline_create_new_cline();
        graf_par_redraw_the_graph();
        return;
    }

    if (com >= M_WW && com <= M_WS) {
        graf_par_window_zoom_com(com - M_WW);
        return;
    }

    if (com >= M_AA && com <= M_AC) {
        do_torus_com(com - M_AA);
        return;
    }

    if (com >= M_KC && com <= M_KM) {
        kinescope_do_movie_com(com - M_KC);
        return;
    }

    if (com >= M_GA && com <= M_GC) {
        graf_par_add_a_curve_com(com - M_GA);
        return;
    }

    if (com >= M_GFF && com <= M_GFO) {
        graf_par_freeze_com(com - M_GFF);
        return;
    }

    if (com >= M_GCN && com <= M_GCU) {
        graf_par_change_cmap_com(com - M_GCN);
        nullcline_redraw_dfield();

        return;
    }

    if (com == M_GFKK || com == M_GFKN) {
        graf_par_key_frz_com(com - M_GFKN);
        return;
    }
    if (com == M_UKE || com == M_UKV) {
        tabular_new_lookup_com(com - M_UKE);
        return;
    }
    if (com == M_R) {
        menudrive_drw_all_scrns();
        return;
    }

    if (com == M_EE) {
        menudrive_clr_all_scrns();
        DF_FLAG = 0;
        return;
    }

    if (com == M_X) {
        graf_par_xi_vs_t();
        return;
    }

    if (com == M_3) {
        graf_par_get_3d_com();
        return;
    }

    if (com == M_P) {
        init_conds_new_parameter();
        return;
    }

    if (com >= M_MC && com <= M_MS) {
        many_pops_do_windows_com(com - M_MC);
        return;
    }
     // CLONE 
    if (com >= M_FP && com <= M_FL) {
         // menudrive do file com 
        switch (com) {
        case M_FT:
            adjoints_do_transpose();
            break;
        case M_FG:
            many_pops_get_intern_set();
            break;
        case M_FI:
            flag_tips = 1 - flag_tips;
            break;
        case M_FP:
            txt_make_view();
            break;
        case M_FW:
            do_lunch(0);
            break;
        case M_FS:
            lunch_file_inf();
            break;
        case M_FA:
#ifdef AUTO
            auto_nox_win();
#endif
            break;
        case M_FC:
            calc_q_calc();
            break;
        case M_FR:
            do_lunch(1);
            break;
        case M_FB:
            tfBell = 1 - tfBell;
            break;
        case M_FH:
             //menudrive_xpp_hlp();
            break;
        case M_FX:
            menudrive_edit_xpprc();
            break;
        case M_FU:
            menudrive_do_tutorial();
            break;
        case M_FQ:
            if (pop_list_yes_no_box()) {
                main_bye_bye();
            }
            break;
        case M_FER:
            edit_rhs();
            break;
        case M_FEF:
            edit_rhs_functions();
            break;
        case M_FES:
            edit_rhs_save_as();
            break;
        case M_FEL:
            extra_load_new_dll();
            break;
        case M_FL:
            init_conds_clone_ode();
            break;
        default:
            break;
        }
        return;
    }

    if (com >= M_TT && com <= M_TS) {
        many_pops_do_gr_objs_com(com - M_TT);
        return;
    }
    if (com >= M_TEM && com <= M_TED) {
        many_pops_edit_object_com(com - M_TEM);
        return;
    }
    if (com >= M_BR && com <= M_BH) {
        pp_shoot_find_bvp_com(com - M_BR);
        return;
    }

    if (com >= M_V2 && com <= M_VT) {
        graf_par_change_view_com(com - M_V2);
    }
    if (com >= M_UAN && com <= M_UAR) {
        adjoints_make_adj_com(com - M_UAN);
    }
    if (com >= M_UCN && com <= M_UCA) {
        numerics_set_col_par_com(com - M_UCN);
    }
    if (com >= M_UPN && com <= M_UPP) {
        numerics_get_pmap_pars_com(com - M_UPN);
    }
    if (com >= M_UHN && com <= M_UH2) {
        markov_do_stochast_com(com - M_UHN);
    }
    if (com >= M_UT && com <= M_UC) {
        numerics_quick_num(com - M_UT);
    }
    return;
}

void
menudrive_do_stochast(void) {
    static char *n[] = {"New seed",  "Compute",     "Data",     "Mean",
                        "Variance",  "Histogram",   "Old hist", "Fourier",
                        "Power",     "fIt data",    "Stat",     "Liapunov",
                        "stAutocor", "Xcorrel etc", "spEc.dns", "2D-hist"};
    static char key[] = "ncdmvhofpislaxe2";
    Window temp = main_win;
    char ch;
    int32 i;
    ch = (char)pop_up_list(&temp, "Stochastic", n, key, 16, 10, 0, 10,
                           2*DCURY + 8, stoch_hint, info_pop, info_message);
    for (i = 0; i < 16; i++) {
        if (ch == key[i]) {
            break;
        }
    }

    if (i >= 0 && i < 16) {
        menudrive_run_the_commands(M_UHN + i);
    }
    return;
}

void
menudrive_get_pmap_pars(void) {
    static char *map[] = {"(N)one", "(S)ection", "(M)ax/min", "(P)eriod"};
    static char mkey[] = "nsmp";
    char ch;
    Window temp = main_win;
    int32 i;

    ch = (char)pop_up_list(&temp, "Poincare map", map, mkey, 4, 13, POIMAP, 10,
                           6*DCURY + 8, map_hint, info_pop, info_message);

    for (i = 0; i < 4; i++) {
        if (ch == mkey[i]) {
            break;
        }
    }

    if (i >= 0 && i < 4) {
        menudrive_run_the_commands(M_UPN + i);
    }
    return;
}

void
menudrive_set_col_par(void) {
    char ch;
    Window tempw = main_win;
    static char *n[] = {"(N)o color", "(V)elocity", "(A)nother quantity"};
    static char key[] = "nva";
    int32 i;
    ch = (char)pop_up_list(&tempw, "Color code", n, key, 3, 11, 0, 10,
                           12*DCURY + 8, color_hint, info_pop, info_message);
    for (i = 0; i < 3; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < 3) {
        menudrive_run_the_commands(i + M_UCN);
    }
    return;
}

void
menudrive_make_adj(void) {
    Window temp = main_win;
    static char *n[] = {"(N)ew adj", "(M)ake H",     "(A)djoint", "(O)rbit",
                        "(H)fun",    "(P)arameters", "(R)ange"};
    static char key[] = "nmaohpr";
    char ch;
    int32 i;
    ch = (char)pop_up_list(&temp, "Adjoint", n, key, 7, 10, 0, 10,
                           11*DCURY + 8, adj_hint, info_pop, info_message);
    for (i = 0; i < 7; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < 7) {
        menudrive_run_the_commands(M_UAN + i);
    }
    return;
}

void
menudrive_do_gr_objs(void) {
    char ch;
    int32 i;
    static char *list[] = {"(T)ext", "(A)rrow",      "(P)ointer", "(M)arker",
                           "(E)dit", "(D)elete all", "marker(S)"};
    static char key[] = "tapmeds";
    static char title[] = "Text,etc";
    static char *elist[] = {"(M)ove", "(C)hange", "(D)elete"};
    static char ekey[] = "mcd";
    static char etitle[] = "Edit";
    Window temp = main_win;
    ch = (char)pop_up_list(&temp, title, list, key, 7, 10, 0, 10,
                           10*DCURY + 8, text_hint, info_pop, info_message);
    if (ch == PAUSE_NUMBER) {
        return;
    }
    if (ch == 'e') {
        ch = (char)pop_up_list(&temp, etitle, elist, ekey, 3, 9, 0, 10,
                               10*DCURY + 8,

                               edit_hint, info_pop, info_message);

        if (ch == PAUSE_NUMBER) {
            return;
        }
        for (i = 0; i < 3; i++) {
            if (ch == ekey[i]) {
                break;
            }
        }
        if (i >= 0 && i < 3) {
            menudrive_run_the_commands(M_TEM + i);
        }
        return;
    }
    for (i = 0; i < 7; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < 7) {
        menudrive_run_the_commands(M_TT + i);
    }
    return;
}

void
menudrive_new_lookup(void) {
    static char *n[] = {"(E)dit", "(V)iew"};
    static char key[] = "ev";
    char ch;
    Window temp = main_win;
    if (NTable == 0) {
        return;
    }
    ch = (char)pop_up_list(&temp, "Tables", n, key, 2, 12, 1, 10,
                           11*DCURY + 8, tab_hint, info_pop, info_message);
    if (ch == key[0]) {
        menudrive_run_the_commands(M_UKE);
    }
    if (ch == key[1]) {
        menudrive_run_the_commands(M_UKV);
    }
    return;
}

void
menudrive_find_bvp(void) {
    static char *n[] = {"(R)ange", "(N)o show", "(S)how", "(P)eriodic"};
    static char key[] = "rnsp";
    char ch;
    int32 i;
    Window temp = main_win;
    ch = (char)pop_up_list(&temp, "Bndry Value Prob", n, key, 4, 16, 1, 10,
                           6*DCURY + 8, bvp_hint, info_pop, info_message);
    if (ch == PAUSE_NUMBER) {
        return;
    }
    for (i = 0; i < 4; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < 4) {
        menudrive_run_the_commands(M_BR + i);
    }
    return;
}

void
menudrive_change_view(void) {
    Window temp = main_win;
    static char *n[] = {"2D", "3D", "Array", "Toon"};
    static char key[] = "23at";
    char ch;
    int32 i;

    ch = (char)pop_up_list(&temp, "Axes", n, key, 4, 5, 0, 10, 13*DCURY + 8,
                           view_hint, info_pop, info_message);
    for (i = 0; i < 4; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < 4) {
        menudrive_run_the_commands(M_V2 + i);
    }
    return;
}

void
menudrive_do_windows(void) {
    int32 i;
    char ch;
    static char *list[] = {"(C)reate", "(K)ill all", "(D)estroy",   "(B)ottom",
                           "(A)uto",   "(M)anual",   "(S)imPlot On"};
    static char *list2[] = {"(C)reate",     "(K)ill all", "(D)estroy",
                            "(B)ottom",     "(A)uto",     "(M)anual",
                            "(S)imPlot Off"};
    static char key[] = "ckdbams";
    static char title[] = "Make window";
    Window temp = main_win;
    if (SimulPlotFlag == 0) {
        ch = (char)pop_up_list(&temp, title, list, key, 7, 11, 0, 10,
                               14*DCURY + 8, half_hint, info_pop,
                               info_message);
    } else {
        ch = (char)pop_up_list(&temp, title, list2, key, 7, 11, 0, 10,
                               14*DCURY + 8, half_hint, info_pop,
                               info_message);
    }
    for (i = 0; i < 7; i++) {
        if (ch == key[i]) {
            break;
        }
    }

    if (i >= 0 && i < 7) {
        menudrive_run_the_commands(M_MC + i);
    }
    return;
}

void
menudrive_add_a_curve(void) {
    int32 com = -1;
    static char *na[] = {"(A)dd curve",  "(D)elete last", "(R)emove all",
                         "(E)dit curve", "(P)ostscript",  "S(V)G",
                         "(F)reeze",     "a(X)es opts",   "exp(O)rt data",
                         "(C)olormap"};
    static char *nc[] = {"(N)ormal",   "(P)eriodic", "(H)ot",      "(C)ool",
                         "(B)lue-red", "(G)ray",     "c(U)behelix"};
    static char *nf[] = {"(F)reeze", "(D)elete",   "(E)dit",    "(R)emove all",
                         "(K)ey",    "(B)if.Diag", "(C)lr. BD", "(O)n freeze"};
    static char *nf2[] = {"(F)reeze",     "(D)elete",    "(E)dit",
                          "(R)emove all", "(K)ey",       "(B)if.Diag",
                          "(C)lr. BD",    "(O)ff freeze"};
    static char *nk[] = {"(N)o key", "(K)ey"};
    static char keya[] = "adrepvfxoc";
    static char keyc[] = "nphcbgu";
    static char keyf[] = "fderkbco";
    static char keyk[] = "nk";
    Window temp = main_win;
    char ch;
    int32 i;
    int32 j;
    ch = (char)pop_up_list(&temp, "Curves", na, keya, 10, 15, 0, 10,
                           8*DCURY + 8, graf_hint, info_pop, info_message);
    for (i = 0; i < 10; i++) {
        if (ch == keya[i]) {
            break;
        }
    }
    if (i == 6) {
        if (AutoFreezeFlag == 0) {
            ch = (char)pop_up_list(&temp, "Freeze", nf, keyf, 8, 15, 0, 10,
                                   8*DCURY + 8, frz_hint, info_pop,
                                   info_message);
        } else {
            ch = (char)pop_up_list(&temp, "Freeze", nf2, keyf, 8, 15, 0, 10,
                                   8*DCURY + 8, frz_hint, info_pop,
                                   info_message);
        }
        for (j = 0; j < 8; j++) {
            if (ch == keyf[j]) {
                break;
            }
        }
        if (j == 4) {
            ch = (char)pop_up_list(&temp, "Key", nk, keyk, 2, 9, 0, 10,
                                   8*DCURY + 8, no_hint, info_pop,
                                   info_message);
            if (ch == keyk[0]) {
                com = M_GFKN;
            }
            if (ch == keyk[1]) {
                com = M_GFKK;
            }
        } else {
            if (j >= 0 && j < 8) {
                com = M_GFF + j;
            }
        }
    } else {
        if (i == 9) {
            ch = (char)pop_up_list(&temp, "Colormap", nc, keyc, 7, 15, 0, 10,
                                   8*DCURY + 8, cmap_hint, info_pop,
                                   info_message);
            for (j = 0; j < 7; j++) {
                if (ch == keyc[j]) {
                    break;
                }
            }
            if (j >= 0 && j < 7) {
                com = M_GCN + j;
            }
        } else {
            if (i >= 0 && i < 10) {
                com = M_GA + i;
            }
        }
    }
    menudrive_run_the_commands(com);
}

void
menudrive_do_movie(void) {
    int32 i;
    char ch;
    int32 nkc = 6;
    static char *list[] = {"(C)apture",  "(R)eset", "(P)layback",
                           "(A)utoplay", "(S)ave",  "(M)ake AniGif",
                           "(X)tra"};
    static char key[] = "crpasmx";
    Window temp = main_win;
    ch = (char)pop_up_list(&temp, "Kinescope", list, key, nkc, 11, 0, 10,
                           8*DCURY + 8, kin_hint, info_pop, info_message);
    for (i = 0; i < nkc; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < nkc) {
        menudrive_run_the_commands(i + M_KC);
    }
    return;
}

void
menudrive_do_torus(void) {
    Window temp = main_win;
    static char *n[] = {"(A)ll", "(N)one", "(C)hoose"};
    static char key[] = "anc";
    char ch;
    int32 i;
    ch = (char)pop_up_list(&temp, "Torus", n, key, 3, 9, 1 - TORUS, 10,
                           4*DCURY + 8, phas_hint, info_pop, info_message);
    for (i = 0; i < 3; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < 3) {
        menudrive_run_the_commands(M_AA + i);
    }
    return;
}

void
menudrive_window_zoom(void) {
    static char *n[] = {"(W)indow", "(Z)oom In", "Zoom (O)ut",
                        "(F)it",    "(D)efault", "(S)croll"};
    static char key[] = "wzofds";
    char ch;
    int32 i;
    Window temp = main_win;
    ch = (char)pop_up_list(&temp, "Window", n, key, 6, 13, 0, 10,
                           13*DCURY + 8, wind_hint, info_pop, info_message);
    for (i = 0; i < 6; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < 6) {
        menudrive_run_the_commands(M_WW + i);
    }
    return;
}

void
menudrive_direct_field(void) {
    int32 i;
    static char *n[] = {"(D)irect Field", "(F)low", "(N)o dir. fld.",
                        "(C)olorize", "(S)caled Dir.Fld"};
    static char key[] = "dfncs";
    char ch;
    Window temp = main_win;
    ch = (char)pop_up_list(&temp, "Two-D Fun", n, key, 5, 18, 0, 10,
                           6*DCURY + 8, flow_hint, info_pop, info_message);

    for (i = 0; i < 5; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < 5) {
        menudrive_run_the_commands(M_DD + i);
    }
    return;
}

void
menudrive_new_clines(void) {
    int32 i;
    Window temp = main_win;
    static char *n[] = {"(N)ew",    "(R)estore", "(A)uto",
                        "(M)anual", "(F)reeze",  "(S)ave"};
    static char key[] = "nramfs";
    char ch;
    ch = (char)pop_up_list(&temp, "Nullclines", n, key, 6, 10, 0, 10,
                           6*DCURY + 8, null_hint, info_pop, info_message);
    for (i = 0; i < 6; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < 6) {
        menudrive_run_the_commands(M_NN + i);
    }
    return;
}

void
menudrive_froz_cline_stuff(void) {
    Window temp = main_win;
    static char *n[] = {"(F)reeze", "(D)elete all", "(R)ange", "(A)nimate"};
    static char key[] = "fdra";
    char ch;
    int32 i;
    ch = (char)pop_up_list(&temp, "Freeze cline", n, key, 4, 10, 0, 10,
                           6*DCURY + 8, null_freeze, info_pop, info_message);
    for (i = 0; i < 4; i++) {
        if (ch == key[i]) {
            break;
        }
    }
    if (i >= 0 && i < 4) {
        menudrive_run_the_commands(M_NFF + i);
    }
    return;
}

void
menudrive_find_equilibrium(void) {
    int32 i;
    static char *n[] = {"(G)o", "(M)ouse", "(R)ange", "monte(C)ar"};
    static char key[] = "gmrc";
    char ch;
    Window temp = main_win;
    ch = (char)pop_up_list(&temp, "Equilibria", n, key, 4, 12, 1, 10,
                           6*DCURY + 8, sing_hint, info_pop, info_message);
    if (ch == PAUSE_NUMBER) {
        return;
    }
    for (i = 0; i < 4; i++) {
        if (ch == key[i]) {
            break;
        }
    }

    if (i > -1 && i < 4) {
        menudrive_run_the_commands(i + M_SG);
    }
    return;
}

void
menudrive_ini_data_menu(void) {
    int32 i;
    Window temp = main_win;
    static char *n[] = {"(R)ange",   "(2)par range", "(L)ast",    "(O)ld",
                        "(G)o",      "(M)ouse",      "(S)hift",   "(N)ew",
                        "s(H)oot",   "(F)ile",       "form(U)la", "m(I)ce",
                        "DAE guess", "(B)ackward"};
    static char key[] = "r2logmsnhfuidb";
    char ch;
    ch = (char)pop_up_list(&temp, "Integrate", n, key, 14, 13, 3, 10,
                           3*DCURY + 8, ic_hint, info_pop, info_message);

    if (ch == PAUSE_NUMBER) {
        return;
    }

    for (i = 0; i < 14; i++) {
        if (ch == key[i]) {
            break;
        }
    }

    menudrive_run_the_commands(i);
}

void
menudrive_new_param(void) {
    menudrive_run_the_commands(M_P);
}

void
menudrive_clear_screens(void) {
    menudrive_run_the_commands(M_EE);
}

void
menudrive_x_vs_t(void) {
    menudrive_run_the_commands(M_X);
}

void
menudrive_redraw_them_all(void) {
    menudrive_run_the_commands(M_R);
}

void
menudrive_get_3d_par(void) {
    menudrive_run_the_commands(M_3);
}
