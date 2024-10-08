#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "integers.h"
#include "functions.h"
#include "xmalloc.h"

#define NCMD 45  // add new commands as needed

#define MAKEC 0
#define XORFX 1
#define SILENT 2
#define CONVERT 3
#define NOICON 4
#define NEWSEED 5
#define ALLWIN 6
#define SETFILE 7
#define MSSTYLE 8
#define PWHITE 9
#define RUNNOW 10
#define BIGF 11
#define SMALLF 12
#define PARFILE 13
#define OUTFILE 14
#define ICFILE 15
#define FCOLOR 16
#define BCOLOR 17
#define BBITMAP 18
#define GRADS 19
#define MINWIDTH 20
#define MINHEIGHT 21
#define MWCOLOR 22
#define DWCOLOR 23
#define BELL 24
#define ITRNSETS 25
#define USET 26
#define RSET 27
#define INCLUDE 28
#define QSETS 29
#define QPARS 30
#define QICS 31
#define QUIET 32
#define LOGFILE 33
#define ANIFILE 34
#define VERSION 35
#define MKPLOT 36
#define PLOTFMT 37
#define NOOUT 38
#define DFDRAW 39
#define NCDRAW 40
#define READSET 42
#define WITH 43
#define EQUIL 44

static int32 use_intern_sets = 1;
static char setfilename[XPP_MAX_NAME];
static char parfilename[XPP_MAX_NAME];
static char icfilename[XPP_MAX_NAME];
char includefilename[MAX_INCLUDE_FILES][XPP_MAX_NAME];

static char readsetfile[XPP_MAX_NAME];
static int32 externaloptionsflag = 0;
static char externaloptionsstring[1024];
int32 NincludedFiles = 0;

static int32 select_intern_sets = 0;

int32 Nintern_2_use = 0;

static struct SetName {
    char *name;
    struct SetName *next;
} *sets2use, *setsNOTuse;

static int32 cli_is_set_name(struct SetName *set, char *nam);
static struct SetName *cli_add_set(struct SetName *set, char *nam);
static int32 cli_parse_it(char *com);

static int32 loadsetfile = 0;
static int32 loadparfile = 0;
static int32 loadicfile = 0;
int32 loadincludefile = 0;
int32 querysets = 0;
int32 querypars = 0;
int32 queryics = 0;
int32 dryrun = 0;
int32 noicon = 1;
int32 newseed = 0;

typedef struct Vocab {
    char name[10];
    int32 len;
} Vocab;

static Vocab my_cmd[NCMD] = {
    {"-m", 3},          {"-xorfix", 7},     {"-silent", 7},     {"-convert", 8},    {"-iconify", 7},
    {"-newseed", 7},    {"-allwin", 6},     {"-setfile", 7},    {"-ee", 3},         {"-white", 6},
    {"-runnow", 7},     {"-bigfont", 8},    {"-smallfont", 10}, {"-parfile", 8},    {"-outfile", 8},
    {"-icfile", 7},     {"-forecolor", 10}, {"-backcolor", 10}, {"-backimage", 10}, {"-grads", 6},
    {"-width", 6},      {"-height", 7},     {"-mwcolor", 8},    {"-dwcolor", 8},    {"-bell", 4},
    {"-internset", 10}, {"-uset", 5},       {"-rset", 5},       {"-include", 8},    {"-qsets", 6},
    {"-qpars", 6},      {"-qics", 5},       {"-quiet", 6},      {"-logfile", 8},    {"-anifile", 8},
    {"-version", 8},    {"-mkplot", 7},     {"-plotfmt", 8},    {"-noout", 6},      {"-dfdraw", 7},
    {"-ncdraw", 7},     {"-def", 4},        {"-readset", 8},    {"-with", 5},       {"-equil", 6}};

int32
cli_is_set_name(struct SetName *set, char *nam) {
    struct SetName *curr;
    if (set == NULL) {
        return 0;
    }

    curr = set;

    while (curr) {
        if (strcmp(curr->name, nam) == 0) {
            return 1;
        }
        curr = (struct SetName *)curr->next;
    }

    return 0;
}

struct SetName *
cli_add_set(struct SetName *set, char *nam) {
    if (!cli_is_set_name(set, nam)) {
        struct SetName *curr;
        curr = xmalloc(sizeof(struct SetName));
        curr->name = (char *)nam;
        curr->next = (struct SetName *)set;
        set = curr;
    }

    return set;
}

void
cli_do(int32 argc, char **argv) {
    int32 i;
    int32 k;

    silent = 0;
    got_file = 0;
    xorfix = 1;
    setfilename[0] = 0;
    parfilename[0] = 0;
    icfilename[0] = 0;
    for (i = 1; i < argc; i++) {
        k = cli_parse_it(argv[i]);
        if (k == 1) {
            strcpy(setfilename, argv[i + 1]);
            i++;
            loadsetfile = 1;
        }
        if (k == 2) {
            if (not_already_set.SMALL_FONT_NAME) {
                strcpy(font_name_small, argv[i + 1]);
                not_already_set.SMALL_FONT_NAME = 0;
            }
            i++;
        }
        if (k == 3) {
            if (not_already_set.BIG_FONT_NAME) {
                strcpy(font_name_big, argv[i + 1]);
                not_already_set.BIG_FONT_NAME = 0;
            }
            i++;
        }
        if (k == 4) {
            strcat(parfilename, "!load ");
            strcat(parfilename, argv[i + 1]);
            i++;
            loadparfile = 1;
        }
        if (k == 5) {
            ggets_plintf(argv[i + 1]);
            strncpy(batch_out, argv[i + 1], sizeof(batch_out));
            strncpy(user_out_file, argv[i + 1], sizeof(user_out_file));
            i++;
        }
        if (k == 6) {
            strcat(icfilename, argv[i + 1]);
            i++;
            loadicfile = 1;
        }
        if (k == 7) {
            if (strlen(argv[i + 1]) != 6) {
                ggets_plintf("Color must be given as hexadecimal string.\n");
                exit(-1);
            }
            load_eqn_set_option("FORECOLOR", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 8) {
            if (strlen(argv[i + 1]) != 6) {
                ggets_plintf("Color must be given as hexadecimal string.\n");
                exit(-1);
            }
            load_eqn_set_option("BACKCOLOR", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 9) {
            load_eqn_set_option("BACKIMAGE", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 10) {
            load_eqn_set_option("GRADS", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 11) {
            load_eqn_set_option("WIDTH", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 12) {
            load_eqn_set_option("HEIGHT", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 13) {
            if (strlen(argv[i + 1]) != 6) {
                ggets_plintf("Color must be given as hexadecimal string.\n");
                exit(-1);
            }
            load_eqn_set_option("MWCOLOR", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 14) {
            if (strlen(argv[i + 1]) != 6) {
                ggets_plintf("Color must be given as hexadecimal string.\n");
                exit(-1);
            }
            load_eqn_set_option("DWCOLOR", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 15) {
            load_eqn_set_option("BELL", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 16) {
            use_intern_sets = atoi(argv[i + 1]);
            select_intern_sets = 1;
            i++;
        }
        if (k == 17) {
            sets2use = cli_add_set(sets2use, argv[i + 1]);
            i++;
            select_intern_sets = 1;
        }
        if (k == 18) {
            setsNOTuse = cli_add_set(setsNOTuse, argv[i + 1]);
            i++;
            select_intern_sets = 1;
        }
        if (k == 19) {
            if (NincludedFiles > MAX_INCLUDE_FILES) {
                printf("Max number of include files exceeded.\n");
            }
            strcpy(includefilename[NincludedFiles], argv[i + 1]);
            NincludedFiles++;
            i++;
            loadincludefile = 1;
        }
        if (k == 20) {
            load_eqn_set_option("QUIET", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 21) {
            load_eqn_set_option("LOGFILE", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 22) {
            strcpy(anifile, argv[i + 1]);
            use_ani_file = 1;
            i++;
        }
        if (k == 23) {
            printf("XPPAUT Version %g.%g\nCopyright 2015 Bard Ermentrout\n", (double)MAJOR_VERSION,
                   (double)MINOR_VERSION);
            exit(0);
        }
        if (k == 24) {
            load_eqn_set_option("PLOTFMT", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 25) {
            suppress_out = 1;
        }
        if (k == 26) {
            load_eqn_set_option("DFDRAW", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 27) {
            load_eqn_set_option("NCDRAW", argv[i + 1], 1, NULL);
            i++;
        }
        if (k == 28) {  // -readset
            strcpy(readsetfile, argv[i + 1]);
            i++;
            externaloptionsflag = 1;
        }
        if (k == 29) {  // -with
            strcpy(externaloptionsstring, argv[i + 1]);
            i++;
            externaloptionsflag = 2;
        }
        if (k == 30) {  // -equil
            batch_equil = atoi(argv[i + 1]);
            i++;
            printf(" Batch equilibria %d \n", batch_equil);
        }
    }
    return;
}

int32
cli_if_needed_load_ext_options(void) {
    FILE *fp;
    char myopts[1024];
    char myoptsx[1026];
    if (externaloptionsflag == 0) {
        return 1;
    }

    if (externaloptionsflag == 1) {
        fp = fopen(readsetfile, "r");
        if (fp == NULL) {
            ggets_plintf("%s external set not found\n", readsetfile);
            return 0;
        }
        fgets(myopts, 1024, fp);
        sprintf(myoptsx, "$ %s", myopts);
        ggets_plintf("Got this string: {%s}\n", myopts);
        load_eqn_extract_action(myoptsx);
        fclose(fp);
        return 1;
    }

    if (externaloptionsflag == 2) {
        sprintf(myoptsx, "$ %s", externaloptionsstring);
        load_eqn_extract_action(myoptsx);
        return 1;
    }
    return -1;
}

int32
cli_if_needed_select_sets(void) {
    if (!select_intern_sets) {
        return 1;
    }
    for (int32 j = 0; j < Nintern_set; j++) {
        intern_set[j].use = (uint32)use_intern_sets;
        Nintern_2_use += use_intern_sets;

        if (cli_is_set_name(sets2use, intern_set[j].name)) {
            ggets_plintf("Internal set %s was included\n", intern_set[j].name);
            if (intern_set[j].use == 0) {
                Nintern_2_use++;
            }
            intern_set[j].use = 1;
        }

        if (cli_is_set_name(setsNOTuse, intern_set[j].name)) {
            ggets_plintf("Internal set %s was excluded\n", intern_set[j].name);
            if (intern_set[j].use == 1) {
                Nintern_2_use--;
            }
            intern_set[j].use = 0;
        }
    }

    ggets_plintf("A total of %d internal sets will be used\n", Nintern_2_use);

    return 1;
}

int32
cli_if_needed_load_set(void) {
    FILE *fp;
    if (!loadsetfile) {
        return 1;
    }
    fp = fopen(setfilename, "r");
    if (fp == NULL) {
        ggets_plintf("Couldn't load %s\n", setfilename);
        return 0;
    }
    lunch_read(fp);
    fclose(fp);
    return 1;
}

int32
cli_if_needed_load_par(void) {
    if (!loadparfile) {
        return 1;
    }
    ggets_plintf("Loading external parameter file: %s\n", parfilename);
    lunch_io_parameter_file(parfilename, 1);
    return 1;
}

int32
cli_if_needed_load_ic(void) {
    if (!loadicfile) {
        return 1;
    }
    ggets_plintf("Loading external initial condition file: %s\n", icfilename);
    lunch_io_ic_file(icfilename, 1);
    return 1;
}

int32
cli_parse_it(char *com) {
    int32 j;
    for (j = 0; j < NCMD; j++) {
        if (strncmp(com, my_cmd[j].name, (size_t)my_cmd[j].len) == 0) {
            break;
        }
    }

    if (j < NCMD) {
        switch (j) {
        case MAKEC:
            ggets_plintf(" C files are no longer part of this version. \n Sorry \n");
            break;
        case MKPLOT:
            make_plot_flag = 1;
            break;
        case SILENT:
            xpp_batch = 1;
            break;
        case XORFX:
            xorfix = 0;
            break;
        case CONVERT:
            convert_style = 1;
            break;
        case NOICON:
            noicon = 0;
            break;
        case NEWSEED:
            ggets_plintf("Random number seed changed\n");
            newseed = 1;
            break;
        case ALLWIN:
            all_win_vis = 1;
            break;
        case MSSTYLE:
            ms_style = 1;
            break;
        case PWHITE:
            ggets_plintf("-white option is no longer part of this version. \n Sorry \n");
            break;
        case RUNNOW:
            run_immediately = 1;
            break;
        case SETFILE:
            return 1;
        case SMALLF:
            return 2;
        case BIGF:
            return 3;
        case PARFILE:
            return 4;
        case OUTFILE:
            return 5;
        case ICFILE:
            return 6;
        case FCOLOR:
            return 7;
        case BCOLOR:
            return 8;
        case BBITMAP:
            return 9;
        case GRADS:
            return 10;
        case MINWIDTH:
            return 11;
        case MINHEIGHT:
            return 12;
        case MWCOLOR:
            return 13;
        case DWCOLOR:
            return 14;
        case BELL:
            return 15;
        case ITRNSETS:
            return 16;
        case USET:
            return 17;
        case RSET:
            return 18;
        case INCLUDE:
            return 19;
        case QUIET:
            return 20;
        case LOGFILE:
            return 21;
        case ANIFILE:
            return 22;
        case VERSION:
            return 23;
        case PLOTFMT:
            return 24;
        case NOOUT:
            return 25;
        case DFDRAW:
            return 26;
        case NCDRAW:
            return 27;
        case READSET:
            return 28;
        case WITH:
            return 29;
        case EQUIL:
            return 30;
        case QSETS:
            xpp_batch = 1;
            querysets = 1;
            dryrun = 1;
            break;
        case QPARS:
            xpp_batch = 1;
            querypars = 1;
            dryrun = 1;
            break;
        case QICS:
            xpp_batch = 1;
            queryics = 1;
            dryrun = 1;
            break;
        default:
            break;
        }
    } else {
        if (com[0] == '-' || got_file == 1) {
            ggets_plintf("Problem reading option %s\n", com);
            ggets_plintf("\nUsage: xppaut filename [options ...]\n\n");
            ggets_plintf("Options:\n");
            ggets_plintf("  -silent                Batch run without the interface and "
                         "dump solutions to a file\n");
            ggets_plintf("  -xorfix                Work-around for exclusive Or with "
                         "X on "
                         "some monitors/graphics setups\n");
            ggets_plintf("  -convert               Convert old style ODE files (e.g. "
                         "phaseplane) to new ODE style\n");
            ggets_plintf("  -newseed               Randomizes the random number "
                         "generator "
                         "which will often use the same seed\n");
            ggets_plintf("  -ee                    Emulates shortcuts of Evil Empire "
                         "style "
                         "(MS)\n");
            ggets_plintf("  -allwin                Brings XPP up with all the windows "
                         "visible\n");
            ggets_plintf("  -white                 Uses white screen instead of "
                         "black\n");
            ggets_plintf("  -setfile <filename>    Loads the set file before "
                         "starting up\n");
            ggets_plintf("  -runnow                Runs ode file immediately upon "
                         "startup "
                         "(implied by -silent)\n");
            ggets_plintf("  -bigfont <font>        Use the big font whose filename is "
                         "given\n");
            ggets_plintf("  -smallfont <font>      Use the small font whose filename is "
                         "given\n");
            ggets_plintf("  -parfile <filename>    Load parameters from the named "
                         "file\n");
            ggets_plintf("  -outfile <filename>    Send output to this file (default is "
                         "output.dat)\n");
            ggets_plintf("  -icfile <filename>     Load initial conditions from the "
                         "named "
                         "file\n");
            ggets_plintf("  -forecolor <######>    Hexadecimal color (e.g. 000000) for "
                         "foreground\n");
            ggets_plintf("  -backcolor <######>    Hexadecimal color (e.g. EDE9E3) for "
                         "background\n");
            ggets_plintf("  -backimage <filename>  Name of bitmap file (.xbm) to "
                         "load in "
                         "background\n");
            ggets_plintf("  -mwcolor <######>      Hexadecimal color (e.g. 808080) for "
                         "main window\n");
            ggets_plintf("  -dwcolor <######>      Hexadecimal color (e.g. FFFFFF) for "
                         "drawing window\n");
            ggets_plintf("  -grads < 1 | 0 >       Color gradients will | won't be "
                         "used\n");
            ggets_plintf("  -width N               Minimum width in pixels of main "
                         "window\n");
            ggets_plintf("  -height N              Minimum height in pixels of main "
                         "window\n");
            ggets_plintf("  -bell < 1 | 0 >        Events will | won't trigger "
                         "system bell\n");
            ggets_plintf("  -internset < 1 | 0 >   Internal sets will | won't be run "
                         "during batch run\n");
            ggets_plintf("  -uset <setname>        Named internal set will be run "
                         "during "
                         "batch run\n");
            ggets_plintf("  -rset <setname>        Named internal set will not be run "
                         "during batch run\n");
            ggets_plintf("  -include <filename>    Named file will be included (see "
                         "#include directive)\n");
            ggets_plintf("  -qsets                 Query internal sets (output saved to "
                         "OUTFILE)\n");
            ggets_plintf("  -qpars                 Query parameters (output saved to "
                         "OUTFILE)\n");
            ggets_plintf("  -qics                  Query initial conditions (output "
                         "saved "
                         "to OUTFILE)\n");
            ggets_plintf("  -quiet <1 |0>          Do not print *anything* out to "
                         "console\n");
            ggets_plintf("  -logfile <filename>    Print console output to specified "
                         "logfile \n");
            ggets_plintf("  -anifile <filename>    Load an animation code file "
                         "(.ani) \n");
            ggets_plintf("  -plotfmt <svg|ps>       Set Batch plot format\n");
            ggets_plintf("  -mkplot                Do a plot in batch mode \n");
            ggets_plintf(" -ncdraw 1|2               Draw nullclines in batch (1) to "
                         "file "
                         "(2) \n");
            ggets_plintf(" -dfdraw 1-5       Draw dfields in batch (1-3) to file "
                         "(4-5)  \n");
            ggets_plintf("  -version               Print XPPAUT version and exit \n");
            ggets_plintf("  -readset <filename>   Read in a set file like the internal "
                         "sets\n");
            ggets_plintf("  -with string   String must be surrounded with quotes; "
                         "anything "
                         "that is in an internal set is valid\n");
            ggets_plintf("  -equil <0|1>    Write equilibria to equil.dat and if <1> "
                         "manifolds um1.dat,...,sm2.dat\n");
            ggets_plintf("\n");

            ggets_plintf("Environment variables:\n");
            ggets_plintf("  XPPHELP                Path to XPPAUT documentation file "
                         "<xpphelp.html>\n");
            ggets_plintf("  XPPBROWSER             Web browser (e.g. "
                         "/usr/bin/firefox)\n");
            ggets_plintf("  XPPSTART               Path to start looking for ODE "
                         "files\n");
            ggets_plintf("\n");
            exit(0);
        } else {
            strcpy(this_file, com);
            got_file = 1;
        }
    }
    return 0;
}
