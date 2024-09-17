#define SETVAR(i, x)                                                           \
    do {                                                                       \
        if ((i) < NVAR)                                                        \
            variables[(i)] = (x);                                              \
    } while (0)
#define GETVAR(i) (i) < NVAR ? variables[(i)] : 0.0
