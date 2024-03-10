#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* .C calls */
extern void my_RtailsMSS(
        double *Rlocation, double *x, double *d, double *logd, double *F,
        double *logF, double *cF, double *logcF);

static const R_CMethodDef CEntries[] = {
    {"my_RtailsMSS",    (DL_FUNC) &my_RtailsMSS,    14},
    {NULL, 0}
};

void R_init_FMStable(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
