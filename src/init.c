#include <R.h>
#include <Rinternals.h>
#include <Rversion.h>
#include <R_ext/Rdynload.h>

// Declare your C functions
extern SEXP RtailsMSS(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"RtailsMSS", (DL_FUNC) &RtailsMSS, 2},
    {NULL, NULL, 0}
};

void R_init_yourpackage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

