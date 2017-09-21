#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP cpbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP clbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cmbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP cpwbart(SEXP, SEXP, SEXP);
extern SEXP cwbart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mc_cores_openmp();

static const R_CallMethodDef CallEntries[] = {
    {"cpbart",  (DL_FUNC) &cpbart,  20},
    {"clbart",  (DL_FUNC) &clbart,  20},
    {"cmbart",  (DL_FUNC) &cmbart,  10},
    {"cpwbart", (DL_FUNC) &cpwbart,  3},
    {"cwbart",  (DL_FUNC) &cwbart,  23},
    {"mc_cores_openmp",(DL_FUNC) &mc_cores_openmp,0},
    {NULL, NULL, 0}
};

void R_init_BART(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
