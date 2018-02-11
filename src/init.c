#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP rCholWishart(SEXP, SEXP, SEXP);
extern SEXP rInvCholWishart(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"rCholWishart",    (DL_FUNC) &rCholWishart,    3},
    {"rInvCholWishart", (DL_FUNC) &rInvCholWishart, 3},
    {NULL, NULL, 0}
};

void R_init_matrixdist(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
