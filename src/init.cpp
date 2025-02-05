#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "STNNGP.h"

static const R_CallMethodDef CallEntries[] = {
    {"rSTNNGP", (DL_FUNC) &rSTNNGP, 35},
    {"rSTNNGP_NS", (DL_FUNC) &rSTNNGP_NS, 39},
    {"sSTNNGP", (DL_FUNC) &sSTNNGP, 44},
    {"sSTNNGP_NS", (DL_FUNC) &sSTNNGP_NS, 48},
    {"sSTNNGP_misalign", (DL_FUNC) &sSTNNGP_misalign, 46},
    {"sSTNNGP_NS_misalign", (DL_FUNC) &sSTNNGP_NS_misalign, 51},
    {"sSTNNGPPredict", (DL_FUNC) &sSTNNGPPredict, 24},
    {"sSTNNGPPredict_NS", (DL_FUNC) &sSTNNGPPredict_NS, 27},
    {"sSTNNGPLogit", (DL_FUNC) &sSTNNGPLogit, 45},
    {"sSTNNGPLogit_NS", (DL_FUNC) &sSTNNGPLogit_NS, 47},
    {NULL, NULL, 0}
};

void R_init_STNNGP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
