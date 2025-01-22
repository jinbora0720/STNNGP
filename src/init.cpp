#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "STNNGP.h"

static const R_CallMethodDef CallEntries[] = {
    {"rSTNNGP", (DL_FUNC) &rSTNNGP, 35},
    {"rSTNNGP_NS", (DL_FUNC) &rSTNNGP_NS, 39},
    {"sSTNNGP", (DL_FUNC) &sSTNNGP, 41},
    {"sSTNNGP_NS", (DL_FUNC) &sSTNNGP_NS, 46},
    {"sSTNNGP_misalign", (DL_FUNC) &sSTNNGP_misalign, 44},
    {"sSTNNGP_NS_misalign", (DL_FUNC) &sSTNNGP_NS_misalign, 49},
    {"sSTNNGPPredict", (DL_FUNC) &sSTNNGPPredict, 23},
    {"sSTNNGPPredict_NS", (DL_FUNC) &sSTNNGPPredict_NS, 26},
    {"sSTNNGPLogit", (DL_FUNC) &sSTNNGPLogit, 44},
    {"sSTNNGPLogit_NS", (DL_FUNC) &sSTNNGPLogit_NS, 46},
    {NULL, NULL, 0}
};

void R_init_STNNGP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
