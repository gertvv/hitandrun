#include "har.h"
#include <R_ext/Visibility.h>

static const R_CallMethodDef callMethods[] = {
  { "hitandrun_har", (DL_FUNC) &hitandrun_har, 5 },
  { "hitandrun_bbReject", (DL_FUNC) &hitandrun_bbReject, 5 },
  { "hitandrun_simplexSample", (DL_FUNC) &hitandrun_simplexSample, 3 },
  { "hitandrun_hypersphereSample", (DL_FUNC) &hitandrun_hypersphereSample, 2 },
  { "hitandrun_sab", (DL_FUNC) &hitandrun_sab, 6 },
  { NULL, NULL, 0 }
};

void attribute_visible R_init_hitandrun(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
