// Generated with tools::package_native_routine_registration_skeleton()

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void auc_uno(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP auc_SZ(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Cham_Diao(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"auc_uno", (DL_FUNC)&auc_uno, 15},
    {NULL, NULL, 0}};

static const R_CallMethodDef CallEntries[] = {
    {"auc_SZ", (DL_FUNC)&auc_SZ, 13},
    {"Cham_Diao", (DL_FUNC)&Cham_Diao, 11},
    {NULL, NULL, 0}};

void R_init_hdnom(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
