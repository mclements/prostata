#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* INFO:
   https://stackoverflow.com/questions/42313373/r-cmd-check-note-found-no-calls-to-r-registerroutines-r-usedynamicsymbols
   https://github.com/quanteda/quanteda/issues/749
   https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#index-Registering-native-routines
*/

/* .C calls */
extern void r_get_user_random_seed(void *);
extern void r_next_rng_substream();
extern void r_set_user_random_seed(void *);

/* .Call calls */
extern SEXP r_create_current_stream();
extern SEXP r_remove_current_stream();

static const R_CMethodDef CEntries[] = {
    {"r_get_user_random_seed", (DL_FUNC) &r_get_user_random_seed, 1},
    {"r_next_rng_substream",   (DL_FUNC) &r_next_rng_substream,   0},
    {"r_set_user_random_seed", (DL_FUNC) &r_set_user_random_seed, 1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"r_create_current_stream",   (DL_FUNC) &r_create_current_stream,   0},
    {"r_remove_current_stream",   (DL_FUNC) &r_remove_current_stream,   0},
    {NULL, NULL, 0}
};

void R_init_microsimulation(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
