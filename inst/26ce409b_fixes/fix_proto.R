## This script fixes C prototypes to pass CRAN checks
## Because osqp is slow in accepting my PR, I have to resort to this

## We will prepend a comment line to C source to flag that a fix has been done
FIXED_FLAG <- "// CRAN PROTOTYPE FIXES DONE"
cran_prototypes_fixed <- function(lines, signal = FIXED_FLAG) {
  grepl(FIXED_FLAG, lines[1])
}

# Fix pardiso_interface.c
fname <- "osqp_sources/lin_sys/direct/pardiso/pardiso_interface.c"
lines <- readLines(fname)
## Skip if already fixed
if (!cran_prototypes_fixed(lines)) {
  ## Fix line 35
  lines[35] <- "c_int mkl_get_max_threads(void);"
  lines <- c("// CRAN PROTOTYPES FIXED", lines)
  writeLines(text = lines, con = fname)
}

# Fix pardiso_loader.h
fname <- "osqp_sources/lin_sys/direct/pardiso/pardiso_loader.h"
lines <- readLines(fname)
## Skip if already fixed
if (!cran_prototypes_fixed(lines)) {
  ## Fix line 22
  lines[22] <- "c_int lh_unload_pardiso(void);"
  lines <- c("// CRAN PROTOTYPES FIXED", lines)
  writeLines(text = lines, con = fname)
}

# Fix pardiso_loader.c
fname <- "osqp_sources/lin_sys/direct/pardiso/pardiso_loader.c"
lines <- readLines(fname)
## Skip if already fixed
if (!cran_prototypes_fixed(lines)) {
  ## Fix lines 25, 58, 90
  lines[25] <- "typedef int (*mkl_get_mt_t)(void);"
  lines[58] <- "c_int mkl_get_max_threads(void) {"
  lines[90] <- "c_int lh_unload_pardiso(void) {"
  lines <- c("// CRAN PROTOTYPES FIXED", lines)
  writeLines(text = lines, con = fname)
}
