## This script fixes C prototypes to pass CRAN checks
## It also renames some header files and edits files that use those
## headers appropriately to avoid collisions with well-know library
## headers that may be instead found by CMAKE (e.g. proj.h)

#' Replace C/C++ source file lines with replacements
#' @param path relative file path of source
#' @param line_no line numbers to be replaced (integer vector)
#' @param replacement the corresponding replacements (character vector)
#' @param comment_prefix the comment prefix, default C/C++
replace_lines <- function(path, line_no, replacement, comment_prefix = "//") {
  stopifnot(length(line_no) == length(replacement))
  FIXED_FLAG <- paste(comment_prefix, "CRAN FIXES DONE")
  lines <- readLines(path)
  ## Skip if already fixed
  if (!grepl(pattern = FIXED_FLAG, x = lines[1])) {
    lines[line_no] <- replacement
    lines <- c(FIXED_FLAG, lines)
    writeLines(text = lines, con = path)
  }
  invisible(TRUE)
}


## # Fix pardiso_interface.c
## replace_lines(path = "osqp_sources/lin_sys/direct/pardiso/pardiso_interface.c",
##               line_no = 35,
##               replacement = "c_int mkl_get_max_threads(void);")

## # Fix pardiso_loader.h
## replace_lines(path = "osqp_sources/lin_sys/direct/pardiso/pardiso_loader.h",
##               line_no = 22,
##               replacement = "c_int lh_unload_pardiso(void);")

## # Fix pardiso_loader.c
## replace_lines(path = "osqp_sources/lin_sys/direct/pardiso/pardiso_loader.c",
##               line_no = c(25, 58, 90),
##               replacement = c("typedef int (*mkl_get_mt_t)(void);",
##                               "c_int mkl_get_max_threads(void) {",
##                               "c_int lh_unload_pardiso(void) {"))
# Fix example.c
replace_lines(path = "osqp_sources/lin_sys/direct/qdldl/qdldl_sources/examples/example.c",
              line_no = 17,
              replacement = "int main(int argc, char* argv[])")

## Rename proj.h to osqp_proj.h and edit files that use it.
file.rename(from = "osqp_sources/include/proj.h", to = "osqp_sources/include/osqp_proj.h")

## Fix auxil.c
replace_lines("osqp_sources/src/auxil.c", 3, '#include "osqp_proj.h"')

## Fix polish.c
replace_lines("osqp_sources/src/polish.c", 7, '#include "osqp_proj.h"')

## Fix proj.c
replace_lines("osqp_sources/src/proj.c", 1, '#include "osqp_proj.h"')

## We don't use CMAKE anymore, so not needed
## Fix CMakeLists.txt
## replace_lines("osqp_sources/include/CMakeLists.txt", 12, '    "${CMAKE_CURRENT_SOURCE_DIR}/osqp_proj.h"',
##               comment_prefix = "#")

## replace_lines("osqp_sources/lin_sys/direct/qdldl/qdldl_sources/CMakeLists.txt", 2, 'cmake_minimum_required (VERSION 3.5)',
##               comment_prefix = "#")

## replace_lines("osqp_sources/CMakeLists.txt", 2, 'cmake_minimum_required (VERSION 3.5)',
##               comment_prefix = "#")

## Drop lines 242-255  as they are not needed
## lines <- readLines("osqp_sources/CMakeLists.txt")[-(242:255)]
##writeLines(lines, "osqp_sources/CMakeLists.txt")

## Replace findR.cmake

## replace_lines("osqp_sources/configure/cmake/FindR.cmake", 11, 'find_program(R_EXEC NAMES R R.exe PATHS ${R_HOME}/bin)',
##               comment_prefix='#')


