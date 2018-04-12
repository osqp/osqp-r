#ifndef OSQP_CONFIGURE_H
# define OSQP_CONFIGURE_H

# ifdef __cplusplus
extern "C" {
# endif /* ifdef __cplusplus */

// Operating system (no cmake)
#ifdef _WIN32
#define IS_WINDOWS
#elif __APPLE__
#define IS_MAC
#else
#define IS_LINUX
#endif

// EMBEDDED
/* #undef EMBEDDED */

// PRINTING
#define PRINTING

// PROFILING
#define PROFILING

// CTRLC
#define CTRLC

// DFLOAT
/* #undef DFLOAT */

// DLONG
#define DLONG

// ENABLE_MKL_PARDISO
#define ENABLE_MKL_PARDISO


# ifdef __cplusplus
}
# endif /* ifdef __cplusplus */

#endif /* ifndef OSQP_CONFIGURE_H */
