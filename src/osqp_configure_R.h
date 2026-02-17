#ifndef OSQP_CONFIGURE_H
# define OSQP_CONFIGURE_H

/* Operating system */
#ifdef _WIN32
#define IS_WINDOWS
#elif __APPLE__
#define IS_MAC
#else
#define IS_LINUX
#endif

/* Algebra backend */
#define OSQP_ALGEBRA_BUILTIN

/* Printing */
#define OSQP_ENABLE_PRINTING

/* Profiling — MUST be defined at compile time for timing info
   (setup_time, solve_time, etc.). Without this, all timing fields return 0. */
#define OSQP_ENABLE_PROFILING

/* Interrupt handling — safe for R: saves/restores previous SIGINT handler */
#define OSQP_ENABLE_INTERRUPT

/* No embedded mode */
/* #undef OSQP_EMBEDDED_MODE */

/* No float — use double */
/* #undef OSQP_USE_FLOAT */

/* No long long — use int (R's IntegerVector is 32-bit int) */
/* #undef OSQP_USE_LONG */

/* No code generation */
/* #undef OSQP_CODEGEN */

/* No derivatives */
/* #undef OSQP_ENABLE_DERIVATIVES */

/* No external profiler annotations */
/* #undef OSQP_PROFILER_ANNOTATIONS */

/* No debug mode */
/* #undef OSQP_ENABLE_DEBUG */

/* No custom memory */
/* #undef OSQP_CUSTOM_MEMORY */

/* No custom printing */
/* #undef OSQP_CUSTOM_PRINTING */

/* No settings packing */
/* #undef OSQP_PACK_SETTINGS */

#endif /* ifndef OSQP_CONFIGURE_H */
