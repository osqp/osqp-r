# Changelog

## Version 1.0.0

### Breaking changes

- Replaced R6 with S7 classes. Model methods are now accessed via `@`
  (e.g., `model@Solve()`) instead of `$`.
- Renamed settings to match OSQP 1.0: `warm_start` -\> `warm_starting`,
  `polish` -\> `polishing`. The old names are accepted with deprecation
  warnings.

### New features

- Upgraded OSQP C library from 0.6.3 to 1.0.0.
- Added `WarmStart()` and `ColdStart()` methods for primal/dual warm
  starting.
- New OSQP 1.0 settings: `rho_is_vec`, `check_dualgap`, `cg_max_iter`,
  `cg_tol_reduction`, `cg_tol_fraction`, `cg_precond`, `profiler_level`.
- New solver info fields: `dual_obj_val`, `duality_gap`.

### Other changes

- Added vignette with examples for lasso regression, portfolio
  optimization, and warm starting.
- User-facing messages now use the `cli` package.
- Vendored QDLDL v0.1.8 for CRAN compatibility.
- OSQP C library printing now uses `Rprintf` instead of `printf`.
- Declared GNU make in `SystemRequirements`.
- `R CMD check` passes with `Status: OK`.

## Version 0.6.3.2

CRAN release: 2023-10-20

- CMAKE file fixes per R-exts

## Version 0.6.3.1

CRAN release: 2023-10-04

- Bug fix: typo in sparse coercion routines (`dgCMatrix` where
  `dtCMatrix` should be).

## Version 0.6.3

CRAN release: 2023-10-03

- Sync up to version [0.6.3 of OSQP
  release](https://github.com/osqp/osqp/releases/tag/v0.6.3)
- Fix `params.R` ([issue](https://github.com/osqp/osqp-r/issues/18)
  [\#18](https://github.com/osqp/osqp-r/issues/18))
- Added `time-limit` settings parameter per [PG’s
  code](https://github.com/osqp/osqp-r/pull/24)
- Added check for lower bounds not exceeding upper bounds ([issue
  29](https://github.com/osqp/osqp-r/issues/29))

## Version 0.6.0.8

CRAN release: 2023-01-31

- Fix prototype of `main` in `qdldl_sources/examplaes/example.c`
- Fix use of `cmake`; was previously bombing due to `proj.h` header
  being found from system areas, rather than locally
- Fix `configure` script so that R-exts prescriptions are followed
- Update requirements and files to default of C++17 (`src/Makevars`,
  `src/Makevars.win`, `src/osqp/Makefile`)

## Version 0.6.0.7

CRAN release: 2022-11-09

- Fix use of bitwise `&` with boolean operands (file
  `osqp_interface.cpp`)
- Fix function declaration without prototype (files `paradiso_loader.h`,
  `paradiso_interface.c`, `paradiso_loader.c`)

## Version 0.6.0.6

CRAN release: 2022-09-29

- Change maintainer to Balasubramanian Narasimhan

## Version 0.6.0.5

CRAN release: 2021-12-07

- (Missing news)

## Version 0.6.0.4

CRAN release: 2021-12-06

- (Missing news)

## Version 0.6.0.3 (28 October 2019) OSQP 0.6.0.3

CRAN release: 2019-10-28

- Patches to OSQP-R interface to handle infinite-valued bounds

## Version 0.6.0.2 (9 September 2019) OSQP 0.6.0.2

CRAN release: 2019-09-11

- Patches to OSQP-R build for solaris and fedora

## Version 0.6.0 (2 September 2019) OSQP 0.6.0

CRAN release: 2019-09-04

- Updated OSQP to version 0.6.0

## Version 0.5.0 (10 December 2018) OSQP 0.5.0

CRAN release: 2018-12-11

- Updated OSQP to version 0.5.0

## Version 0.4.1 (25 September 2018) OSQP 0.4.1

CRAN release: 2018-09-27

- Updated OSQP to version 0.4.1

## Version 0.4.0 (25 July 2018) OSQP 0.4.0

CRAN release: 2018-07-23

- Updated OSQP version

## Version 0.3.1 (11 June 2018) OSQP 0.3.1

CRAN release: 2018-07-01

- First release on CRAN with new name
