# Version 0.6.3.2

* CMAKE file fixes per R-exts

# Version 0.6.3.1

* Bug fix: typo in sparse coercion routines (`dgCMatrix` where
  `dtCMatrix` should be).

# Version 0.6.3

* Sync up to version [0.6.3 of OSQP release](https://github.com/osqp/osqp/releases/tag/v0.6.3)
* Fix `params.R` ([issue #18](https://github.com/osqp/osqp-r/issues/18))
* Added `time-limit` settings parameter per [PG's code](https://github.com/osqp/osqp-r/pull/24)
* Added check for lower bounds not exceeding upper bounds ([issue 29](https://github.com/osqp/osqp-r/issues/29))

# Version 0.6.0.8

* Fix prototype of `main` in `qdldl_sources/examplaes/example.c`
* Fix use of `cmake`; was previously bombing due to `proj.h` header
  being found from system areas, rather than locally
* Fix `configure` script so that R-exts prescriptions are followed
* Update requirements and files to default of C++17 (`src/Makevars`,
  `src/Makevars.win`, `src/osqp/Makefile`)

# Version 0.6.0.7

* Fix use of bitwise `&` with boolean operands (file `osqp_interface.cpp`)
* Fix function declaration without prototype (files
  `paradiso_loader.h`, `paradiso_interface.c`, `paradiso_loader.c`)

# Version 0.6.0.6

* Change maintainer to Balasubramanian Narasimhan

# Version 0.6.0.5

* (Missing news)

# Version 0.6.0.4

* (Missing news)

# Version 0.6.0.3 (28 October 2019) OSQP 0.6.0.3

* Patches to OSQP-R interface to handle infinite-valued bounds

# Version 0.6.0.2 (9 September 2019) OSQP 0.6.0.2

* Patches to OSQP-R build for solaris and fedora

# Version 0.6.0 (2 September 2019) OSQP 0.6.0

* Updated OSQP to version 0.6.0

# Version 0.5.0 (10 December 2018) OSQP 0.5.0

* Updated OSQP to version 0.5.0

# Version 0.4.1 (25 September 2018) OSQP 0.4.1

* Updated OSQP to version 0.4.1

# Version 0.4.0 (25 July 2018) OSQP 0.4.0

* Updated OSQP version

# Version 0.3.1 (11 June 2018) OSQP 0.3.1

* First release on CRAN with new name
