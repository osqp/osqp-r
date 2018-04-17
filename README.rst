R interface for OSQP
====================

.. image:: https://travis-ci.org/oxfordcontrol/osqp-r.svg?branch=master
    :target: https://travis-ci.org/oxfordcontrol/osqp-r
.. image:: https://ci.appveyor.com/api/projects/status/bx1navxa474nhlpd/branch/master?svg=true
    :target: https://ci.appveyor.com/project/goulart-paul/osqp-r/branch/master
.. image:: https://www.r-pkg.org/badges/version/rosqp
    :target: https://www.r-pkg.org/pkg/rosqp
.. image:: https://cranlogs.r-pkg.org/badges/rosqp
    :target: https://www.r-pkg.org/pkg/rosqp

Provides R-bindings to `OSQP <http://osqp.readthedocs.io/>`__: the Operator
Splitting QP Solver.

The OSQP (Operator Splitting Quadratic Program) solver is a numerical
optimization package for solving problems in the form

::

    minimize        0.5 x' P x + q' x

    subject to      l <= A x <= u

where ``x in R^n`` is the optimization variable. The objective function
is defined by a positive semidefinite matrix ``P in S^n_+`` and vector
``q in R^n``. The linear constraints are defined by matrix
``A in R^{m x n}`` and vectors ``l in R^m U {-inf}^m``,
``u in R^m U {+inf}^m``.


Documentation
-------------

The interface is documented `here <http://osqp.readthedocs.io/>`__.
