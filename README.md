[![Build Status](https://travis-ci.org/felipeZ/Haskell-abinitio.svg?branch=master)](https://travis-ci.org/felipeZ/Haskell-abinitio)

Haskell-abinitio
================

 See http://felipez.github.io/Haskell-abinitio for documentation.

Description
-----------
This package can calculate the Hartree Fock energy
of a given molecule geometry and a basis set solving the
Roothaan Hall equations through a self consistent field
procedure. It uses the Harris Functional as an initial
density guess and the DIIS method to greatly improve
the convergence.
The entire code is written using the [Repa](https://hackage.haskell.org/package/repa)
library and focusing our efforts on efficiency, parallelism and code readability.
Using Haskell’s higher order abstraction we are
trying to develop an EDSL appropriate for quantum
mechanics problems, creating code operators able to
fairly mimic the physical ones.

The original idea of the project can be found at
[Haskell ab initio](https://themonadreader.files.wordpress.com/2013/03/issue214.pdf)
