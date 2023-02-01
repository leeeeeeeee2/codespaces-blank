# Software dependencies

## Compiler

As MPR relies on (relatively) recent Fortran features, we require one of those compilers
- `gfortran` versions >v6.0 (routinely checked versions: 7.3 and 8.3)
- `nag` compiler version >v6.2 (routinely checked versions: 6.2)
- `intel` compiler version >v18.0 (routinely checked versions: 18.0, 19.0)

Other compilers have not been tested and likely do not work out of the box.

## Libraries

### FORCES

We rely on a set of independently maintained Fortran utility routines ([FORtran library for Computational Environmental Systems](https://git.ufz.de/chs/forces)). 
It in turn requires the NetCDF4 Fortran library for input/output parsing. See its documentation for installation details.

### flogging

This logging system aims to be simple to use, similar to existing write statements and therefore easy to implement in existing codes, while providing some convenient extra features not present in fortran by default.
MPR uses a [fork](https://github.com/gitporst/flogging) from the [original source](https://github.com/Exteris/flogging) with slight modifications of so it nicely integrates with all compilers and CMake.

## Helper tools

In addition to the compiler and some libraries, we rely on the following software:

- [cmake](https://cmake.org/download/) 
  (For compiling MPR; versions >3.5)
- [Python](https://www.python.org/downloads/)
  (For supporting scripts; >3.5)
    - Python library [f90nml](https://f90nml.readthedocs.io/en/latest/) (for interfacing the Fortran nml format), 

For running the complete stack of CI/CD scripts, we rely on: 
- Fortran library [pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit) (for UnitTesting the Fortran code)
- Python library [python-graphviz](https://pygraphviz.github.io/) (for supporting scripts)
- Python library [pytest](https://f90nml.readthedocs.io/en/latest/) (for unit testing of supporting scripts)
- [doxygen](https://www.doxygen.nl/index.html) and [graphviz](https://graphviz.org/) (for creating the documentation website)

