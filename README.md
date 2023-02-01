# Multiscale parameter regionalization -- MPR

<div align="center">
<img src="https://git.ufz.de/chs/logos/-/raw/master/MPR.png" alt="MPR-LOGO" width="251px" style="width:251px;"/>
</div>

This is the MPR tool of the

> Department Computational Hydrosystems<br/>
> Helmholtz Centre for Environmental Research - UFZ<br/>
> Permoserstr. 15<br/>
> 04318 Leipzig, Germany

- The current release is **[MPR v1.0.4][1]**.
- The latest MPR release notes can be found in the file [RELEASES][3] or [online][4].
- General information can be found on the [MPR website](https://www.ufz.de/index.php?en=40126).
- The MPR comes with a [LICENSE][6] agreement, this includes also the GNU Lesser General Public License.

**Please note:** The [GitLab repository](https://git.ufz.de/chs/MPR) grants read access to the code.
If you would like to contribute to the code, please contact [mhm-admin@ufz.de](mailto:mhm-admin@ufz.de).

## Documentation

The online documentation for mHM can be found here (pdf versions are provided there as well):
- [latest](https://chs.pages.ufz.de/MPR/latest)
- [stable](https://chs.pages.ufz.de/MPR/stable)

## Dependencies and Requirements

* Fortran compiler: We support [gfortran](https://gcc.gnu.org/fortran/), [nagfor](https://www.nag.com/content/nag-fortran-compiler) and [ifort](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
* Build system: We support [make](https://www.gnu.org/software/make/) and [ninja](https://ninja-build.org/)
* [cmake](https://cmake.org/): Software for build automation
* [NetCDF-Fortran](https://github.com/Unidata/netcdf-fortran): NetCDF I/O for Fortran
* (optional) [fypp](https://github.com/aradi/fypp): Fortran pre-processor written in Python

It is recommended to have a clean installation at a custom location
for a C compiler, a Fortran compiler, the NetCDF C library and the
NetCDF Fortran library with consistent compilers.

We recommend to use a [conda](https://docs.conda.io/en/latest/) environment by using [Miniconda](https://docs.conda.io/en/latest/miniconda.html) to get all dependencies easily:

```bash
conda create -y --prefix ./mpr_env
conda activate ./mpr_env
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install -y cmake make fortran-compiler netcdf-fortran fypp f90nml python-graphviz
```

## Libraries bundled with MPR

### Doxygen Awesome

Doxygen Awesome is included as subrepo under `doc/doxygen-awesome-css/` from https://github.com/jothepro/doxygen-awesome-css. It is released under the MIT license.

### Cmake Fortran Scripts

The CHS Cmake Fortran Scripts repository is included as subrepo under `cmake/` from https://git.ufz.de/chs/cmake-fortran-scripts. It is released under the GNU LGPLv3 license.

### HPC Fortran Module Loads

The CHS HPC Fortran Module Loads repository is included as subrepo under `hpc-module-loads/` from https://git.ufz.de/chs/HPC-Fortran-module-loads. It is released under the GNU LGPLv3 license.


## Cite as

Please cite the MPR tool ([Schweppe et al. (2022)](https://doi.org/10.5194/gmd-15-859-2022)) and the MPR framework ([Samaniego et al. (2010)](https://doi.org/10.1029/2008WR007327), [Kumar et al. (2013](https://doi.org/10.1029/2012WR012195)) in your scientific work.

You can reference the code by its [Zenodo ID](http://doi.org/10.5281/zenodo.4650513).

[1]: https://git.ufz.de/chs/mpr/tree/1.0.4
[3]: doc/src/07_RELEASES.md
[4]: https://git.ufz.de/chs/mpr/tags/
[5]: https://chs.pages.ufz.de/MPR
[6]: LICENSE
[7]: doc/mpr_papers.md

