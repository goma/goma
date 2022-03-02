
## Build Instructions

Goma now uses `CMake` to build, Makefile build system has been removed

Common configuration options can be configured from the command line
or through the `CMake` GUI interfaces
 
    cmake -B <build directory> <path to goma src>
    # cmake -B <build directory> <path to goma src> -DMDE=10 -DMAX_PROB_VAR=23 ...more options
    make -C <build directory>

Install should be handled by setting CMAKE_INSTALL_PREFIX

Options that might be useful:

    DISABLE_COLOR_ERROR_PRINT=[ON|OFF]
    CHECKFINITE=[ON|OFF]
    FP_EXCEPT=[ON|OFF] # Floating point exceptions, faster than 
                         CHECKFINITE but will abort on failure
    ENABLE_UMFPACK=[ON|OFF]
    ENABLE_ARPACK=[ON|OFF]
    ENABLE_OMEGA_H=[ON|OFF]
    ENABLE_SPARSE=[ON|OFF]
    ENABLE_METIS=[ON|OFF]
    MDE=<number>
    MAX_CONC=<number>
    MAX_EXTERNAL_FIELD=<number>
    MAX_PROB_VAR=<number>

## Third party libraries

### Configure options

Goma will do its best to find third party libraries, hints can be given through cmake:

    -DTrilinos_DIR=<path to trilinos>/lib/cmake/Trilinos
    -DSEACASExodus_DIR=<path to seacas>/lib/cmake/SEACASExodus
    -DUMFPACK_DIR=<path to suitesparse install>
    -DSparse_PREFIX=<path to sparse install>
    -DARPACK_PREFIX=<path to arpack install>
    -DMETIS_PREFIX=<path to metis install>

PETSc is configured through `pkg-config` so `pkg-config` should be available on the command line and `PETSC_DIR` (and maybe `PETSC_ARCH` depending on build) should be set environment variables

### Building third party libraries

#### Spack

Third party libraries can be built with spack for development of goma using:

    spack env create goma <path to goma src>/spack.yaml
    spack env activate goma
    spack install

Currently expects the spack develop branch

Then Goma can be configured using that environment:

    spack env activate goma
    cmake -B <build directory> <path to goma src>
    make -C <build directory>

A spack package is in progress and soon one might be able to:

   spack install goma

#### Legacy build script

The legacy build script has been updated slightly with newer libraries.
This is missing Omega_h and PETSc support

We are considering this script deprecated in favor of the `spack` tool

See [scripts/README.md](scripts/README.md)

#### Docker build

A docker build with libraries pre-built is available at [https://hub.docker.com/r/westonortiz/goma-libs](https://hub.docker.com/r/westonortiz/goma-libs)

Source for the docker build is here: [https://github.com/wortiz/docker-goma-libs](https://github.com/wortiz/docker-goma-libs)



