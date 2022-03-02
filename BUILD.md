
## Build Instructions

Goma now uses `CMake` to build, Makefile build system has been removed

Common configuration options can be configured from the command line
or through the `CMake` GUI interfaces

    cmake -B <build directory>
    # cmake -B <build directory> <path to goma> -DMDE=10 -DMAX_PROB_VAR=23 ...more options
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

### Third party libraries

Goma will do its best to find third party libraries, hints can be given through cmake:

    -DTrilinos_DIR=<path to trilinos>/lib/cmake/Trilinos
    -DSEACASExodus_DIR=<path to seacas>/lib/cmake/SEACASExodus
    -DUMFPACK_DIR=<path to suitesparse install>
    -DSparse_PREFIX=<path to sparse install>
    -DARPACK_PREFIX=<path to arpack install>
    -DMETIS_PREFIX=<path to metis install>

PETSc is configured through `pkg-config` so `pkg-config` should be available on the command line and `PETSC_DIR` (and maybe `PETSC_ARCH` depending on build) should be set environment variables