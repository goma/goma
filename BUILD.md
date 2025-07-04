## Build Instructions

Goma uses CMake to build it's source code. A number of third party libraries are required
to use goma. A script to build third party libraries is available in the `tpls` folder.

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

#### TPL build script

A third party library builder is available in the `tpls` folder named `tpls/install-tpls.py`. 
This tool should work with any semi-recent python3 (RedHat 8 or later)

Some basic build dependencies should already be present on your system:

For Ubuntu this will install the necessary packages to run the script:

    sudo apt-get install git build-essential m4 zlib1g-dev libx11-dev gfortran pkg-config autoconf

For CentOS / Fedora

 .  sudo dnf install git patch gcc gcc-c++ gcc-gfortran m4 make wget bzip2 tar zlib-devel libX11-devel pkgconfig

or using yum

    sudo yum install git patch gcc gcc-c++ gcc-gfortran m4 make wget bzip2 tar zlib-devel libX11-devel pkgconfig

X11 is used to build `blot` and is optional and should be autodetected by the build script if not present.

See [tpls/README.md] for more details but a simple example would be:

    /path/to/goma/tpls/install-tpls.py -j <num procs> /path/to/library/install

This will install to `/path/to/library/install`

After completed you would then do:

    source /path/to/library/install/config.sh
    cd /path/to/goma
    cmake -Bbuild
    cmake --build build -j <num proc>


#### Spack

Goma is available in spack though many users have found it difficult to use:

    spack install goma@release


#### Docker build

A docker build with libraries pre-built is available at [https://hub.docker.com/r/westonortiz/goma-libs](https://hub.docker.com/r/westonortiz/goma-libs)

Source for the docker build is here: [https://github.com/goma/goma/tree/main/tpls](https://github.com/goma/goma/tree/main/tpls)



