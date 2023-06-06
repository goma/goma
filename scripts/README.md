# Goma dependencies build scripts

Before building Goma, a number of packages and libraries must be available.

* OpenMPI
* CMake
* HDF5
* MatIO
* NetCDF
* ARPACK
* BLAS
* LAPACK
* ParMETIS
* Sparse
* SuperLU_DIST
* Y12M
* ScaLAPACK
* MUMPS
* SuiteSparse
* Trilinos w/SEACAS

An example build script is available at `build-goma-dependencies.sh`

Dependencies for these scripts are:

* gcc
* g++
* gfortran
* m4
* zlib library
* libX11 library             

# Additional Requirements

## Environment Variables
OpenMPI should be added to the path and library path:

    export PATH="/[path to gomalibs]/openmpi-4.1.1/bin:$PATH"

SEACAS tools from Trilinos (e.g. aprepro and blot) should be added to
your path

    export PATH="/[path to gomalibs]/trilinos-13.0.1/bin:$PATH"  

An example configuration for bash will be written to

    [path to gomalibs]/config.sh

One should be able to build goma with (assuming a bash like shell):

	source [path to gomalibs]/config.sh
	cmake -B build-goma [path to goma src]
	make -C build-goma

See [BUILD.md](../BUILD.md) for more build options

## Example dependencies for common operating systems:

The build-goma-dependencies.sh script relies on several packages readily available in many repositories.

For Ubuntu this will install the necessary packages to run the script:

`sudo apt-get install git build-essential m4 zlib1g-dev libx11-dev gfortran pkg-config autoconf`

For CentOS / Fedora

`sudo [yum|dnf] install git patch gcc gcc-c++ gcc-gfortran m4 make wget bzip2 tar zlib-devel libX11-devel pkgconfig`

## Goma Dependencies build script usage

Name

	build-goma-dependencies.sh

Synopsis

	build-goma-dependencies.sh -jN [library install path]

Description

	This script attemps to download, compile and organize all the
	dependencies of goma.  They are downloaded into [library
	install path]/tars and extracted into [library install path].
	No attempt is made by the script to place any file outside of
	[library install path].

	The minumum requirements for this process to succeed are:
	    
	    - gcc
	    - gcc-g++
	    - gcc-gfortran
	    - m4
	    - zlib
	    - libX11

	A command such as 'yum install gcc-c++ gcc-gfortran m4
	zlib-static libX11-devel' would ensure these packages are
	available on an enterprise linux distribution (e.g. CentOS,
	Redhat, Fedora, etc.).

Options

        -jN  N : Number of make jobs to run.
                 Several packages do not possess sufficient Make
                 configurations to run in parallel, so these are made
                 with one make job.

