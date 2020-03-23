# Goma dependencies build scripts

Before building Goma, a number of packages and libraries must be available.

The build scripts in this directory exist to alleviate the pain of building-by-hand the libraries goma needs to run. They are used regularly in Ubuntu 14.04+ and CentOS 6+.

They will build the following packages and libraries 

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

`build-goma-dependencies.sh`: This script builds Trilinos 12.14.1
 with support for c++11 and expects gcc to be version 4.8.1 or greater.

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

    export LD_LIBRARY_PATH="/[path to gomalibs]/openmpi-4.0.2/lib:$LD_LIBRARY_PATH"
    export PATH="/[path to gomalibs]/openmpi-4.0.2/bin:$PATH"

SEACAS tools from Trilinos (e.g. aprepro and blot) should be added to
your path

    export PATH="/[path to gomalibs]/trilinos-12.18.1/bin:$PATH"  

## settings.mk

Copy settings.mk-example to settings.mk and adjust `GOMA_LIBS` to point to the location where the libraries were built. If the build script was not used, locations for all the libraries referenced in  (examples available in the main goma directory).

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

