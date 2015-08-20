# Goma dependencies build scripts

These build scripts exist to alleviate the pain of building the
libraries goma needs to run.

There are 3 options to choose from:

`build-goma-dependencies.sh`: This script is known to work with the
most platforms and older compiler versions.

`build-goma-dep-trilinos-12.sh`: This script builds with the newest
trilinos and expects gcc 4.8 or gcc 4.9.

`build-goma-dep-trilinos-12-noc++11.sh`: This script builds with the
newest trilinos but disables support for c++11 to support older
compilers.

Dependencies for these scripts are:

* gcc
* g++
* gfortran
* m4
* zlib library
* libX11 library             

# Additional Requirements

openmpi should be added to your path and library path:

    export LD_LIBRARY_PATH="/[path to gomalibs]/openmpi-1.6.4/lib:$LD_LIBRARY_PATH"
    export PATH="/[path to gomalibs]/openmpi-1.6.4/bin:$PATH"

SEACAS tools from Trilinos (e.g. aprepro and blot) should be added to
your path

    export PATH="/[path to gomalibs]/trilinos-11.8.1-Build/bin:$PATH"  

Or for Trilinos 12:

    export PATH="/[path to gomalibs]/trilinos-12.0.1-Build/bin:$PATH"

After that adjust your settings.mk to point to the proper gomalibs
and Trilinos location (examples available in the main goma directory).

If you built Trilinos 12 with c++11 (not the noc++11 script) you also
need to declare the CXXSTD in your settings.mk:

    CXXSTD=-std=c++11


## Goma Dependencies build script usage

Name

	build-goma-dependencies.sh

Synopsis

	build-goma-dependincies.sh -jN [library install path]

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

