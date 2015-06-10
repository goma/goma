# Goma Dependincies build script

Name

	build-goma-dependencies.sh

Synopsis

	build-goma-dependincies.sh -jN [library install path]

Description

	This script attemps to download, compile and organize all the dependencies of goma.
	They are downloaded into [library install path]/tars and extracted into [library install path].
	No attempt is made by the script to place any file outside of [library install path].

	The minumum requirements for this process to succeed are:
	    
	    - gcc
	    - gcc-g++
	    - gcc-gfortran
	    - zlib
	    - libX11

	A command such as 'yum install gcc-c++ gcc-gfortran zlib-static libX11-devel' would ensure these
	packages are available on an enterprise linux distribution (e.g. CentOS, Redhat, Fedora, etc.).

Options

        -jN  N : Number of make jobs to run.
                 Several packages do not possess sufficient Make configurations
                 to run in parallel, so these are made with one make job.

