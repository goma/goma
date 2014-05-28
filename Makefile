#/************************************************************************ *
#* Goma - Multiphysics finite element software                             *
#* Sandia National Laboratories                                            *
#*                                                                         *
#* Copyright (c) 2014 Sandia Corporation.                                  *
#*                                                                         *
#* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
#* the U.S. Government retains certain rights in this software.            *
#*                                                                         *
#* This software is distributed under the GNU General Public License.      *
#\************************************************************************/

# The main Makefile is located in the src directory

# Available targets are

# DEFAULT (all): makes the goma executable and libs
# debug: makes the gomad executable
# install: installs goma executable, libraries, and include files to PREFIX
# clean: removes the libraries and object files
# realclean: performs clean, and removes executables

# Make changes to settings.mk, see settings.mk-example

default: all

.DEFAULT:
	cd src && ${MAKE} $@

