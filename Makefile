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
	cd src && ${MAKE}  $@  TARGET=goma  


tags:
	find . -name "*.[ch]" -exec etags -a '{}' \;

compile_commands.json:
	${MAKE} -C src realclean
	bear ${MAKE} -C src all

cppcheck: compile_commands.json
	cppcheck --enable=all --project=compile_commands.json

realclean:
	rm -f compile_commands.json
	cd src && ${MAKE} $@
guts:
	cd src && ${MAKE}  TARGET=goma_guts "DEFINES = -Dlinux -DCOMPILER_64BIT -DENABLE_AMESOS -DTRILINOS -DCHECK_FINITE -DMDE=27 -DCOUPLED_FILL -DPARALLEL -DHAVE_STRATIMIKOS -DEIGEN_SERIAL -DNO_CHEBYSHEV_PLEASE -DMAX_PROB_VAR=16 -DMAX_EXTERNAL_FIELD=4 -DMAX_CONC=4 "

ve:
	cd src && ${MAKE} TARGET=goma_ve "DEFINES = -Dlinux -DCOMPILER_64BIT -DENABLE_AMESOS -DTRILINOS -DCHECK_FINITE -DMDE=27 -DMAX_PROB_VAR=44 -DMAX_EXTERNAL_FIELD=6 -DMAX_CONC=8 -DCOUPLED_FILL -DPARALLEL -DHAVE_STRATIMIKOS -DEIGEN_SERIAL -DPHASE_COUPLED_FILL -DNEW_HOFFMAN_FCN_PLEASE -DALE_DCA_INFO_PLEASE -DRELAX_ON_TRANSIENT_PLEASE -DRESET_TRANSIENT_RELAXATION_PLEASE  -DALLOW_NEGATIVE_TIMES_PLEASE  -DDM_COORD_SCALE_PLEASE -DVAR_UPDATE_UNITY_SCALE"

