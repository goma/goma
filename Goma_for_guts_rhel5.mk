# /************************************************************************\
#  * copyright (c) 2006 Sandia Corporation.                               *
#  *					         			*
#  * Under the terms of Contract DE-AC04-94AL85000, there is a            *
#  * non-exclusive license for use of this work by or on behalf of the    *
#  * U.S. Government. Export of this program may require a license from   *
#  * the United States Government.                                        *
#  *					         			*
#  * This software is the property of Sandia Corporation and discloses    *
#  * material protectable under copyright laws of the United States.      *
#  * Use, Duplication, or Disclosure is prohibited, unless allowed        *
#  * subject to the terms of a separate license agreement.                *
#  *					         			*
#  * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES          *
#  * DEPARTMENT OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR       *
#  * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY    *
#  * LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,    *
#  * OR USEFULNESS OF ANY INFORMATION, APPARATUS OR PROCESS DISCLOSED,    *
#  * OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.*
#  *					         			*
# \************************************************************************/
# Makefile for goma
# $Id: Goma_for_guts_rhel5.mk,v 5.2 2010-04-07 22:27:00 prschun Exp $
#
# Instructions:	Read the comments below and modify accordingly for your
#		platform.
#
# Note:		Don't forget to do a "make depend" on most machines.
#

# Select the ARCH and MACH option appropriate for your platform,
# This checked-in makefile is set for
# 880 LAN machines with Red Hat Enterprise 4 Linux:
#            ARCH = aix4
#            ARCH = hpux9
#            MACH = sgi
#            ARCH = irix6
            MACH = linux
            ARCH = rhel4
# Use this for Red Hat 9 and previous releases:
#            ARCH = linux 
#            MACH = sun
#            ARCH = sun5
#            MACH = dec
#            ARCH = tru64

#            PURE = purify -cache-dir="/tmp"
#            PURE = purify 
#            PURE = quantify -cache-dir="/tmp"
#PURE_PATH = /usr/local/rational/releases/purify.i386_linux2.2003a.06.00
#            PURE = $(PURE_PATH)/purify -cache-dir="/tmp"
            PURE =

          INSTALL = cp
#       INSTFLAGS =  

#=======================================================================
#=======================================================================
# Internal (engsci LAN in Blgd 880) Users:
#       
# External (to Sandia) user sites (e.g., CRMPC)
#       INSTALL_DIR is your base location for all goma-related software

        GOMA_LIBS = /home/goma/libraries
        GOMA_DIST = /home/goma/production/$(ARCH)
#      INSTALL_DIR = [to be set by external users]
      INSTALL_DIR = 

# Other-than-default compiler needed for CGM on linux
# Use:                      /usr/local/gcc-2.95.3/bin/gcc, g++ and g77
#
# On Sun & sgi (cc,CC,f77); on tru64 (cc,cxx,f77); on linux (gcc,g++,g77)
#         PGI_PATH = /usr/local/pgi_5.0/linux86/5.0/bin
#          MPI_TOP = /net/sass1809/opt/SUNWhpc
# This is for Red Hat Enterprise 4:
#          MPI_TOP = /usr/local/mpi/mpich/32Bit/1.2.7/gcc-3.4.3
#          MPI_TOP = /usr/local/mpi/sierra/32Bit/1.2.7/gcc-3.4.3
# This is for Red Hat Enterprise 5: 1/15/2010
           MPI_TOP = /home/sntools/extras/mpi/mpich-1.2.7p1-gcc-3.4.6-32Bit
# This is for Red Hat 9 and below:
#          MPI_TOP = /usr/local/mpi/mpich/1.2.5.2/gcc-3.2.2
               CC = $(PURE) $(MPI_TOP)/bin/mpicc
              CPP = $(PURE) $(MPI_TOP)/bin/mpiCC
               FC = $(PURE) $(MPI_TOP)/bin/mpif77

#               CC = $(PURE) gcc
#              CPP = $(PURE) g++
#               FC = $(PURE) g77
# Some Linux platforms may require this:
#              FC = $(PURE) f77 -fno-second-underscore

#              TOP = $(INSTALL_DIR)/goma
#         HOME_DIR = $(HOME)/.linux
#               LD = $(PURE) $(CC)
#     TRILINOS or CGM
              LD = $(PURE) $(CPP)
              LEX = lex
            SHELL = /bin/sh
             MAKE = make
#       MAKEFLAGS = b

          APREPRO = aprepro
         EXO1EXO2 = ex1ex2v2
            FASTQ = fastq
            GJOIN = gjoin
            GEN3D = gen3d
           GREPOS = grepos
            GROPE = grope
               RM = rm -f
              LEX = flex
             YACC = bison

               AR = ar
          ARFLAGS = srv
           RANLIB = touch
#sun4      RANLIB = ranlib
#
###################################################################
#
# LOCAL DEFINES SECTION
#
#  USE_CHEMKIN: If you want to use chemkin, you need to define this
#               compiler definition during the compile.
#  SENKIN_OUTPUT: If you want to use senkin postprocessor, define
#                 this along with the USE_CHEMKIN compiler directive.
#  PARALLEL:    If you want to run parallel applications, define
#               This. This works even for serial applications, if
#               you add the prefix "mpirun -np 1" to your command
#               line. Therefore, I have made this the default.
#  HAVE_FRONT:  Has been deprecated.
#  HAVE_AZTEC:  Has been deprecated.
#  MATRIX_DUMP: Define this if you want to enable the dumping of
#               matrices to disk.
#  TRILINOS:    Define this if you want to link to the new Trilinos
#              linear solver package, including its version of Aztec.
#
#  'platform:'  Several 'ifdef's on platform type still exist.
#               The following platforms must contain a -D'platform'
#               DEFINE on the compile line:
#                      SGI         -Dsgi
#                      Linux       -Dlinux
#                      Tru64       -Ddec_osf1
#                      HP          -Dhpux
#
#  NEW_PARSER	To make a version that uses the new parser, activate 
#  		the -DNEW_PARSER line below and comment out
#  		the blank LY_SRC= and LY_OBJ= lines way below.

           DEFINES = -Dlinux -DENABLE_AMESOS \
                     -DTRILINOS -DCHECK_FINITE \
                     -DNO_CHEBYSHEV_PLEASE \
                     -DMDE=27 \
                     -DMAX_PROB_VAR=15 \
                     -DMAX_EXTERNAL_FIELD=4 \
                     -DMAX_CONC=4 \
                     -DCOUPLED_FILL \
                     -DPARALLEL -DEIGEN_SERIAL

#                     -DMAX_CONC=1
#                     -DUSE_CGM
#                     -DEIGEN_SERIAL \
#                     -DEIGEN_PARALLEL \
#                     -DNO_CHEBYSHEV_PLEASE \
#
#                     -DNEW_PARSER
#                     -DSUN_SUNOS_57
#
#                     -D_RWSTD_COMPILE_INSTANTIATE \
#                     -DCOUPLED_FILL
#                     -DXCODE
#                     -DLIBRARY_MODE
#                     -DTRILINOS 
#                     -DUSE_CHEMKIN -DSENKIN_OUTPUT -DDEBUG_HKM
#
###########################################################################
#
#        CHEMKIN SECTION
#       -----------------
#
#
#ARCH_NAME=sol
ARCH_NAME=linux
#
CHEMKIN_DEV_PATH    = /home/goma/chemkin_dev
#
#
CHEMKIN_ARCH_PATH   = $(CHEMKIN_DEV_PATH)/arch/$(ARCH_NAME)
#
#
CHEMKINIII_LIB_PATH = $(CHEMKIN_ARCH_PATH)/lib
CHEMKINIII_LIBS     =  -L$(CHEMKINIII_LIB_PATH) -lsklib -llklib \
                         -ltranlib -lcklib
CHEMKINIII_LIBS_DEP = $(CHEMKINIII_LIB_PATH)/libsklib.a \
                      $(CHEMKINIII_LIB_PATH)/liblklib.a \
                      $(CHEMKINIII_LIB_PATH)/libtranlib.a \
                      $(CHEMKINIII_LIB_PATH)/libcklib.a
#
CPC_LIB_PATH = $(CHEMKINIII_LIB_PATH)
CPC_INC_PATH = $(CHEMKIN_ARCH_PATH)/include
CPC_LIBS     = -L$(CPC_LIB_PATH) -lcpc -lckwrapper
CPC_LIBS_DEP = $(CPC_LIB_PATH)/libcpc.a \
               $(CPC_LIB_PATH)/libckwrapper.a
#
#  Default is to not link chemkin into Goma.
#  This is done by zeroing out these defines below
#  If you want to include the chemkin libraries, take out '#chm'
#  in the lines below and added USE_CHEMKIN to the list of
#  preprocessor defines.
#
CHEMKIN_INC  = # -I$(CPC_INC_PATH)
CHEMKIN_DEPS = # $(CPC_LIBS_DEP) $(CHEMKINIII_LIBS_DEP)
CHEMKIN_LIBS = # $(CPC_LIBS) $(CHEMKINIII_LIBS)
#
###################################################################
#
#  MPI SECTION
#               Should work without any knowledge of MPI...
#  NO MPI
#     MPI_LIB = 
#
##  mpi version (1.2.*)
#linux IP version, RH9 (desktop)
#         MPI_TOP = /usr/local/mpi/mpich/1.2.5.2/gcc-3.2.2
#linux GM version, RH9 (cluster)
#         MPI_TOP = /usr/local/mpi/mpich/1.2.5..10-gm-2.0.9-2.4.20-28.9smp.gpfs-i686/gcc-3.2.2
#dec      MPI_TOP = /usr/local/mpi/mpich-1.2.1
#sgi      MPI_TOP = /usr/lib64/mips4
#sun      MPI_TOP = /net/sass1809/usr/local/mpi/mpich-1.2.4
#
## path to include/mpi.h
#sgi      
#      MPI_INC = -I/usr/include
#sun/dec/linux
      MPI_INC = -I$(MPI_TOP)/include
#
## path to mpi libraries
#sgi
#      MPI_LIB = -L$(MPI_TOP) -lmpi
#dec/linux
      MPI_LIB = -L$(MPI_TOP)/lib -lmpe -lmpich -lpmpich
#sun
#      MPI_LIB = -L$(MPI_TOP)/lib -lmpe -lmpich -lpmpich -lthread -lsocket
#
#
############################################################################
#
#  ARPACK SECTION
#    To use ARPACK, either EIGEN_SERIAL or EIGEN_PARALLEL must be defined above.
#    Default is not to link in ARPACK.
#
      ARPACK_TOP = $(GOMA_LIBS)/ARPACK
#
#    Default: no arpack
#      ARPACK_LIB =
#
#    EIGEN_SERIAL: arpack only
      ARPACK_LIB = -L$(ARPACK_TOP)/lib -larpack_$(MACH)
#
#    EIGEN_PARALLEL: parpack and arpack
#      ARPACK_LIB = -L$(ARPACK_TOP)/lib -lparpack_MPI_$(MACH) -larpack_$(MACH)
#
############################################################################
#
# SEAMS LIBRARIES
#
#    The engsci lan uses the SEAMS distribution defined by the link to
#    /usr/local/eng_sci/current
#
#        ==================================================
# ==>NETCDF
#
#       NETCDF_INC need not be defined if EXODUSII_INC is defined below.
#
# For engsci LAN:
#       SEAMS_DIR  = /usr/local/eng_sci/struct/i686/current-gcc
#      For RHEL5
        SEAMS_DIR  = /home/goma/production/linux/SEAMS_2006-03-16-32-GCC
#
# For external customers:
#       SEAMS_DIR = $(GOMA_DIST)/SEAMS
       NETCDF_LIB = -L$(SEAMS_DIR)/lib -lnetcdf
#
#  Users External to Sandia (e.g., CRMPC )
#
#  NB: use SEAMS_64 in place of SEAMS for 64-bit sgi
#       NETCDF_LIB = -L$(GOMA_DIST)/SEAMS/lib  -lnetcdf
#       NETCDF_LIB = -L$(INSTALL_DIR)/SEAMS/lib  -lnetcdf
#
#        ==================================================
# ==>EXODUS LIBRARY 
#
       EXODUSII_INC = -I$(SEAMS_DIR)/inc
#
       EXODUSII_LIB = -L$(SEAMS_DIR)/lib -lexoIIv2c
#
#  Users External to Sandia (e.g., CRMPC )
#
#  NB: use SEAMS_64 in place of SEAMS for 64-bit sgi
#       EXODUSII_INC = -I$(GOMA_DIST)/SEAMS/inc
#       EXODUSII_LIB = -L$(GOMA_DIST)/SEAMS/lib  -lexoIIv2c
#       EXODUSII_INC = -I$(INSTALL_DIR)/SEAMS/inc
#       EXODUSII_LIB = -L$(INSTALL_DIR)/SEAMS/lib  -lexoIIv2c
#
# ==> CHACO LIBRARY
#
        CHACO_LIB = 
#        CHACO_LIB = -L$(SEAMS_DIR)/lib -lchaco
#
#  Users External to Sandia (e.g., CRMPC )
#
#  NB: use SEAMS_64 in place of SEAMS for 64-bit sgi
#       CHACO_LIB = -L$(GOMA_DIST)/SEAMS/lib  -lchaco
#       CHACO_LIB = -L$(INSTALL_DIR)/SEAMS/lib  -lchaco
#

#
############################################################################
#
# AZTEC LIBRARY
#
#   As of 7/8/2002, Aztec version 2.0 and all previous versions are no
#   longer supported for Goma. Version 2.1 or later is now a requirement
#   to build Goma.  Once the object-oriented Aztec release (aztecoo)
#   is made available in Trilinos, it can be used by defining TRILINOS
#   above and using INC and LIB lines which will be provided at that time.
#   Until then, Aztec must be taken from the version 2.1 distribution.
#   The AZTEC_2 and AZTEC_2_1 define flags are no longer necessary and
#   have been expunged from Goma.
#
#   NOTE: If ARPACK is linked in, it must use the LAPACK that comes in its
#         distribution.  This is taken care of in the Version 2.1 lines.
#
# LOCA Hopf tracking requires the Aztec2.1 Komplex library
# and "KOMPLEX" must be defined above.
# Currently, this is included with Trilinos.
#       KOMPLEX_TOP = $(GOMA_LIBS)/trilinos
#
# Default: Do not include KOMPLEX
#       KOMPLEX_INC =
#       KOMPLEX_LIB = -lkomplex
        KOMPLEX_LIB =
#
# Not using Trilinos? Use these lines to include KOMPLEX:
#       KOMPLEX_INC = -I$(KOMPLEX_TOP)/packages/komplex/src
# Serial version:
#       KOMPLEX_LIB = -L$(KOMPLEX_TOP)/lib -lkomplex
# Parallel version:
#       KOMPLEX_LIB = -L$(KOMPLEX_TOP)/lib-mpi -lkomplex
#
# Version 1.1 is now deprecated.
#
# Version 2.0 is now deprecated.
#
# Version 2.1
#       AZTEC_INC = -I/usr/local/Aztec-2.1/include $(KOMPLEX_INC)
# Serial version:
#      AZTEC_LIB = -L/usr/local/Aztec-2.1/lib -laztec -ly12m \
#                   $(ARPACK_LIB) -llapack $(KOMPLEX_LIB)
# Parallel version:
#       AZTEC_LIB = -L/usr/local/Aztec-2.1/lib -laztec-mpi -ly12m \
#                   $(ARPACK_LIB) -llapack $(KOMPLEX_LIB)
#

#
# TRILINOS Version -- now available!
#   Current path
#    TRILINOS_DIR = $(GOMA_DIST)/trilinos-4.0
    TRILINOS_DIR = $(GOMA_DIST)/trilinos-6.0.14
#   Several ARCH options are available - choose yours.
#    TRILINOS_ARCH = LINUX_RHE_SERIAL
#    TRILINOS_ARCH = LINUX_RHE4_MPI
    TRILINOS_ARCH = RHEL4_PARALLEL
    TRILINOS_INC  = $(TRILINOS_DIR)/$(TRILINOS_ARCH)/include
    TRILINOS_LIB  = $(TRILINOS_DIR)/$(TRILINOS_ARCH)/lib
#
       AZTEC_INC  = -I$(TRILINOS_INC) -I$(TRILINOS_DIR)/include
       AZTEC_LIB  = -L$(TRILINOS_LIB) -laztecoo -lifpack -ly12m \
                     $(ARPACK_LIB) -llapack $(KOMPLEX_LIB)
#      AZTEC_LIB  = -L$(TRILINOS_DIR)/lib-mpi -laztecoo \
#                     -lifpack -llapack -ly12m -lspblas
#
# Aztec version 1.1 is now deprecated.
#
# Aztec version 2.0 is now deprecated.
#
# Aztec version 2.1R (still available but must define AZTEC_2_1 vice TRILINOS):
#       AZTEC_INC = -I$(GOMA_DIST)/Aztec-2.1/include
#       AZTEC_LIB = -L$(GOMA_DIST)/Aztec-2.1/lib  -laztec
#       AZTEC_LIB = -L$(GOMA_DIST)/Aztec-2.1/lib  -laztec-mpi
#       AZTEC_INC = -I$(INSTALL_DIR)/Aztec-2.1/include
#       AZTEC_LIB = -L$(INSTALL_DIR)/Aztec-2.1/lib  -laztec
#       AZTEC_LIB = -L$(INSTALL_DIR)/Aztec-2.1/lib  -laztec-mpi
#
# Amesos solver package
#	AMESOS_LIB = 
	AMESOS_LIB =  -lamesos  -lepetraext -lepetra  -lteuchos
#
#######################################################################
#
# SPARSE LIBRARY
#
#      SPARSE_LIB = -lsparse
#
#  Users External to Sandia (e.g., CRMPC )
#
      SPARSE_LIB = -L$(GOMA_DIST)/sparse/lib  -lsparse
#      SPARSE_LIB = -L$(INSTALL_DIR)/sparse/lib  -lsparse
#
##############################################################################
#
# FRONTAL SOLVER LIBRARY
#
# The default is to include the frontal solver in the default compile.
# To exclude the frontal solver, take out the HAVE_FRONT define from
# the local defines section above and then comment out all of the
# entries below.
#
#      FRONT_LIB =  -L$(GOMA_LIBS)/front/lib -lmpfront_$(ARCH) -lparlib_$(ARCH)
#
#     FRONT_DEP = $(FRONT_TOP)/snl_parlib/lib_snl_parlib_sun5.5.a \
#                    $(FRONT_TOP)/snl_mpfront/lib_snl_mpfront_sun5.5.a
#
#  Users External to Sandia (e.g., CRMPC )
#
      FRONT_LIB = -L$(GOMA_DIST)/front/lib  -lmpfront_$(MACH) -lparlib_$(MACH)
#      FRONT_LIB = -L$(INSTALL_DIR)/front/lib  -lmpfront -lparlib
#FRONT_LIB = $(INSTALL_DIR)/front/snl_mpfront/*.o  $(INSTALL_DIR)/front/snl_parlib/*.o
# FRONT_LIB = $(GOMA_LIBS)/front/snl_mpfront/Objects_linux/*.o  $(GOMA_LIBS)/front/snl_parlib/Objects_linux/*.o
#
##############################################################################
#
# UMFPACK SOLVER LIBRARY
#
#      UMFPACK_INC = -I$(GOMA_LIBS)/umfpack/include
#      UMFPACK_LIB = -L$(GOMA_LIBS)/umfpack/lib  -lumfpack_$(ARCH) -lamd_$(ARCH)  
#
#
#  Users External to Sandia (e.g., CRMPC )
#
      UMFPACK_INC = -I$(GOMA_DIST)/umfpack/include
      UMFPACK_LIB = -L$(GOMA_DIST)/umfpack/lib  -lumfpack  -lamd
#      UMFPACK_INC = -I$(INSTALL_DIR)/umfpack/include
#      UMFPACK_LIB = -L$(INSTALL_DIR)/umfpack/lib  -lumfpack  -lamd
#
##############################################################################
#
# BLAS and LAPACK  LIBRARIES
# dec:     no specification required (included in dxml & cxml libs)
# sun:     no specification required (included in sunperf lib)
# sgi:     -L/usr/lib64/mips4 -lblas
# linux:   -L/usr/lib -llapack & -L/usr/lib -lblas

       LAPACK_LIB =
         BLAS_LIB =
#       LAPACK_LIB = -L/usr/lib -llapack
#         BLAS_LIB = -L/usr/lib -lblas
#
##############################################################################
#
#  LIBRARY for brkfix for linking these objects directly to Goma
       BRK_LIB = 
#       BRK_LIB = -L$(GOMA_DIST)/brkfix/lib -lbrk
#
##############################################################################
#
# OTHER SOLVER LIBRARIES
#
         Y12M_LIB = 
      LINPACK_LIB = 
      HARWELL_LIB =
#     SUPERLU_LIB = 
      SUPERLU_LIB =  -L$(GOMA_DIST)/SuperLU/2.0/build_RHE4 -lsuperludist
#      HARWELL_LIB = -L/home/pasacki/.sun5/lib -lma28
#      HARWELL_LIB = -L/usr/local/lib -lharwell
#
##############################################################################
#
#  CGM (Common Geometry Module) LIBRARY
#
#  Don't forget to set up the correct linker ($LD), above, and the
#  correct linker options $(LD_FLAGS), below.  These settings are
#  required in your .cshrc because they are accessed at run-time, not
# just for compile-time. 
#
#  The easiest way to get LD_LIBRARY_PATH and other variables defined
#  properly here is to have the following sequence in your .cshrc:
#    setenv CUBIT_STEP_PATH Path to location ot the STEP library
#    setenv ACIS_DIR = root directory for ACIS
#    setenv ACIS_ARCH = ACIS' name for the machine architecture you're on
#    setenv ACIS_LIB_DIR $(ACIS_DIR)/lib/$(ACIS_ARCH)_so
#    setenv LD_LIBRARY_PATH $(LD_LIBRARY_PATH):$(ACIS_LIB_DIR)
#
#  Specifically for engsci LAN on the following platforms:
#  For sass2889 (Solaris 2.8):
#    setenv CUBIT_STEP_PATH /net/sass2248/usr/local/eng_sci/cubit/acis/acis7.0.5/step/tools/solaris
#    setenv ACIS_DIR /net/sass2248/usr/local/eng_sci/cubit/acis/acis7.0.5
#    setenv ACIS_ARCH solaris
#  For sada290 (Tru64):
#    setenv CUBIT_STEP_PATH /usr/local/eng_sci/cubit/acis/acis7.0.5/step/tools/osf
#    setenv ACIS_DIR /usr/local/eng_sci/cubit/acis/acis7.0.5
#    setenv ACIS_ARCH osf1_64
#  For linux desktops:
#    setenv CUBIT_STEP_PATH /usr/local/eng_sci/cubit/acis/acis7.0.5/step/tools/linux
#    setenv ACIS_DIR /usr/local/eng_sci/cubit/acis/acis7.0.5
#    setenv ACIS_ARCH linux
#For all platforms
#    setenv ACIS_LIB_DIR ${ACIS_DIR}/lib/${ACIS_ARCH}_so
#    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ACIS_LIB_DIR}
#
         ACIS_INC = -I${ACIS_DIR} -I${ACIS_DIR}/base \
           -I${ACIS_DIR}/bool -I${ACIS_DIR}/br -I${ACIS_DIR}/covr \
           -I${ACIS_DIR}/cstr -I${ACIS_DIR}/fct -I${ACIS_DIR}/ga \
           -I${ACIS_DIR}/gi -I${ACIS_DIR}/intr -I${ACIS_DIR}/kern \
           -I${ACIS_DIR}/oper -I${ACIS_DIR}/part -I${ACIS_DIR}/rbase \
           -I${ACIS_DIR}/scm -I${ACIS_DIR}/swp -I${ACIS_DIR}/ofst \
           -I${ACIS_DIR}/law -I${ACIS_DIR}/clr
         ACIS_LIB = -L$(ACIS_LIB_DIR) -lhealhusk -lblend -ligeshusk \
           -lstephusk -lcaselib -llop_husk -lsweep -lsbool -lskin \
           -lcover -loffset -ltransutl -loperator -lrem_husk \
           -lrbi_husk -llopt_husk -lboolean -lskin -lga_husk \
           -lrnd_husk -lct_husk -leuler -lconstrct -lfaceter \
           -lintersct -lkernel -llawutil \
           -lbaseutil -lXt -lX11 -liostream 
          CGM_ROOT = $(GOMA_DIST)/cgm
          CGM_INC = -I$(CGM_ROOT)/util -I$(CGM_ROOT)/geom \
            -I$(CGM_ROOT)/list -I$(CGM_ROOT)/geom/ACIS \
            -I$(CGM_ROOT)/geom/virtual
          CGM_LIB = -L${CGM_ROOT}/geom -lcubit_geom \
            -L${CGM_ROOT}/util -lcubit_util \
            -L${CGM_ROOT}/list -lcubit_list \
            -L${CGM_ROOT}/geom/virtual -lcubit_virtual \
            -L${CGM_ROOT}/geom/facet -lcubit_facet \
            -L${CGM_ROOT}/geom/facetbool -lcubit_facetboolstub \
            -L${CGM_ROOT}/geom -lcubit_geom \
            -L${CGM_ROOT}/geom/ACIS -lcubit_ACIS \
            -L${CGM_ROOT}/list -lcubit_list \
            -L${CGM_ROOT}/util -lcubit_util
# Uncomment the following to clobber the above if you don't want
# CGM/ACIS libraries:
        ACIS_INC =
        ACIS_LIB =
         CGM_INC =
         CGM_LIB =

#
############################################################################
#
#          SEACAS = /usr/local/eng_sci/struct
#          SEACAS = /usr/local
#
             _LIB = 
#            _LIB = -L/usr/local/lib -L/usr/lib
# This specification required for 64-bit Libraries for SGI
#            _LIB = -L/usr/lib64
#
############################################################################
#
# FORTRAN LIBRARIES
#
# SUN/Solaris:      -L/net/sass1809/opt/SUNWspro/lib  -lM77 -lF77 -lsunmath -lsunperf -lfsu
# SGI/IRIX64:       -lftn
# HP/HPUX:          -lf
# RH/Linux 6.1/7.x: -lc -lgcc -lg2c
# RH/Linux 9.0:     -lc  -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.2 -lgcc -lg2c
# IBM/AIX4.1.3:     -lxlf90 -lxlf -lm
# DEC/OSF1/TRU64:   -lfor -lFutil -lUfor
#
# ----------------------------------------------------------------------
# For Red Hat Enterprise 4:
#      FORTRAN_LIB =  -lc  -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 -lgcc -lg2c
      FORTRAN_LIB =  -lc  -L/usr/local/gcc/32Bit/3.4.3/lib -lgcc -lg2c -lstdc++
# For Red Hat 9 and below:
#      FORTRAN_LIB =  -lc  -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 -lgcc -lg2c
#
############################################################################
#
# SYSTEM LIBRARIES
#
# SUN/Solaris:       -lm -lrpcsvc -lnsl -lposix4
# SGI/IRIX64:        -lm
# HP/HPUX:           -lm
# RH/Linux 6.1/7.x:  -lm
# RH/Linux 9.0:      -lm
# IBM/AIX4.1.3:      -lm
# DEC/OSF1:          -lm -lrpcsvc -lots -lots2
# CPQ/DEC/TRU64:     -lm -ldxml -lcxml
#
# ----------------------------------------------------------------------
# Non-CGM:
          SYS_LIB = -lm
# CGM: (For some reason this is now necessary for the C++ to get the C
# exit_code() routine.  It doesn't hurt anything so I'm making it the
# default. 
#          SYS_LIB = -lm libCrun.so.1
# for sun, libCrun.so.1 leads to loader warning for non-CGM builds.
# dal note: above disabled 1103 since linux is default system
#
############################################################################
#
# GOMA and SPARSE take include files from current working directory...
#
         GOMA_INC =
#
       SPARSE_INC =
#
#######################################################################
#
        _INC = -I.

        INCLUDES = $(_INC) $(GOMA_INC) $(SPARSE_INC) \
                   $(AZTEC_INC) $(MPI_INC) $(EXODUSII_INC) $(NETCDF_INC) \
                   $(CHEMKIN_INC) $(CGM_INC) $(ACIS_INC) $(UMFPACK_INC)

#######################################################################
#
#  COMPILER FLAG SECTION
#
# Select flags appropriate for platform and purpose; enter below ---- line.
# Flags used in building software for various platforms at SNL shown below.
#
# SUN/Solaris (SNL Default)
# 
#  CCFLAGS =    -xO3   optimized; do not use more than level 3 for production
#          =    -msupersparc -O3    use with gcc
#          =    -g -fnonstd         -f nonstd makes NaNs/etc. crash.
#          =    -g -pg              for profiling 
#
#          =    -fast <=== RRL Note (8/25/98) can cause core dumps on
#                          otherwise working code; recommend it never be used
#          =    -v                More lint-like information
#
# SGI/IRIX64
#          =    -64 -O2 -OPT:Olimit=5200 -common -woff 1174,1552
# HP/HPUX
#          =    -Ae +O1 (or +O2)
# RH/Linux 6.1/7.x
#          =    -O
# IBM/AIX4.1.3
#          =    -O3 -qstrict
# DEC/OSF1
#          =    -O3 -migrate -ansi_alias -N32
#
# CPQ/Tru64
#          =    -O -std1
#
# profiling
#          =    -g -pg
#
# ALL PLATFORMS
#          =    -g   debug (can use '-O -g' on linux to debug optimized code)
#
#
# FFLAGS & LD_FLAGS
#          =    $(CCFLAGS)    for SUN and Linux
#
#          =    -O -64 & -64  for SGI/IRIX64
#
#          =    +O1  &  $(CCFLAGS)  for HP/HPUX
# ----------------------------------------------------------------------------
#
          CCFLAGS = -O
#	  CCFLAGS = -g 
#
           FFLAGS = $(CCFLAGS) $(DEFINES) $(INCLUDES)
#          FFLAGS = -64 -O2 $(DEFINES) $(INCLUDES)
#
           CFLAGS = $(CCFLAGS) $(DEFINES) $(INCLUDES)
         CPPFLAGS = $(CCFLAGS) $(DEFINES) $(INCLUDES) -Wno-deprecated
#cgm_sun         CPPFLAGS = $(CCFLAGS) $(DEFINES) $(INCLUDES) -DCUBIT_ACIS_VERSION=705 -instances=static -z muldefs

# If you're using the old non-C++ linker, then this is probably OK:
         LD_FLAGS = $(CCFLAGS)
#         LD_FLAGS = $(FFLAGS)
# Otherwise LD = CC (or some other C++ compiler), and special flags
# are necessary:
#sun         LD_FLAGS = $(CPPFLAGS) -instances=static -z muldefs
#other       LD_FLAGS = $(CPPFLAGS)
#
##############################################################################
#
#LIBRARY MACRO
#
             LIBS = $(_LIB) \
                    $(FRONT_LIB) \
                    $(SPARSE_LIB) \
                    $(LAPACK_LIB) \
                    $(Y12M_LIB) \
                    $(LINPACK_LIB) \
                    $(AZTEC_LIB) \
                    $(AMESOS_LIB) \
                    $(SUPERLU_LIB) \
                    $(ARPACK_LIB) \
                    $(HARWELL_LIB) \
                    $(UMFPACK_LIB) \
                    $(BLAS_LIB) \
                    $(EXODUSII_LIB) \
                    $(NETCDF_LIB) \
                    $(CGM_LIB) \
                    $(ACIS_LIB) \
                    $(BRK_LIB) \
                    $(CHACO_LIB) 
                   

         END_LIBS = \
                    $(MPI_LIB) \
                    $(FORTRAN_LIB) \
                    $(SYS_LIB)



#
# The following is a library dependence macro
#
          LIB_DEP = $(FRONT_DEP)
#
##############################################################################
##############################################################################
#
# old stuff down here ...
#
# IBM AIX 3.2 (spnode0?)
# LIBS	=	-L/usr/local/lib -lsparse -lharwell -lexoIIv2c -lnetcdf \
#                 -lxlf90 -lxlf -lm -bI:/usr/lpp/xlf/lib/lowsys.exp
#
################################################################
#
# The binaries and libraries will be difft for each platform, but the
# include files will be the same.

LIBDIR = $(TOP)/lib

BINDIR = $(TOP)/bin

INCDIR = $(TOP)/include

G_LIB = libgoma.a

            U_LIB = libgomau.a

           TARGET = goma

U_LIB = libgomau.a

TARGET = goma
#
#######################################################################
#######################################################################
#
#______ Users should not normally need to modify below this line_______
#
#######################################################################
#######################################################################
#
# _____ Element routines "el_" prefix _______________________________________

EL_SRC=	el_elm_info.c el_quality.c exo_conn.c

EL_OBJ=	el_elm_info.o el_quality.o exo_conn.o

EL_INC=	el_elm.h\
        el_elm_info.h\
        el_geom.h\
        el_quality.h\
        exo_conn.h\
        exo_struct.h


MAIN_SRC = main.c

MAIN_OBJ = main.o

MAIN_INC = std.h\
           goma.h

# _____ Geometry Model routines "gm_" prefix _______________________

GM_SRC = gm_cgm_c_interface.C \
         gm_cgm_util.c

GM_OBJ = gm_cgm_c_interface.o \
         gm_cgm_util.o

GM_INC = gm_cgm_c_interface.h \
         gm_cgm_util.h \
         gm_cgm_typedefs.h

# _____ Advanced capabilities routines "ac_" prefix _______________________

AC_SRC = ac_conti.c\
         ac_hunt.c\
         ac_update_parameter.c\
         ac_display_parameter.c\
         ac_loca_interface.c\
         ac_stability.c\
         ac_stability_util.c\
         ac_particles.c

AC_OBJ = ac_conti.o\
         ac_hunt.o\
         ac_update_parameter.o\
         ac_display_parameter.o\
         ac_loca_interface.o\
         ac_stability.o\
         ac_stability_util.o\
         ac_particles.o

AC_INC = ac_conti.h\
         ac_hunt.h\
         ac_update_parameter.h\
         ac_stability.h\
         ac_stability_util.h\
         ac_particles.h

# _____ Library of continuation algorithms "loca_" prefix ____________________

LOCA_SRC = loca_lib.c\
           loca_bord.c\
           loca_util.c\
           loca_eigenvalue.c \
           loca_eigen_c2f.F \
           loca_eigen_cayley.F

LOCA_OBJ = loca_lib.o\
           loca_bord.o\
           loca_util.o\
           loca_eigenvalue.o \
           loca_eigen_c2f.o \
           loca_eigen_cayley.o

LOCA_INC = loca_const.h\
           loca_util_const.h

# _____	Reacting flow routines "rf_" prefix __________________________________

RF_SRC= rf_allo.c\
        rd_dpi.c\
        rf_element_storage.c\
        rd_exo.c\
        rd_mesh.c\
        rf_node.c\
        rf_node_vars.c\
        rf_pre_proc.c\
        rf_setup_problem.c\
        rf_shape.c\
        rf_solve.c\
        rf_util.c\
        rf_vars.c\
        wr_dpi.c\
        wr_side_data.c\
        wr_soln.c\
        wr_exo.c

RF_OBJ= rf_allo.o\
        rd_dpi.o\
        rf_element_storage.o\
        rd_exo.o\
        rd_mesh.o\
        rf_node.o\
        rf_node_vars.o\
        rf_pre_proc.o\
        rf_setup_problem.o\
        rf_shape.o\
        rf_solve.o\
        rf_util.o\
        rf_vars.o\
        wr_dpi.o\
        wr_side_data.o\
        wr_soln.o\
        wr_exo.o

RF_INC= rf_allo.h\
        rf_bc.h\
        rf_bc_const.h\
        rf_element_storage_const.h\
        rf_element_storage_struct.h\
        rf_fem.h\
        rf_fem_const.h\
        rf_fill_const.h\
        rf_io.h\
        rf_io_const.h\
        rf_io_defn.h\
        rf_io_structs.h\
        rf_masks.h\
        rf_masks_const.h\
        rf_masks_globs.h\
        rf_mp.h\
        rf_node.h\
        rf_node_const.h\
        rf_pre_proc.h\
        rf_shape.h\
        rf_solve.h\
        rf_solver.h\
        rf_solver_const.h\
        rf_util.h\
        rf_vars_const.h \
        rf_vars_defn.h\
        rd_dpi.h\
        rd_exo.h\
        rd_mesh.h\
        wr_dpi.h\
        wr_exo.h\
        wr_side_data.h\
        wr_soln.h


# _____ Solver routines "sl_" prefix __________________________________________

SL_SRC= sl_lu.c\
        sl_lustat.c\
        sl_lu_fill.c\
        sl_ma28.c\
        sl_matrix_dump.c\
        sl_matrix_util.c\
        sl_squash.c\
        sl_util.c\
        sl_aux.c\
        sl_auxutil.c\
        sl_umf.c\
        sl_front_setup.c\
        sl_eggrollwrap.c\
        sl_eggrollutil.c\
        sl_eggroll01.c\
        sl_eggroll02.c\
        sl_eggroll03.c\
        sl_eggroll04.c\
        sl_eggroll05.c\
        sl_amesos_interface.C

SL_OBJ= sl_lu.o\
        sl_lustat.o\
        sl_lu_fill.o\
        sl_ma28.o\
        sl_matrix_dump.o\
        sl_matrix_util.o\
        sl_squash.o\
        sl_util.o\
        sl_aux.o\
        sl_auxutil.o\
        sl_umf.o\
        sl_front_setup.o\
        sl_eggrollwrap.o\
        sl_eggrollutil.o\
        sl_eggroll01.o\
        sl_eggroll02.o\
        sl_eggroll03.o\
        sl_eggroll04.o\
        sl_eggroll05.o\
        sl_amesos_interface.o

SL_INC= spConfig.h\
        spDefs.h\
        spmatrix.h\
        sl_aux.h\
        sl_auxutil.h\
        sl_eggroll.h\
        sl_eggroll_def.h\
        sl_lu.h\
        sl_matrix_util.h\
        sl_umf.h\
        sl_util.h\
        sl_util_structs.h \
        sl_rcm.h


# _____	Machine dependent routines "md_" prefix ______________________________

MD_SRC=	md_timer.c \
        md_ieee.c

MD_OBJ=	md_timer.o \
        md_ieee.o

MD_INC= md_timer.h

# _____ Distributed processing routines "dp_" prefix __________________________

DP_SRC= dp_comm.c\
        dp_map_comm_vec.c\
        dp_utils.c\
        dp_vif.c

DP_OBJ= dp_comm.o\
        dp_map_comm_vec.o\
        dp_utils.o\
        dp_vif.o

DP_INC= dp_comm.h\
        dp_map_comm_vec.h\
        dp_types.h\
        dp_utils.h\
        dp_vif.h\
        dpi.h


# _____	Moving mesh routines "mm_" prefix _____________________________________


BC_SRC= bc_colloc.c\
        bc_curve.c\
        bc_dirich.c\
        bc_integ.c\
        bc_rotate.c\
        bc_special.c\
        bc_contact.c \
        bc_surfacedomain.c

BC_OBJ= bc_colloc.o\
        bc_curve.o\
        bc_dirich.o\
        bc_integ.o\
        bc_rotate.o\
        bc_special.o\
        bc_contact.o \
        bc_surfacedomain.o

BC_INC= bc_colloc.h\
        bc_curve.h\
        bc_dirich.h\
        bc_integ.h\
        bc_rotate.h\
        bc_special.h\
        bc_contact.h \
        bc_surfacedomain.h


# _____ Moving mesh routines "mm_" prefix _____________________________________

MM_SRC= mm_as_alloc.c\
        mm_augc_util.c\
        mm_bc.c\
        mm_bc_conflict.c\
        mm_chemkin.c\
        mm_eh.c\
        mm_fill.c\
        mm_fill_aux.c\
        mm_fill_fill.c\
	mm_fill_ls.c \
        mm_fill_interface.c\
        mm_fill_jac.c\
        mm_fill_porous.c\
	mm_porous_EOS.c\
        mm_fill_potential.c\
        mm_fill_pthings.c\
        mm_fill_ptrs.c\
        mm_fill_rs.c\
        mm_fill_shell.c\
        mm_fill_solid.c\
        mm_fill_species.c\
        mm_fill_stress.c\
        mm_fill_terms.c\
        mm_fill_common.c\
        mm_fill_util.c\
        mm_flux.c\
        mm_input.c\
        mm_input_bc.c\
        mm_input_mp.c\
        mm_input_particles.c\
        mm_input_util.c\
        mm_interface.c\
        mm_matrl.c\
        mm_more_utils.c\
        mm_ns_bc.c\
        mm_shell_bc.c\
        mm_numjac.c\
        mm_placid.c\
        mm_post_proc.c\
        mm_post_proc_util.c\
        mm_prob_def.c\
        mm_propertyJac.c\
        mm_qp_storage.c\
        mm_qtensor_model.c\
        mm_shell_util.c\
        mm_sol_nonlinear.c\
        mm_species.c\
        mm_std_models.c \
        mm_unknown_map.c\
        mm_viscosity.c \
        mm_dil_viscosity.c

MM_OBJ= mm_as_alloc.o\
        mm_bc.o\
        mm_bc_conflict.o\
        mm_chemkin.o\
        mm_eh.o\
        mm_fill.o\
        mm_fill_aux.o\
        mm_fill_fill.o\
	mm_fill_ls.o\
        mm_fill_interface.o\
        mm_fill_jac.o\
        mm_fill_porous.o\
	mm_porous_EOS.o\
        mm_fill_potential.o\
        mm_fill_pthings.o\
        mm_fill_ptrs.o\
        mm_fill_terms.o\
        mm_fill_common.o\
        mm_fill_shell.o\
        mm_fill_solid.o\
        mm_fill_species.o\
        mm_fill_stress.o\
        mm_fill_rs.o\
        mm_fill_util.o\
        mm_flux.o\
        mm_input.o\
        mm_input_bc.o\
        mm_input_mp.o\
        mm_input_particles.o\
        mm_input_util.o\
        mm_interface.o\
        mm_matrl.o\
        mm_more_utils.o\
        mm_numjac.o\
        mm_ns_bc.o\
        mm_shell_bc.o\
        mm_placid.o\
        mm_post_proc.o\
        mm_post_proc_util.o\
        mm_prob_def.o\
        mm_propertyJac.o\
        mm_qp_storage.o\
        mm_shell_util.o\
        mm_sol_nonlinear.o\
        mm_species.o\
        mm_augc_util.o\
        mm_std_models.o\
        mm_qtensor_model.o\
        mm_unknown_map.o\
        mm_viscosity.o \
        mm_dil_viscosity.o

MM_INC=	mm_as.h\
        mm_as_alloc.h\
        mm_as_const.h\
        mm_as_structs.h\
        mm_augc_util.h\
        mm_bc.h\
        mm_cards.h\
        mm_chemkin.h\
        mm_eh.h\
        mm_elem_block.h\
        mm_elem_block_structs.h\
        mm_fill.h\
        mm_fill_aux.h\
        mm_fill_fill.h\
        mm_fill_jac.h\
        mm_fill_porous.h\
        mm_fill_potential.h\
        mm_fill_pthings.h\
        mm_fill_ptrs.h\
        mm_fill_rs.h\
        mm_fill_shell.h\
        mm_fill_solid.h\
        mm_fill_species.h\
        mm_fill_stress.h\
        mm_fill_terms.h\
        mm_fill_common.h\
        mm_fill_util.h\
        mm_flux.h\
        mm_input.h\
        mm_interface.h\
        mm_more_utils.h\
        mm_mp.h\
        mm_mp_const.h\
        mm_mp_structs.h\
        mm_names.h\
        mm_ns_bc.h\
        mm_shell_bc.h\
        mm_numjac.h\
        mm_parser.h\
        mm_post_def.h\
        mm_post_proc.h\
        mm_prob_def.h\
        mm_qp_storage.h\
        mm_qtensor_model.h\
        mm_shell_util.h\
        mm_sol_nonlinear.h\
        mm_species.h\
        mm_std_models.h\
        mm_unknown_map.h\
        mm_viscosity.h\
        mm_dil_viscosity.h


# _____ User routines "user_" prefix _________________________________________

U_SRC=	user_ac.c\
        user_bc.c\
        user_continuation.c\
        user_mp.c \
        user_mp_gen.c \
        user_out_hkm.c\
        user_post.c\
        user_pre.c\
        user_print.c \
        user_senkin.F

U_OBJ=	user_ac.o\
        user_bc.o\
        user_continuation.o\
        user_mp.o\
        user_mp_gen.o\
        user_out_hkm.o\
        user_post.o\
        user_pre.o\
        user_print.o \
        user_senkin.o

U_INC=  user_ac.h\
        user_bc.h\
        user_continuation.h\
        user_mp.h\
        user_mp_gen.h\
        user_out_hkm.h\
        user_post.h\
        user_pre.h\
        usr_print.h

LY_SRC=	mm_parser.tab.c\
        t.tab.c\
        lex.t.c\
        lex.yy.c

LY_OBJ=	mm_parser.tab.o\
        t.tab.o\
        lex.t.o\
        lex.yy.o

# To build a version that uses the lex/yacc parser, commnet out the following two lines
# and activate the -DNEW_PARSER line in the defines section:
LY_SRC=	
LY_OBJ=

G_SRC=  $(BC_SRC) $(EL_SRC) $(RF_SRC) $(MAIN_SRC) $(SL_SRC) $(MD_SRC) $(MM_SRC)\
        $(DP_SRC) $(AC_SRC) $(LOCA_SRC) $(GM_SRC) $(LY_SRC)

G_OBJ=  $(BC_OBJ) $(EL_OBJ) $(RF_OBJ) $(MAIN_OBJ) $(SL_OBJ) $(MD_OBJ) $(MM_OBJ)\
        $(DP_OBJ) $(AC_OBJ) $(LOCA_OBJ) $(GM_OBJ) $(LY_OBJ)

SRCS=	$(G_SRC) $(U_SRC)

OBJS=	$(G_OBJ) $(U_OBJ)

INCS=   $(EL_INC) $(MAIN_INC) $(RF_INC) $(SL_INC) $(MD_INC) $(DP_INC) $(BC_INC)\
        $(MM_INC) $(U_INC) $(AC_INC) $(LOCA_INC) $(GM_INC)


all: goma

#  Make all objects depend on all includes. This is reasonable for this project
#  with its interdepencies
$(OBJS) : $(INCS)

# Make the object files depend on this Makefile. We use static defines. So, any
# change in the Makefile where the defines are located should cause a global recompile.
$(OBJS): Goma_for_guts.mk

#
#  GOMA target with no CHEMKIN linking
#
goma:	$(G_LIB) $(U_LIB) $(LIB_DEP)
	$(RM) $@;\
	$(LD) $(LD_FLAGS) -t -o $@ $(G_LIB) $(U_LIB) $(LIBS) $(END_LIBS)
#
# GOMA target with CHEMKIN linked in
#
goma_c:	$(G_LIB) $(U_LIB) $(LIB_DEP) $(CHEMKIN_DEPS)
	$(RM) $@;\
	$(LD) $(LD_FLAGS) -o $@ $(G_LIB) $(U_LIB) $(LIBS) $(CHEMKIN_LIBS) \
               $(END_LIBS)
# To build a version with the lex/yacc parser, uncomment the lines below:
mm_parser.tab.c: mm_parser.y mm_table.l mm_table.y mm_table.l lex.yy.c t.tab.c lex.t.c mm_cards.h mm_parser.h
	bison -dvtl  mm_parser.y
t.tab.c: mm_table.y mm_table.l mm_parser.h
	bison -dvtl -pt -bt mm_table.y 
lex.yy.c: mm_parser.l mm_parser.h
	flex -l mm_parser.l
lex.t.c: mm_table.l mm_parser.h
	flex -ls -t mm_table.l | sed -f yy_to_t_lsed > lex.t.c
# This is the end of the block that needs to be uncommneted to build the lex/yacc parser.



gomad:	$(G_OBJ) $(U_OBJ) $(LIB_DEP)
	$(RM) $@;\
	$(LD) $(LD_FLAGS) -o $@ $(G_OBJ) $(U_OBJ) $(LIBS) $(END_LIBS)

debug:
	$(MAKE) CCFLAGS="-g" gomad


$(G_LIB): $(G_OBJ)
	$(AR) $(ARFLAGS) $@ $?;\
	$(RANLIB) $@ 

$(U_LIB): $(U_OBJ)
	$(AR) $(ARFLAGS) $@ $?;\
	$(RANLIB) $@ 


install: $(TARGET) $(G_LIB) $(U_LIB)
	@if [ ! -d $(BINDIR) ]; then\
		echo "Making $(BINDIR)";\
		mkdir -p $(BINDIR);\
	fi;\
	echo "Installing $(TARGET) in $(BINDIR)...";\
	$(INSTALL) $(INSTFLAGS) $(TARGET) $(BINDIR);\
	if [ ! -d $(LIBDIR) ]; then\
		echo "Making $(LIBDIR)";\
		mkdir -p $(LIBDIR);\
	fi;\
	echo "Installing $(G_LIB) $(U_LIB) in $(LIBDIR)...";\
	$(INSTALL) $(INSTFLAGS) $(G_LIB) $(LIBDIR);\
	$(INSTALL) $(INSTFLAGS) $(U_LIB) $(LIBDIR);\
	if [ ! -d $(INCDIR) ]; then\
		echo "Making $(INCDIR)";\
		mkdir -p $(INCDIR);\
	fi;\
	echo "Installing includes in $(INCDIR)...";\
	for i in $(INCS); do \
		echo "Installing $$i";\
		$(INSTALL) $(INSTFLAGS) $$i $(INCDIR);\
	done


clean:
	-rm -f *.o *~ *.BAK

realclean:
	-rm -f *.o *~ *.BAK goma libgoma.a libgomau.a

tags:
	find . -name "*.[ch]" -exec etags -a '{}' \;

depend:
	makedepend -- $(INCLUDES) $(CFLAGS) -- $(SRCS)

.c.o:
	$(CC) -c $(CFLAGS) $<

.F.o:
	$(FC) -c $(FFLAGS) $<

.C.o:
	$(CPP) -c $(CPPFLAGS) $< 

.fas.exoII:
	$(RM) $@
	$(FASTQ) -aprepro -mesh=$*.gen $*.fas
	$(EXO1EXO2) $*.gen $*.exoII
	$(RM) $*.gen

LINT = lint -Aa -Dlint -D__lint
FLINT = lintfor -Dlint $(COMPFLAGS) $(F77FLAGS) $(LDEFINES)

.c.ln:
	$(LINT) $(CFLAGS) $(CPPFLAGS) -c $<

%.ln:%.F
	$(FLINT) $(CPPFLAGS) -c $<

LN_GOMA = $(CP_LN) $(CK_LN) $(MD_LN) $(BC_LN)

goma_lint: $(LN_GOMA)
	$(LINT) $(LN_GOMA) -lm

.SUFFIXES: .exoII .gen .fas .ln


# DO NOT DELETE THIS LINE -- make depend depends on it.
