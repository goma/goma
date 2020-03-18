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

# Main goma makefile, specific user settings should go in
# the USER_SETTINGS_FILE

MAKEFILE_NAME = Makefile

# Goma relative source directories
BIN_DIR = ../bin
OBJ_DIR = ../build
LIB_DIR = ../lib
INCLUDE_DIR = ../include
BRKFIX_DIR = brkfix
USER_SETTINGS_FILE = ../settings.mk

# _MACH is used if libraries were compiled with machine specific suffixes
# (libblas_linux)

_MACH =

# Useful for standard goma libary install location

ARCH = rhel6_ompi 

# Goma library path, by default we assume all major non-standard libraries are under
# the goma library path

GOMA_LIBS = /home/goma/production/$(ARCH)

# MPI install location, used for include, lib, and bin directories

MPI_TOP = /usr/bin/openmpi-2.1.1

# GCC_TOP is the corresponding location of the gcc compilers and libs that are
# used by the MPI compilers
GCC_TOP = /usr/lib64

# PURIFY
# ------
# Examples:
# 	PURE = purify
#	PURE = purify -cache-dir="/tmp"
#	PURE = quantify -cache-dir="/tmp"
PURE =

# MPI compilers

CC = $(MPI_TOP)/bin/mpicc
CXX = $(MPI_TOP)/bin/mpicxx
FC = $(MPI_TOP)/bin/mpif77

# USER OVERRIDES from settings.mk
-include $(USER_SETTINGS_FILE)

# Applications and settings
# =========================

LD = $(PURE) $(CXX)

SHELL = /bin/sh

# install
INSTALL ?= cp
INSTFLAGS ?= 

# clean
RM ?= rm -f

# archiver
AR ?= ar
ARFLAGS ?= srv
RANLIB ?= touch

# LOCAL DEFINES SECTION
# =====================
#
#  USE_CHEMKIN: If you want to use chemkin, you need to define this
#               compiler definition during the compile.
#  SENKIN_OUTPUT: If you want to use senkin postprocessor, define
#                 this along with the USE_CHEMKIN compiler directive.
#  PARALLEL:    If you want to run parallel applications, define
#               This. This works even for serial applications, if
#               you add the prefix "mpirun -np 1" to your command
#               line. Therefore, I have made this the default.
#  MATRIX_DUMP: Define this if you want to enable the dumping of
#               matrices to disk.
#  TRILINOS:    Define this if you want to link to the new Trilinos
#               linear solver package, including its version of Aztec.
#
#  'platform:'  Several 'ifdef's on platform type still exist.
#               The following platforms must contain a -D'platform'
#               DEFINE on the compile line:
#                      Linux       -Dlinux
#                      (if 64bit compiler)     -DCOMPILER_64BIT

DEFINES ?= -Dlinux \
           -DCOMPILER_64BIT \
           -DENABLE_AMESOS \
           -DTRILINOS \
           -DCHECK_FINITE\
           -DNO_CHEBYSHEV_PLEASE \
           -DMDE=27 \
           -DMAX_PROB_VAR=15 \
           -DMAX_EXTERNAL_FIELD=4 \
           -DMAX_CONC=4 \
           -DCOUPLED_FILL \
           -DPARALLEL \
           -DHAVE_STRATIMIKOS \
           -DHAVE_SPARSE \
           -DEIGEN_SERIAL

# Other defines:
#          -DCHECKFINITE
#          -DMAX_CONC=1 \
#          -DEIGEN_SERIAL \
#          -DEIGEN_PARALLEL \
#          -DNO_CHEBYSHEV_PLEASE \
#          -DNEW_PARSER \
#          -DSUN_SUNOS_57 \
#          -D_RWSTD_COMPILE_INSTANTIATE \
#          -DCOUPLED_FILL \
#          -DXCODE \
#          -DLIBRARY_MODE \
#          -DUSE_CHEMKIN -DSENKIN_OUTPUT -DDEBUG_HKM \

# Git Version information
# check for executable
HAS_GIT := $(shell command -v git 2> /dev/null)



ifdef HAS_GIT
	GIT_VERSION=$(shell git describe --dirty --always --tags 2> /dev/null)

ifneq "$(GIT_VERSION)" ""
	VERSION = -DGIT_VERSION=\"$(GIT_VERSION)\"
endif

endif

# External Libraries
# ==================

# MPI
# -----------
# MPI includes
GOMA_MPI_INC ?= -I $(MPI_TOP)/include

# MPI libraries
GOMA_MPI_LIB ?= -L$(MPI_TOP)/lib  -lmpi  -lmpi_mpifh

# For MPI, the libraries need to be explicitly set in the path:
GCC_LIB = $(GCC_TOP)/lib
GCC_LIB64 = $(GCC_TOP)/lib64
LD_LIBRARY_PATH := $(GCC_LIB64):$(GCC_LIB):${LD_LIBRARY_PATH}

# TRILINOS
# --------
# Trilinos install directory (include, lib)
# If built with build-goma-dependencies.sh, includes flags for SEACAS, UMFpack, and all other necessary
# libraries and includes otherwise not listed below.

TRILINOS_TOP ?= $(GOMA_LIBS)/trilinos-12.12.1-Built

include $(TRILINOS_TOP)/include/Makefile.export.Trilinos

# Trilinos includes and libraries
TRILINOS_INC_REMOVE_USR_INCLUDE = $(subst -I/usr/include, ,  $(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS))
TRILINOS_INC ?= $(subst -I, -isystem, $(TRILINOS_INC_REMOVE_USR_INCLUDE))
TRILINOS_LIB ?= $(Trilinos_LIBRARY_DIRS) $(Trilinos_TPL_LIBRARY_DIRS) \
$(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES)

# ARPACK
# ------
#
# To use ARPACK, either EIGEN_SERIAL or EIGEN_PARALLEL must be defined above.
# Default is to link to ARPACK

ARPACK_TOP ?= $(GOMA_LIBS)/ARPACK

# EIGEN_SERIAL: arpack only

ARPACK_LIB ?= -L$(ARPACK_TOP)/lib -larpack$(_MACH)

# FRONTAL SOLVER LIBRARY
# ----------------------
#
# The default is not to include the frontal solver in compile. 
# It no longer functions properly under 64-bit builds and should not 
# be utilized as a solver.

#FRONT_LIB ?= -L$(GOMA_LIBS)/front/lib  -lmpfront$(_MACH) \
              -lparlib$(_MACH)


# SPARSE LIBRARY
# --------------

SPARSE_TOP ?= $(GOMA_LIBS)/sparse

SPARSE_LIB ?= -L$(SPARSE_TOP)/lib  -lsparse$(_MACH)

# Sparse includes are located in the goma directory
SPARSE_INC ?= -I$(SPARSE_TOP)/src

# FORTRAN LIBRARIES
# -----------------

FORTRAN_LIB ?= -lc  -lgcc  -lstdc++ -lgfortran

# SYSTEM LIBRARIES
# ----------------

SYS_LIB ?= -lm  -lz

# INCLUDES
# --------

GOMA_INC = -I$(INCLUDE_DIR)

INCLUDES = $(GOMA_INC) $(SPARSE_INC) $(TRILINOS_INC) \
           $(GOMA_MPI_INC) $(SEACAS_INC) $(UMFPACK_INC)



# Compiler settings
# =================

# override CCFLAGS with DEBUG_CCFLAGS if building debug target
DEBUG_FLAGS ?= -g -gdwarf-2

# All compiler options
CCFLAGS ?= -O2 -DNDEBUG

# C/C++ warning flags
C_WARN_FLAGS ?= -Wall -Wextra -Wno-unused-parameter
CXX_WARN_FLAGS ?= -Wall -Wextra

# Fortran warning flags
F_WARN_FLAGS ?=

# Fortran flags
FFLAGS = $(CCFLAGS) $(F_WARN_FLAGS) $(DEFINES) $(INCLUDES)

# C compiler flags
CSTD ?= -std=gnu99
CFLAGS = $(CCFLAGS) $(C_WARN_FLAGS) $(DEFINES) $(VERSION) $(INCLUDES) $(CSTD)

# C++ compiler flags
CXXSTD ?= -std=c++11
CXXFLAGS = $(CCFLAGS) $(CXX_WARN_FLAGS) $(DEFINES) $(INCLUDES) -Wno-deprecated $(CXXSTD)

# linker flags
LD_EXTRA_FLAGS ?=
LD_FLAGS = $(CCFLAGS) $(Trilinos_EXTRA_LD_FLAGS) $(LD_EXTRA_FLAGS)

# LIBRARIES
# ---------
LIBS = $(FRONT_LIB) \
       $(SPARSE_LIB) \
       $(TRILINOS_LIB) \
       $(Y12M_LIB) \
       $(SUPERLU_LIB) \
       $(MUMPS_LIB) \
       $(PARMETIS_LIB) \
       $(LAPACK_LIB) \
       $(ARPACK_LIB) \
       $(HARWELL_LIB) \
       $(UMFPACK_LIB) \
       $(BLAS_LIB) \
       $(SCALAPACK_LIB) \
       $(SEACAS_LIB) 

END_LIBS = $(GOMA_MPI_LIB) \
           $(FORTRAN_LIB) \
           $(SYS_LIB)


LIB_DEP = $(FRONT_DEP)

# Source files for goma
# =====================

# _____ Brk and Fix ___________________________________________________________

BFU_SRC = brk_utils.c 

BFU_INC = brk_utils.h

BF_SRC = brk_exo_file.c fix_exo_file.c ppi.c wr_coords.c rd_in.c \
	 mk_dm.c sam_perea.c utils.c wr_graph_file.c emuck.c bbb.c

BF_SRC := ${addprefix $(BRKFIX_DIR)/, $(BF_SRC)}

BF_INC = 

# _____ Element routines "el_" prefix _________________________________________

EL_SRC=	el_elm_info.c el_quality.c exo_conn.c

EL_INC=	el_elm.h\
        el_elm_info.h\
        el_geom.h\
        el_quality.h\
        exo_conn.h\
        exo_struct.h


# _____ Main routines _________________________________________________________

MAIN_SRC= main.c

MAIN_INC= std.h

# _____ Advanced capabilities routines "ac_" prefix ___________________________

AC_SRC= ac_conti.c\
        ac_hunt.c\
        ac_update_parameter.c\
        ac_display_parameter.c\
        ac_loca_interface.c\
        ac_stability.c\
        ac_stability_util.c\
        ac_particles.c

AC_INC= ac_conti.h\
        ac_hunt.h\
        ac_update_parameter.h\
        ac_stability.h\
        ac_stability_util.h\
        ac_particles.h


# _____ Library of continuation algorithms "loca_" prefix _____________________

LOCA_SRC= loca_lib.c\
          loca_bord.c\
          loca_util.c\
          loca_eigenvalue.c \
          loca_eigen_c2f.F \
          loca_eigen_cayley.F

LOCA_INC= loca_const.h\
          loca_util_const.h


# _____	Reacting flow routines "rf_" prefix ___________________________________

RF_SRC= rf_allo.c\
        rd_dpi.c\
        rf_element_storage.c\
        rd_exo.c\
        rd_mesh.c\
	rd_pixel_image.c \
	rd_pixel_image2.c \
        rf_node.c\
        rf_node_vars.c\
        rf_pre_proc.c\
        rf_setup_problem.c\
        rf_shape.c\
        rf_solve.c\
	rf_solve_segregated.c\
        rf_util.c\
        rf_vars.c\
        wr_dpi.c\
        wr_side_data.c\
        wr_soln.c\
        wr_exo.c

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
        rf_io_structs.h\
        rf_masks.h\
        rf_mp.h\
        rf_node.h\
        rf_node_const.h\
        rf_pre_proc.h\
        rf_shape.h\
        rf_solve.h\
	rf_solve_segregated.h\
        rf_solver.h\
        rf_solver_const.h\
        rf_util.h\
        rf_vars_const.h \
        rd_dpi.h\
        rd_exo.h\
        rd_mesh.h\
	rd_pixel_image.h\
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
        sl_amesos_interface.C\
        sl_leastsquares.C \
        sl_epetra_interface.cpp \
        sl_epetra_util.cpp \
        sl_aztecoo_interface.cpp \
        sl_stratimikos_interface.cpp

SL_INC= sl_aux.h\
        sl_auxutil.h\
        sl_eggroll.h\
        sl_eggroll_def.h\
        sl_lu.h\
        sl_matrix_util.h\
        sl_umf.h\
        sl_util.h\
        sl_util_structs.h \
        sl_rcm.h \
        sl_epetra_interface.h \
        sl_epetra_util.h \
        sl_aztecoo_interface.h \
        sl_stratimikos_interface.h


# _____	Machine dependent routines "md_" prefix _______________________________

MD_SRC=	md_timer.c \
        md_ieee.c

MD_INC= md_timer.h


# _____ Distributed processing routines "dp_" prefix __________________________

DP_SRC= dp_comm.c\
        dp_map_comm_vec.c\
        dp_utils.c\
        dp_vif.c

DP_INC= dp_comm.h\
        dp_map_comm_vec.h\
        dp_types.h\
        dp_utils.h\
        dp_vif.h\
        dpi.h


# _____	Boundary condition routines "bc_" prefix ______________________________

BC_SRC= bc_colloc.c\
        bc_curve.c\
        bc_dirich.c\
        bc_integ.c\
        bc_rotate.c\
        bc_special.c\
        bc_contact.c \
        bc_surfacedomain.c \
	rotate_util.c

BC_INC= bc_colloc.h\
        bc_curve.h\
        bc_dirich.h\
        bc_integ.h\
        bc_rotate.h\
        bc_special.h\
        bc_contact.h \
        bc_surfacedomain.h \
	rotate_util.h


# _____ Moving mesh routines "mm_" prefix _____________________________________

MM_SRC= mm_as_alloc.c\
        mm_augc_util.c\
        mm_bc.c\
        mm_bc_conflict.c\
        mm_chemkin.c\
        mm_dil_viscosity.c\
        mm_eh.c\
        mm_fill.c\
        mm_fill_aux.c\
        mm_fill_fill.c\
	mm_fill_ls.c \
        mm_fill_interface.c\
        mm_fill_jac.c\
        mm_fill_porous.c\
        mm_fill_population.c\
	mm_porous_EOS.c\
        mm_fill_potential.c\
        mm_fill_pthings.c\
        mm_fill_ptrs.c\
        mm_fill_rs.c\
        mm_fill_shell.c\
        mm_fill_solid.c\
        mm_fill_species.c\
        mm_fill_stress.c\
        mm_fill_common.c\
        mm_fill_terms.c\
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
        mm_std_models_shell.c\
        mm_unknown_map.c\
        mm_viscosity.c

MM_INC=	mm_as.h\
        mm_as_alloc.h\
        mm_as_const.h\
        mm_as_structs.h\
        mm_augc_util.h\
        mm_bc.h\
        mm_chemkin.h\
        mm_dil_viscosity.h\
        mm_eh.h\
        mm_elem_block.h\
        mm_elem_block_structs.h\
        mm_fill.h\
        mm_fill_aux.h\
        mm_fill_fill.h\
        mm_fill_jac.h\
        mm_fill_porous.h\
        mm_fill_population.h\
        mm_fill_potential.h\
        mm_fill_pthings.h\
        mm_fill_ptrs.h\
        mm_fill_rs.h\
        mm_fill_shell.h\
        mm_fill_solid.h\
        mm_fill_species.h\
        mm_fill_stress.h\
        mm_fill_common.h\
        mm_fill_terms.h\
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
        mm_post_def.h\
        mm_post_proc.h\
        mm_prob_def.h\
        mm_qp_storage.h\
        mm_qtensor_model.h\
        mm_shell_util.h\
        mm_sol_nonlinear.h\
        mm_species.h\
        mm_std_models.h\
	mm_std_models_shell.h\
        mm_unknown_map.h\
        mm_viscosity.h

SHELL_SRC= shell_tfmp_util.c

SHELL_INC= shell_tfmp_util.h

# _____ User routines "user_" prefix __________________________________________

USER_SRC= user_ac.c\
          user_bc.c\
          user_continuation.c\
          user_mp.c \
          user_mp_gen.c \
          user_out_hkm.c\
          user_post.c\
          user_pre.c\
          user_print.c \
          user_senkin.F

USER_INC= user_ac.h\
          user_bc.h\
          user_continuation.h\
          user_mp.h\
          user_mp_gen.h\
          user_out_hkm.h\
          user_post.h\
          user_pre.h\
          usr_print.h

# Goma Data Structures library gds

GDS_SRC=gds/gds_vector.c
GDS_INC=gds/gds_vector.h

# _____ Full macros ___________________________________________________________

GOMA_SRC = $(BC_SRC) $(EL_SRC) $(RF_SRC) $(MAIN_SRC) $(SL_SRC) $(MD_SRC) $(MM_SRC)\
           $(SHELL_SRC) $(DP_SRC) $(AC_SRC) $(LOCA_SRC) $(BFU_SRC) 

BRKFIX_SRC = $(BF_SRC)

GOMA_INCS = $(EL_INC) $(MAIN_INC) $(RF_INC) $(SL_INC) $(MD_INC) $(DP_INC) $(BC_INC)\
           $(MM_INC) $(SHELL_INC) $(USER_INC) $(AC_INC) $(LOCA_INC) $(BF_INC) 

GOMA_INCS := ${addprefix $(INCLUDE_DIR)/, $(GOMA_INCS)}

GOMA_OBJ := ${GOMA_SRC:.c=.o}
GOMA_OBJ := ${GOMA_OBJ:.F=.o}
GOMA_OBJ := ${GOMA_OBJ:.C=.o}
GOMA_OBJ := ${GOMA_OBJ:.cpp=.o}

USER_OBJ := ${USER_SRC:.c=.o}
USER_OBJ := ${USER_OBJ:.F=.o}
USER_OBJ := ${USER_OBJ:.C=.o}

BRKFIX_OBJ :=  ${BRKFIX_SRC:.c=.o}
BRKFIX_OBJ :=  ${BRKFIX_OBJ:.F=.o}
BRKFIX_OBJ :=  ${BRKFIX_OBJ:.C=.o}

FIX_OBJ :=  ${FIX_SRC:.c=.o}
FIX_OBJ :=  ${FIX_OBJ:.F=.o}
FIX_OBJ :=  ${FIX_OBJ:.C=.o}

GDS_OBJ := ${GDS_SRC:.c=.o}

RELEASE_OBJ_DIR = $(OBJ_DIR)/release
DEBUG_OBJ_DIR = $(OBJ_DIR)/debug

GOMA_RELEASE_OBJS := ${addprefix $(RELEASE_OBJ_DIR)/,$(GOMA_OBJ)}
BRKFIX_RELEASE_OBJS := ${addprefix $(RELEASE_OBJ_DIR)/,$(BRKFIX_OBJ)}
USER_RELEASE_OBJS := ${addprefix $(RELEASE_OBJ_DIR)/,$(USER_OBJ)}
GDS_RELEASE_OBJS := ${addprefix $(RELEASE_OBJ_DIR)/,$(GDS_OBJ)}

GOMA_DEBUG_OBJS := ${addprefix $(DEBUG_OBJ_DIR)/,$(GOMA_OBJ)}
USER_DEBUG_OBJS := ${addprefix $(DEBUG_OBJ_DIR)/,$(USER_OBJ)}
BRKFIX_DEBUG_OBJS := ${addprefix $(DEBUG_OBJ_DIR)/,$(BRKFIX_OBJ)}
GDS_DEBUG_OBJS := ${addprefix $(DEBUG_OBJ_DIR)/,$(GDS_OBJ)}

# combined sources
SRCS  = $(GOMA_SRC) $(USER_SRC)

# combined objects for dependency compilation
RELEASE_OBJS = $(GOMA_RELEASE_OBJS) $(USER_RELEASE_OBJS)
DEBUG_OBJS = $(GOMA_DEBUG_OBJS) $(USER_DEBUG_OBJS)

# File names used for binaries and libraries
# ==========================================
GOMA_LIB = libgoma.a
GOMA_USER_LIB = libgomau.a
GDS_LIB = libgds.a
TARGET = goma
DEBUG_TARGET = gomad

# Install location
# ================
# will install files in $PREFIX/include/goma, $PREFIX/lib, $PREFIX/bin

PREFIX ?= /usr/local

PREFIX_BIN = $(PREFIX)/bin
PREFIX_LIB = $(PREFIX)/lib
PREFIX_INC = $(PREFIX)/include

# MAKE TARGETS
# ============

# Object files depend on includes and the makefile
RELEASE_DEPENDS = $(GOMA_INCS) $(MAKEFILE_NAME) $(USER_SETTINGS_FILE)
DEBUG_DEPENDS = $(GOMA_INCS) $(MAKEFILE_NAME) $(USER_SETTINGS_FILE)

# Main targets
# ------------
.PHONY: all

all: $(BIN_DIR)/$(TARGET)

# bin/goma executable
$(BIN_DIR)/$(TARGET): $(LIB_DIR)/$(GOMA_LIB) $(LIB_DIR)/$(GOMA_USER_LIB) $(LIB_DIR)/$(GDS_LIB) $(LIB_DEP)
	mkdir -p $(BIN_DIR); \
	$(RM) $@; \
	$(LD) $(LD_FLAGS) -o $@ $(LIB_DIR)/$(GOMA_LIB) $(LIB_DIR)/$(GOMA_USER_LIB) $(LIB_DIR)/$(GDS_LIB) $(LIBS) $(END_LIBS) $(LD_FLAGS)

# goma library
$(LIB_DIR)/$(GOMA_LIB): $(GOMA_RELEASE_OBJS) $(BRKFIX_RELEASE_OBJS)
	mkdir -p $(LIB_DIR); \
	$(AR) $(ARFLAGS) $@ $?; \
	$(RANLIB) $@ 

# goma user library
$(LIB_DIR)/$(GOMA_USER_LIB): $(USER_RELEASE_OBJS)
	@mkdir -p $(LIB_DIR); \
	$(AR) $(ARFLAGS) $@ $?; \
	$(RANLIB) $@ 

# goma data structures library
$(LIB_DIR)/$(GDS_LIB): $(GDS_RELEASE_OBJS)
	@mkdir -p $(LIB_DIR); \
	$(AR) $(ARFLAGS) $@ $?; \
	$(RANLIB) $@ 

# Debug targets
# -------------

.PHONY: debug

debug:
	$(MAKE) CCFLAGS="$(DEBUG_FLAGS)" $(BIN_DIR)/$(DEBUG_TARGET)

# bin/gomad executable
$(BIN_DIR)/$(DEBUG_TARGET): $(GOMA_DEBUG_OBJS) $(BRKFIX_DEBUG_OBJS) $(USER_DEBUG_OBJS) $(GDS_DEBUG_OBJS) $(LIB_DEP)
	mkdir -p $(BIN_DIR); \
	$(RM) $@; \
	$(LD) $(LD_FLAGS) -o $@ $(GOMA_DEBUG_OBJS) $(BRKFIX_DEBUG_OBJS) $(USER_DEBUG_OBJS) $(GDS_DEBUG_OBJS) $(LIBS) $(END_LIBS) $(LD_FLAGS)

# Object pattern rules
# --------------------

# Objects depend on c, F, and C files of the same name

# release objects

$(RELEASE_OBJ_DIR)/%.o: %.c $(RELEASE_DEPENDS)
	@mkdir -p $(RELEASE_OBJ_DIR)/$(BRKFIX_DIR)
	@mkdir -p $(RELEASE_OBJ_DIR)/gds
	$(CC) -c $(CFLAGS) -o $@ $<

$(RELEASE_OBJ_DIR)/%.o: %.F $(RELEASE_DEPENDS)
	@mkdir -p $(RELEASE_OBJ_DIR)
	$(FC) -c $(FFLAGS) -o $@ $<

$(RELEASE_OBJ_DIR)/%.o: %.C $(RELEASE_DEPENDS)
	@mkdir -p $(RELEASE_OBJ_DIR)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

$(RELEASE_OBJ_DIR)/%.o: %.cpp $(RELEASE_DEPENDS)
	@mkdir -p $(RELEASE_OBJ_DIR)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

# debug objects
$(DEBUG_OBJ_DIR)/%.o: %.c $(DEBUG_DEPENDS)
	@mkdir -p $(DEBUG_OBJ_DIR)/$(BRKFIX_DIR)
	@mkdir -p $(DEBUG_OBJ_DIR)/gds
	$(CC) -c $(CFLAGS) -o $@ $<

$(DEBUG_OBJ_DIR)/%.o: %.F $(DEBUG_DEPENDS)
	@mkdir -p $(DEBUG_OBJ_DIR)
	$(FC) -c $(FFLAGS) -o $@ $<

$(DEBUG_OBJ_DIR)/%.o: %.C $(DEBUG_DEPENDS)
	@mkdir -p $(DEBUG_OBJ_DIR)
	$(CXX) -c $(CXXFLAGS) -o $@ $<

$(DEBUG_OBJ_DIR)/%.o: %.cpp $(DEBUG_DEPENDS)
	@mkdir -p $(RELEASE_OBJ_DIR)
	$(CXX) -c $(CXXFLAGS) -o $@ $<



# Installation target
# -------------------

.PHONY: install

install: $(BIN_DIR)/$(TARGET) $(LIB_DIR)/$(GOMA_LIB) $(LIB_DIR)/$(GOMA_USER_LIB)
	@if [ ! -d $(PREFIX_BIN) ]; then\
		echo "Making $(PREFIX_BIN)";\
		mkdir -p $(PREFIX_BIN);\
	fi;\
	echo "Installing $(BIN_DIR)/$(TARGET) in $(PREFIX_BIN)...";\
	$(INSTALL) $(INSTFLAGS) $(BIN_DIR)/$(TARGET) $(PREFIX_BIN);\
	if [ ! -d $(PREFIX_LIB) ]; then\
		echo "Making $(PREFIX_LIB)";\
		mkdir -p $(PREFIX_LIB);\
	fi;\
	echo "Installing $(GOMA_LIB) $(GOMA_USER_LIB) in $(PREFIX_LIB)...";\
	$(INSTALL) $(INSTFLAGS) $(LIB_DIR)/$(GOMA_LIB) $(PREFIX_LIB);\
	$(INSTALL) $(INSTFLAGS) $(LIB_DIR)/$(GOMA_USER_LIB) $(PREFIX_LIB);\
	if [ ! -d $(PREFIX_INC) ]; then\
		echo "Making $(PREFIX_INC)";\
		mkdir -p $(PREFIX_INC);\
	fi;\
	echo "Installing includes in $(PREFIX_INC)...";\
	for i in $(GOMA_INCS); do \
		echo "Installing $$i";\
		$(INSTALL) $(INSTFLAGS) $$i $(PREFIX_INC);\
	done


# Cleanup targets
# --------------

.PHONY: clean

clean:
	$(RM) $(DEBUG_OBJ_DIR)/*.o $(RELEASE_OBJ_DIR)/*.o  $(LIB_DIR)/*.a \
	$(DEBUG_OBJ_DIR)/$(BRKFIX_DIR)/*.o $(RELEASE_OBJ_DIR)/$(BRKFIX_DIR)/*.o

.PHONY: realclean

realclean: clean
	$(RM) $(BIN_DIR)/$(TARGET) $(BIN_DIR)/$(DEBUG_TARGET);
