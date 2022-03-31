cmake_minimum_required(VERSION 3.17)
cmake_policy(VERSION 3.17)
project(goma)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${goma_SOURCE_DIR}/cmake")

set(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})

# Get Trilinos as one entity
message(STATUS "Trilinos_DIR = ${Trilinos_DIR}")
find_package(Trilinos REQUIRED)
# message(STATUS "Trilinos LD Flags = ${Trilinos_EXTRA_LD_FLAGS}")

set(MPI_C_COMPILER ${Trilinos_C_COMPILER})
set(MPI_CXX_COMPILER ${Trilinos_CXX_COMPILER})
set(MPI_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
message(STATUS "MPI_C_COMPILER: ${MPI_C_COMPILER}")
message(STATUS "MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")
message(STATUS "MPI_Fortran_COMPILER: ${MPI_Fortran_COMPILER}")
find_package(MPI REQUIRED)
# Avoid linking to too many libraries by removing duplicates reverse so that
# correct order for linking is preserved
list(REVERSE Trilinos_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
list(REVERSE Trilinos_LIBRARIES)
#
list(REVERSE Trilinos_TPL_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list(REVERSE Trilinos_TPL_LIBRARIES)

if(NOT GOMA_COMPILER_OVERRIDE)
  set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
  set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
  set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
else()
  if(NOT (${CMAKE_C_COMPILER} STREQUAL ${Trilinos_C_COMPILER}))
    message(
      WARNING
        "CMAKE_C_COMPILER ${CMAKE_C_COMPILER} different from Trilinos_C_COMPILER ${Trilinos_C_COMPILER}"
    )
    message(WARNING "force setting compiler with GOMA_COMPILER_OVERRIDE=OFF")
  endif()
  if(NOT (${CMAKE_CXX_COMPILER} STREQUAL ${Trilinos_CXX_COMPILER}))
    message(
      WARNING
        "CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER} different from Trilinos_CXX_COMPILER ${Trilinos_CXX_COMPILER}"
    )
    message(WARNING "force setting compiler with GOMA_COMPILER_OVERRIDE=OFF")
  endif()
  if(NOT (${CMAKE_Fortran_COMPILER} STREQUAL ${Trilinos_Fortran_COMPILER}))
    message(
      WARNING
        "CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER} different from Trilinos_Fortran_COMPILER ${Trilinos_Fortran_COMPILER}"
    )
    message(WARNING "force setting compiler with GOMA_COMPILER_OVERRIDE=OFF")
  endif()
endif()

option(ENABLE_CCACHE OFF)
if(ENABLE_CCACHE)
  message(STATUS "Searching for ccache.")
  find_program(CCACHE_FOUND ccache)
  if(CCACHE_FOUND)
    message(STATUS "Ccache found.")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
  endif(CCACHE_FOUND)
endif()

set(GOMA_COMPILER_OVERRIDE
    OFF
    CACHE
      BOOL
      "allow goma to compile with different compilers than TPLs, use with caution"
)

# Trilinos required packages
foreach(package Amesos AztecOO Epetra)
  list(FIND Trilinos_PACKAGE_LIST ${package} trilinos_package_index)
  if(${trilinos_package_index} LESS 0)
    message(FATAL_ERROR "Trilinos not built with ${package} a required package")
  endif()
endforeach()

# Trilinos required TPLs
foreach(package MPI)
  list(FIND Trilinos_TPL_LIST ${package} trilinos_package_index)
  if(${trilinos_package_index} LESS 0)
    message(FATAL_ERROR "Trilinos not built with TPL ${package} a required TPL")
  endif()
endforeach()

# Trilinos extra TPLs
foreach(package SuperLUDist MUMPS)
  list(FIND Trilinos_TPL_LIST ${package} trilinos_package_index)
  if(${trilinos_package_index} LESS 0)
    message(
      WARNING
        "Trilinos not built with TPL ${package} some features of Goma won't work"
    )
  endif()
endforeach()

# stratimikos and teko
list(FIND Trilinos_PACKAGE_LIST Stratimikos stratimikos_package_index)
if(${stratimikos_package_index} GREATER_EQUAL 0)
  option(ENABLE_STRATIMIKOS "ENABLE_STRATIMIKOS" ON)
  if(ENABLE_STRATIMIKOS)
    message(STATUS "TRILINOS: Stratimikos found, enabling in Goma")
    list(APPEND GOMA_COMPILE_DEFINITIONS GOMA_ENABLE_STRATIMIKOS)
    list(FIND Trilinos_PACKAGE_LIST Teko teko_package_index)
    if(${teko_package_index} GREATER_EQUAL 0)
      option(ENABLE_TEKO "ENABLE_TEKO" ON)
      if(ENABLE_TEKO)
        message(STATUS "TRILINOS: Teko found, enabling in Goma")
        list(APPEND GOMA_COMPILE_DEFINITIONS GOMA_ENABLE_TEKO)
      endif()
    endif()
  endif()
endif()

# Exodus
find_package(SEACASExodus REQUIRED HINTS ${Trilinos_DIR}/../SEACASExodus)
message(STATUS "SEACASExodus_DIR = ${SEACASExodus_DIR}")

option(ENABLE_APREPRO_LIB "ENABLE_APREPRO_LIB" ON)
if(ENABLE_APREPRO_LIB)
  find_package(SEACASAprepro HINTS ${SEACASExodus_DIR}/../SEACASAprepro)
  if(SEACASAprepro_FOUND)
    message(STATUS "Aprepro Library found, enabling")
    message(STATUS "SEACASAprepro_DIR = ${SEACASAprepro_DIR}")
  endif()
endif()

message(STATUS "PETSC_DIR = $ENV{PETSC_DIR}")
option(ENABLE_PETSC "ENABLE_PETSC" ON)
if(ENABLE_PETSC)
  # PETSC search set root of location to find PETSc configuration
  set(PETSC_PATH $ENV{PETSC_DIR}/$ENV{PETSC_ARCH})
  set(ENV{PKG_CONFIG_PATH} ${PETSC_PATH}/lib/pkgconfig)

  find_package(PkgConfig)

  if(PkgConfig_FOUND)
    pkg_search_module(PETSC IMPORTED_TARGET PETSc)

    # Determine PETSc compilers
    if(PETSC_FOUND)
      execute_process(
        COMMAND ${PKG_CONFIG_EXECUTABLE} PETSc --variable=prefix
        COMMAND tr -d '\n'
        OUTPUT_VARIABLE PETSC_PREFIX)
      message(STATUS "Using PETSc prefix ${PETSC_PREFIX}")
      execute_process(
        COMMAND ${PKG_CONFIG_EXECUTABLE} PETSc --variable=ccompiler
        COMMAND tr -d '\n'
        OUTPUT_VARIABLE PETSC_C_COMPILER)
      execute_process(
        COMMAND ${PKG_CONFIG_EXECUTABLE} PETSc --variable=cxxcompiler
        COMMAND tr -d '\n'
        OUTPUT_VARIABLE PETSC_CXX_COMPILER)
      execute_process(
        COMMAND ${PKG_CONFIG_EXECUTABLE} PETSc --variable=fcompiler
        COMMAND tr -d '\n'
        OUTPUT_VARIABLE PETSC_FORTRAN_COMPILER)
      if(ENABLE_PETSC)
        if(NOT (Trilinos_C_COMPILER STREQUAL PETSC_C_COMPILER))
          if(GOMA_COMPILER_OVERRIDE)
            message(
              WARNING "Trilinos and Petsc compiled with different C compilers\n"
                      "Trilinos C: ${Trilinos_C_COMPILER}\n"
                      "PETSc C: ${PETSC_C_COMPILER}")
          else()
            message(
              FATAL_ERROR
                "Trilinos and Petsc compiled with different C compilers\n"
                "Trilinos C: ${Trilinos_C_COMPILER}\n"
                "PETSc C: ${PETSC_C_COMPILER}"
                "Disable error by GOMA_COMPILER_OVERRIDE=ON if you know what you are doing"
            )
          endif()
        endif()
        if(NOT (Trilinos_Fortran_COMPILER STREQUAL PETSC_FORTRAN_COMPILER))
          if(GOMA_COMPILER_OVERRIDE)
            message(
              WARNING
                "Trilinos and Petsc compiled with different Fortran compilers\n"
                "Trilinos Fortran: ${Trilinos_Fortran_COMPILER}\n"
                "PETSc Fortran: ${PETSC_FORTRAN_COMPILER}")
          else()
            message(
              FATAL_ERROR
                "Trilinos and Petsc compiled with different Fortran compilers\n"
                "Trilinos Fortran: ${Trilinos_Fortran_COMPILER}\n"
                "PETSc Fortran: ${PETSC_FORTRAN_COMPILER}"
                "Disable error by GOMA_COMPILER_OVERRIDE=ON if you know what you are doing"
            )
          endif()
        endif()
        if(NOT (Trilinos_CXX_COMPILER STREQUAL PETSC_CXX_COMPILER))
          if(GOMA_COMPILER_OVERRIDE)
            message(
              WARNING
                "Trilinos and Petsc compiled with different CXX compilers\n"
                "Trilinos CXX: ${Trilinos_CXX_COMPILER}\n"
                "PETSc CXX: ${PETSC_CXX_COMPILER}")
          else()
            message(
              FATAL_ERROR
                "Trilinos and Petsc compiled with different CXX compilers\n"
                "Trilinos CXX: ${Trilinos_CXX_COMPILER}\n"
                "PETSc CXX: ${PETSC_CXX_COMPILER}"
                "Disable error by GOMA_COMPILER_OVERRIDE=ON if you know what you are doing"
            )
          endif()
        endif()
      endif()
    endif()

  else()
    message(WARNING "pkg-config not found, not searching for PETSc")
  endif()
else()
  message(STATUS "PETSc is disabled")
endif()

set(GOMA_TPL_LIBRARIES "")
set(GOMA_TPL_LIBRARIES "")

option(ENABLE_ARPACK "ENABLE_ARPACK" ON)
if(ENABLE_ARPACK)
  set(ARPACK_PREFIX
      "${Trilinos_DIR}/../../../"
      CACHE PATH "set ARPACK library location")

  find_package(ARPACK)

  if(ARPACK_FOUND)
    message(STATUS "Using ARPACK ${ARPACK_LIBRARIES}")
    set(GOMA_TPL_LIBRARIES ${GOMA_TPL_LIBRARIES} ${ARPACK_LIBRARIES})
  else()
    message(STATUS "ARPACK library not found")
  endif()

else()
  message(STATUS "ARPACK library disabled")
endif()

option(ENABLE_SPARSE "ENABLE_SPARSE" ON)
if(ENABLE_SPARSE)
  set(Sparse_PREFIX
      "${Trilinos_DIR}/../../../"
      CACHE PATH "set sparse library location")
  find_package(Sparse)
  if(Sparse_FOUND)
    message(STATUS "Using Sparse ${Sparse_LIBRARIES}")
    set(GOMA_TPL_LIBRARIES ${GOMA_TPL_LIBRARIES} ${Sparse_LIBRARIES})
    list(APPEND GOMA_COMPILE_DEFINITIONS GOMA_ENABLE_SPARSE)
  else()
    message(STATUS "Sparse library not found")
    set(Sparse_INCLUDES "")
  endif()
else()
  message(STATUS "Sparse library disabled.")
  set(Sparse_INCLUDES "")
endif()

option(ENABLE_METIS "ENABLE_METIS" ON)
if(ENABLE_METIS)
  set(METIS_PREFIX
      "${Trilinos_DIR}/../../../"
      CACHE PATH "set METIS library location")
  find_package(METIS)
  if(METIS_FOUND)
    list(APPEND GOMA_COMPILE_DEFINITIONS GOMA_ENABLE_METIS)
    message(STATUS "Metis found: ${METIS_INCLUDES}")
    message(STATUS "Metis found: ${METIS_LIBRARIES}")
    set(GOMA_TPL_LIBRARIES ${GOMA_TPL_LIBRARIES} ${METIS_LIBRARIES})
    set(GOMA_TPL_INCLUDES ${GOMA_TPL_INCLUDES} ${METIS_INCLUDES})
  else()
    message(STATUS "metis library not found, disabling")
  endif()

else()
  message(STATUS "metis library disabled")
endif()

option(ENABLE_OMEGA_H "ENABLE_OMEGA_H" ON)
if(ENABLE_OMEGA_H)
  find_package(Omega_h 9)
  if(Omega_h_FOUND)
    set_target_properties(Omega_h::omega_h PROPERTIES INTERFACE_COMPILE_OPTIONS
                                                      "")
    set(GOMA_TPL_LIBRARIES ${GOMA_TPL_LIBRARIES} Omega_h::omega_h)
    list(APPEND GOMA_COMPILE_DEFINITIONS GOMA_ENABLE_OMEGA_H)
  else()
    message(STATUS "Omega_h not found")
  endif()
else()
  message(STATUS "omega_h library disabled")
endif()

if(ENABLE_PETSC AND PETSC_FOUND)
  set(GOMA_TPL_LIBRARIES ${GOMA_TPL_LIBRARIES} PkgConfig::PETSC)
  list(APPEND GOMA_COMPILE_DEFINITIONS GOMA_ENABLE_PETSC)
endif()

# UMFPACK in trilinos
list(FIND Trilinos_TPL_LIST UMFPACK trilinos_package_index)
option(ENABLE_UMFPACK "ENABLE_UMFPACK" ON)
if(ENABLE_UMFPACK)
  if(${trilinos_package_index} LESS 0)
    message(STATUS "Trilinos not built with UMFPACK TPL")
    message(STATUS "Checking for other umfpack")
    find_package(UMFPACK)
    if(UMFPACK_FOUND)
      message(STATUS "UMFPACK found")
    else()
      message(STATUS "UMFPACK library not found")
    endif()
    if(ENABLE_UMFPACK)
      message(STATUS "UMFPACK found")
      list(APPEND GOMA_COMPILE_DEFINITIONS GOMA_ENABLE_UMFPACK)
      set(GOMA_TPL_LIBRARIES ${GOMA_TPL_LIBRARIES} ${UMFPACK_LIBRARIES})
      set(GOMA_TPL_INCLUDES ${GOMA_TPL_INCLUDES} ${UMFPACK_INCLUDES})
    else()
      message(STATUS "UMFPACK library disabled")
    endif()
  else()
    message(STATUS "Trilinos built with UMFPACK TPL")
    message(STATUS "\tEnabling UMFPACK in Goma")
    list(APPEND GOMA_COMPILE_DEFINITIONS GOMA_ENABLE_UMFPACK)
  endif()
else()
  message(STATUS "UMFPACK disabled")
endif()

macro(option_definition OPTION VALUE)
  option(${OPTION} "${OPTION}" ${VALUE})
  if(${OPTION})
    list(APPEND GOMA_COMPILE_DEFINITIONS ${OPTION})
  endif()
endmacro()

if(ARPACK_FOUND)
  if(ENABLE_ARPACK)
    option(EIGEN_SERIAL "EIGEN_SERIAL" ON)
    option(EIGEN_PARALLEL "EIGEN_PARALLEL" OFF)
  else()
    option(EIGEN_SERIAL "EIGEN_SERIAL" OFF)
    option(EIGEN_PARALLEL "EIGEN_PARALLEL" OFF)
  endif()
  if(EIGEN_PARALLEL)
    if(NOT ARPACK_FOUND)
      message(FATAL_ERROR "Cannot use EIGEN_PARALLEL without ARPACK Library")
    endif()
    list(APPEND GOMA_COMPILE_DEFINITIONS EIGEN_PARALLEL)
    if(EIGEN_SERIAL)
      message(WARNING "Using EIGEN_PARALLEL instead of EIGEN_SERIAL")
    endif()
  elseif(EIGEN_SERIAL)
    if(NOT ${ARPACK_FOUND})
      message(FATAL_ERROR "Cannot use EIGEN_SERIAL without ARPACK Library")
    endif()
    list(APPEND GOMA_COMPILE_DEFINITIONS EIGEN_SERIAL)
  endif()
endif()

option(PRINT_STACK_TRACE_ON_EH "PRINT_STACK_TRACE_ON_EH" ON)
if(PRINT_STACK_TRACE_ON_EH)
  list(APPEND GOMA_COMPILE_DEFINITIONS PRINT_STACK_TRACE_ON_EH)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
endif()

if(ENABLE_APREPRO_LIB AND SEACASAprepro_FOUND)
  list(APPEND GOMA_COMPILE_DEFINITIONS GOMA_ENABLE_APREPRO_LIB)
endif()

option_definition(DISABLE_COLOR_ERROR_PRINT OFF)
option_definition(COMPILER_64BIT ON)
option_definition(GOMA_ENABLE_AMESOS ON)
option_definition(TRILINOS ON)
option_definition(CHECK_FINITE OFF)
option_definition(NO_CHEBYSHEV_PLEASE ON)
option_definition(PARALLEL ON)
option_definition(FP_EXCEPT OFF)
option_definition(ANALEIG_PLEASE OFF)
option_definition(DO_NO_UNROLL OFF)
option_definition(USERMAT_INITIALIZATION OFF)
option_definition(ENABLE_LOGGING OFF)
option_definition(ENABLE_KOMPLEX OFF)
set(MDE
    "27"
    CACHE STRING "set MDE")
list(APPEND GOMA_COMPILE_DEFINITIONS MDE=${MDE})
set(MAX_EXTERNAL_FIELD
    "4"
    CACHE STRING "set MAX_EXTERNAL_FIELD")
list(APPEND GOMA_COMPILE_DEFINITIONS MAX_EXTERNAL_FIELD=${MAX_EXTERNAL_FIELD})
set(MAX_CONC
    "4"
    CACHE STRING "set MAX_CONC")
list(APPEND GOMA_COMPILE_DEFINITIONS MAX_CONC=${MAX_CONC})
set(MAX_PROB_VAR
    "15"
    CACHE STRING "set MAX_PROB_VAR")
list(APPEND GOMA_COMPILE_DEFINITIONS MAX_PROB_VAR=${MAX_PROB_VAR})
set(MAX_NUMBER_MATLS
    "9"
    CACHE STRING "set MAX_NUMBER_MATLS")
list(APPEND GOMA_COMPILE_DEFINITIONS MAX_NUMBER_MATLS=${MAX_NUMBER_MATLS})

set(GOMA_INCLUDES
    include/adapt/omega_h_interface.h
    include/adapt/resetup_problem.h
    include/ac_conti.h
    include/ac_hunt.h
    include/ac_particles.h
    include/ac_stability.h
    include/ac_stability_util.h
    include/ac_update_parameter.h
    include/bc_colloc.h
    include/bc_contact.h
    include/bc_curve.h
    include/bc_dirich.h
    include/bc_integ.h
    include/bc/rotate.h
    include/bc_special.h
    include/bc_surfacedomain.h
    include/decomp_interface.h
    include/dp_comm.h
    include/dpi.h
    include/dp_map_comm_vec.h
    include/dp_types.h
    include/dp_utils.h
    include/dp_ghost.h
    include/dp_vif.h
    include/el_elm.h
    include/el_elm_info.h
    include/el_geom.h
    include/el_quality.h
    include/exo_conn.h
    include/exo_struct.h
    include/loca_const.h
    include/loca_util_const.h
    include/md_timer.h
    include/metis_decomp.h
    include/mm_as_alloc.h
    include/mm_as_const.h
    include/mm_as.h
    include/mm_as_structs.h
    include/mm_augc_util.h
    include/mm_bc.h
    include/mm_chemkin.h
    include/mm_dil_viscosity.h
    include/mm_elem_block_structs.h
    include/mm_fill_aux.h
    include/mm_fill_common.h
    include/mm_fill_em.h
    include/mm_fill_fill.h
    include/mm_fill.h
    include/mm_fill_jac.h
    include/mm_fill_ls.h
    include/mm_fill_population.h
    include/mm_fill_porous.h
    include/mm_fill_potential.h
    include/mm_fill_pthings.h
    include/mm_fill_ptrs.h
    include/mm_fill_rs.h
    include/mm_fill_shell.h
    include/mm_fill_solid.h
    include/mm_fill_species.h
    include/mm_fill_split.h
    include/mm_fill_stress.h
    include/mm_fill_terms.h
    include/mm_fill_util.h
    include/mm_fill_stabilization.h
    include/mm_flux.h
    include/mm_input.h
    include/mm_interface.h
    include/mm_more_utils.h
    include/mm_mp_const.h
    include/mm_mp.h
    include/mm_mp_structs.h
    include/mm_names.h
    include/mm_em_bc.h
    include/mm_ns_bc.h
    include/mm_numjac.h
    include/mm_post_def.h
    include/mm_post_proc.h
    include/mm_prob_def.h
    include/mm_qp_storage.h
    include/mm_qtensor_model.h
    include/mm_shell_bc.h
    include/mm_shell_util.h
    include/mm_sol_nonlinear.h
    include/mm_species.h
    include/mm_std_models.h
    include/mm_std_models_shell.h
    include/mm_unknown_map.h
    include/mm_viscosity.h
    include/rd_dpi.h
    include/rd_exo.h
    include/base_mesh.h
    include/rd_mesh.h
    include/rd_pixel_image.h
    include/rf_allo.h
    include/rf_bc_const.h
    include/rf_bc.h
    include/rf_element_storage_const.h
    include/rf_element_storage_struct.h
    include/rf_fem_const.h
    include/rf_fem.h
    include/rf_fill_const.h
    include/rf_io_const.h
    include/rf_io.h
    include/rf_io_structs.h
    include/rf_masks.h
    include/rf_mp.h
    include/rf_node_const.h
    include/rf_pre_proc.h
    include/rf_shape.h
    include/rf_solve.h
    include/rf_solver_const.h
    include/rf_solver.h
    include/rf_solve_segregated.h
    include/rf_util.h
    include/rf_vars_const.h
    include/bc/rotate_coordinates.h
    include/shell_tfmp_struct.h
    include/shell_tfmp_util.h
    include/sl_amesos_interface.h
    include/sl_aux.h
    include/sl_auxutil.h
    include/sl_aztecoo_interface.h
    include/sl_eggroll_def.h
    include/sl_eggroll.h
    include/sl_epetra_interface.h
    include/sl_epetra_util.h
    include/sl_lu.h
    include/sl_petsc.h
    include/sl_petsc_complex.h
    include/sl_matrix_util.h
    include/sl_stratimikos_interface.h
    include/sl_umf.h
    include/sl_util.h
    include/sl_util_structs.h
    include/std.h
    include/wr_dpi.h
    include/wr_exo.h
    include/wr_side_data.h
    include/wr_soln.h
    include/brkfix/bbb.h
    include/brkfix/exo_utils.h
    include/brkfix/fix.h
    include/brkfix/setup_fix_data.h)

set(GOMA_SOURCES
    src/adapt/omega_h_interface.cpp
    src/adapt/resetup_problem.c
    src/ac_conti.c
    src/ac_display_parameter.c
    src/ac_hunt.c
    src/ac_loca_interface.c
    src/ac_particles.c
    src/ac_stability.c
    src/ac_stability_util.c
    src/ac_update_parameter.c
    src/bc_colloc.c
    src/bc_contact.c
    src/bc_curve.c
    src/bc_dirich.c
    src/bc_integ.c
    src/bc/rotate.c
    src/bc_special.c
    src/bc_surfacedomain.c
    src/decomp_interface.c
    src/dp_comm.c
    src/dp_map_comm_vec.c
    src/dp_utils.c
    src/dp_ghost.cpp
    src/dp_vif.c
    src/el_elm_info.c
    src/el_quality.c
    src/exo_conn.c
    src/globals.c
    src/loca_bord.c
    src/loca_eigen_c2f.F
    src/loca_eigen_cayley.F
    src/loca_eigenvalue.c
    src/loca_lib.c
    src/loca_util.c
    src/md_ieee.c
    src/md_timer.c
    src/metis_decomp.c
    src/mm_as_alloc.c
    src/mm_augc_util.c
    src/mm_bc.c
    src/mm_bc_conflict.c
    src/mm_chemkin.c
    src/mm_dil_viscosity.c
    src/mm_fill_aux.c
    src/mm_fill.c
    src/mm_fill_common.c
    src/mm_fill_em.c
    src/mm_fill_fill.c
    src/mm_fill_interface.c
    src/mm_fill_jac.c
    src/mm_fill_ls.c
    src/mm_fill_population.c
    src/mm_fill_porous.c
    src/mm_fill_potential.c
    src/mm_fill_pthings.c
    src/mm_fill_ptrs.c
    src/mm_fill_rs.c
    src/mm_fill_shell.c
    src/mm_fill_solid.c
    src/mm_fill_species.c
    src/mm_fill_stress.c
    src/mm_fill_terms.c
    src/mm_fill_split.c
    src/mm_fill_util.c
    src/mm_fill_stabilization.c
    src/mm_flux.c
    src/mm_input_bc.c
    src/mm_input.c
    src/mm_input_mp.c
    src/mm_input_particles.c
    src/mm_input_util.c
    src/mm_interface.c
    src/mm_matrl.c
    src/mm_more_utils.c
    src/mm_em_bc.c
    src/mm_ns_bc.c
    src/mm_numjac.c
    src/mm_placid.c
    src/mm_porous_EOS.c
    src/mm_post_proc.c
    src/mm_post_proc_util.c
    src/mm_prob_def.c
    src/mm_propertyJac.c
    src/mm_qp_storage.c
    src/mm_qtensor_model.c
    src/mm_shell_bc.c
    src/mm_shell_util.c
    src/mm_sol_nonlinear.c
    src/mm_species.c
    src/mm_std_models.c
    src/mm_std_models_shell.c
    src/mm_unknown_map.c
    src/mm_viscosity.c
    src/rd_dpi.c
    src/rd_exo.c
    src/base_mesh.c
    src/rd_mesh.c
    src/rd_pixel_image2.c
    src/rd_pixel_image.c
    src/rf_allo.c
    src/rf_element_storage.c
    src/rf_node.c
    src/rf_node_vars.c
    src/rf_pre_proc.c
    src/rf_setup_problem.c
    src/rf_shape.c
    src/rf_solve.c
    src/rf_solve_segregated.c
    src/rf_util.c
    src/rf_vars.c
    src/bc/rotate_coordinates.c
    src/shell_tfmp_util.c
    src/sl_amesos_interface.cpp
    src/sl_aux.c
    src/sl_auxutil.c
    src/sl_aztecoo_interface.cpp
    src/sl_eggroll01.c
    src/sl_eggroll02.c
    src/sl_eggroll03.c
    src/sl_eggroll04.c
    src/sl_eggroll05.c
    src/sl_eggrollutil.c
    src/sl_eggrollwrap.c
    src/sl_epetra_interface.cpp
    src/sl_epetra_util.cpp
    src/sl_leastsquares.cpp
    src/sl_lu.c
    src/sl_lu_fill.c
    src/sl_lustat.c
    src/sl_ma28.c
    src/sl_matrix_dump.c
    src/sl_petsc.c
    src/sl_petsc_pcd.c
    src/sl_petsc_complex.c
    src/sl_matrix_util.c
    src/sl_squash.c
    src/sl_stratimikos_interface.cpp
    src/sl_umf.c
    src/sl_util.c
    src/wr_dpi.c
    src/wr_exo.c
    src/wr_side_data.c
    src/wr_soln.c
    src/brkfix/bbb.c
    src/brkfix/setup_fix_data.cpp
    src/brkfix/fix_exo_file.c)

set(GOMA_USER_INCLUDES
    include/user_bc.h
    include/user_continuation.h
    include/user_mp_gen.h
    include/user_mp.h
    include/user_out_hkm.h
    include/user_post.h
    include/user_pre.h
    include/user_shell.h
    include/usr_print.h)

set(GOMA_USER_SOURCES
    src/user_ac.c
    src/user_bc.c
    src/user_continuation.c
    src/user_mp.c
    src/user_mp_gen.c
    src/user_out_hkm.c
    src/user_post.c
    src/user_pre.c
    src/user_print.c
    src/user_senkin.F)

set(GOMA_UTIL_INCLUDES include/bc/rotate_util.h include/mm_eh.h
                       include/util/goma_normal.h include/util/aprepro_helper.h)

set(GOMA_UTIL_SOURCES src/bc/rotate_util.c src/util/goma_normal.c src/mm_eh.c
                      src/util/aprepro_helper.cpp)

set(GDS_INCLUDES include/gds/gds_vector.h)

set(GDS_SOURCES src/gds/gds_vector.c)

# Set a default build type for single-configuration CMake generators if no build
# type is set. https://blog.kitware.com/cmake-and-the-default-build-type/
set(default_build_type "RelWithDebInfo")
get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

if(isMultiConfig)
  if(NOT "Asan" IN_LIST CMAKE_CONFIGURATION_TYPES)
    list(APPEND CMAKE_CONFIGURATION_TYPES Asan)
  endif()
  if(NOT "Gprof" IN_LIST CMAKE_CONFIGURATION_TYPES)
    list(APPEND CMAKE_CONFIGURATION_TYPES Gprof)
  endif()
else()
  set(allowedBuildTypes Asan Gprof Debug Release RelWithDebInfo MinSizeRel)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "${allowedBuildTypes}")

  if(CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE IN_LIST allowedBuildTypes)
    message(FATAL_ERROR "Invalid build type: ${CMAKE_BUILD_TYPE}")
  endif()
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${default_build_type}"
      CACHE STRING "Choose the type of build." FORCE)
endif()

if(${CMAKE_BUILD_TYPE} STREQUAL "Asan")
  set(CMAKE_C_FLAGS_ASAN
      "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer"
      CACHE STRING
            "Flags used by the C compiler for Asan build type or configuration."
            FORCE)
  mark_as_advanced(CMAKE_C_FLAGS_ASAN)

  set(CMAKE_CXX_FLAGS_ASAN
      "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer"
      CACHE
        STRING
        "Flags used by the C++ compiler for Asan build type or configuration."
        FORCE)
  mark_as_advanced(CMAKE_CXX_FLAGS_ASAN)

  set(CMAKE_Fortran_FLAGS_ASAN
      "${CMAKE_Fortran_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer"
      CACHE
        STRING
        "Flags used by the Fortran compiler for Asan build type or configuration."
        FORCE)
  mark_as_advanced(CMAKE_Fortran_FLAGS_ASAN)

  set(CMAKE_EXE_LINKER_FLAGS_ASAN
      "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -fsanitize=address"
      CACHE STRING
            "Linker flags to be used to create executables for Asan build type."
            FORCE)
  mark_as_advanced(CMAKE_EXE_LINKER_FLAGS_ASAN)

  set(CMAKE_SHARED_LINKER_FLAGS_ASAN
      "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -fsanitize=address"
      CACHE
        STRING
        "Linker lags to be used to create shared libraries for Asan build type."
        FORCE)
  mark_as_advanced(CMAKE_SHARED_LINKER_FLAGS_ASAN)
endif()

if(${CMAKE_BUILD_TYPE} STREQUAL "Gprof")
  set(CMAKE_C_FLAGS_GPROF
      "${CMAKE_C_FLAGS_DEBUG} -pg"
      CACHE
        STRING
        "Flags used by the C compiler for Gprof build type or configuration."
        FORCE)
  mark_as_advanced(CMAKE_C_FLAGS_GPROF)

  set(CMAKE_CXX_FLAGS_GPROF
      "${CMAKE_CXX_FLAGS_DEBUG} -pg"
      CACHE
        STRING
        "Flags used by the C++ compiler for Gprof build type or configuration."
        FORCE)
  mark_as_advanced(CMAKE_CXX_FLAGS_GPROF)

  set(CMAKE_Fortran_FLAGS_GPROF
      "${CMAKE_Fortran_FLAGS_DEBUG} -pg"
      CACHE
        STRING
        "Flags used by the Fortran compiler for Gprof build type or configuration."
        FORCE)
  mark_as_advanced(CMAKE_Fortran_FLAGS_GPROF)

  set(CMAKE_EXE_LINKER_FLAGS_GPROF
      "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -pg"
      CACHE
        STRING
        "Linker flags to be used to create executables for Gprof build type."
        FORCE)
  mark_as_advanced(CMAKE_EXE_LINKER_FLAGS_GPROF)

  set(CMAKE_SHARED_LINKER_FLAGS_GPROF
      "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -pg"
      CACHE
        STRING
        "Linker lags to be used to create shared libraries for Gprof build type."
        FORCE)
  mark_as_advanced(CMAKE_SHARED_LINKER_FLAGS_GPROF)
endif()

option(ENABLE_WARNINGS "ENABLE_WARNINGS" ON)

execute_process(
  COMMAND git describe --dirty --always --tags
  OUTPUT_VARIABLE goma_git_version
  RESULT_VARIABLE goma_using_git
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  ERROR_QUIET)

if(goma_using_git EQUAL 0)
  string(REGEX REPLACE "\n$" "" goma_git_version "${goma_git_version}")
else()
  unset(goma_git_version)
endif()
if(goma_git_version)
  list(APPEND GOMA_COMPILE_DEFINITIONS GIT_VERSION=${goma_git_version})
  message(STATUS "Git found setting version: ${goma_git_version}")
endif()

message(STATUS "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
message(STATUS "Trilinos_C_COMPILER=${Trilinos_C_COMPILER}")
message(STATUS "Trilinos_CXX_COMPILER=${Trilinos_CXX_COMPILER}")
message(STATUS "Trilinos_Fortran_COMPILER=${Trilinos_Fortran_COMPILER}")
if(ENABLE_PETSC AND PETSC_FOUND)
  message(STATUS "PETSC_C_COMPILER=${PETSC_C_COMPILER}")
  message(STATUS "PETSC_CXX_COMPILER=${PETSC_CXX_COMPILER}")
  message(STATUS "PETSC_FORTRAN_COMPILER=${PETSC_FORTRAN_COMPILER}")
endif()

message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
if(${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
  message(STATUS "CMAKE_C_FLAGS_RELWITHDEBINFO=${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  message(
    STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  message(
    STATUS
      "CMAKE_Fortran_FLAGS_RELWITHDEBINFO=${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
  )
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Release")
  message(STATUS "CMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE}")
  message(STATUS "CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}")
  message(STATUS "CMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  message(STATUS "CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
  message(STATUS "CMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}")
  message(STATUS "CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Asan")
  message(STATUS "CMAKE_C_FLAGS_ASAN=${CMAKE_C_FLAGS_ASAN}")
  message(STATUS "CMAKE_CXX_FLAGS_ASAN=${CMAKE_CXX_FLAGS_ASAN}")
  message(STATUS "CMAKE_Fortran_FLAGS_ASAN=${CMAKE_Fortran_FLAGS_ASAN}")
elseif(${CMAKE_BUILD_TYPE} STREQUAL "Gprof")
  message(STATUS "CMAKE_C_FLAGS_GPROF=${CMAKE_C_FLAGS_GPROF}")
  message(STATUS "CMAKE_CXX_FLAGS_GPROF=${CMAKE_CXX_FLAGS_GPROF}")
  message(STATUS "CMAKE_Fortran_FLAGS_GPROF=${CMAKE_Fortran_FLAGS_GPROF}")
endif()
message("\n")
message(STATUS "Goma compile defines:")
foreach(d ${GOMA_COMPILE_DEFINITIONS})
  message(STATUS "${d}")
endforeach()
message("\n")

set(GOMA_TPL_INCLUDES
    ${GOMA_TPL_INCLUDES}
    ${Sparse_INCLUDES}
    ${SEACASExodus_INCLUDE_DIRS}
    ${SEACASExodus_TPL_INCLUDE_DIRS}
    ${SEACASAprepro_INCLUDE_DIRS}
    ${SEACASAprepro_TPL_INCLUDE_DIRS}
    ${Trilinos_INCLUDE_DIRS}
    ${Trilinos_TPL_INCLUDE_DIRS}
    ${MPI_C_INCLUDE_PATH})

set(GOMA_TPL_LIBRARY_DIRS
    ${SEACASExodus_LIBRARY_DIRS} ${SEACASExodus_TPL_LIBRARY_DIRS}
    ${SEACASAprepro_LIBRARY_DIRS} ${SEACASAprepro_TPL_LIBRARY_DIRS}
    ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
if(ENABLE_PETSC AND PETSC_FOUND)
  set(GOMA_TPL_LIBRARY_DIRS ${GOMA_TPL_LIBRARY_DIRS}
                            ${PETSC_STATIC_LIBRARY_DIRS})
endif()

set(GOMA_TPL_LIBRARIES
    ${SEACASExodus_LIBRARIES}
    ${SEACASExodus_TPL_LIBRARIES}
    ${SEACASAprepro_LIBRARIES}
    ${SEACASAprepro_TPL_LIBRARIES}
    ${GOMA_TPL_LIBRARIES}
    ${Trilinos_LIBRARIES}
    ${Trilinos_TPL_LIBRARIES}
    ${Trilinos_EXTRA_LD_FLAGS}
    ${MPI_C_LIBRARIES}
    ${MPI_Fortran_LIBRARIES})
if(ENABLE_PETSC AND PETSC_FOUND)
  set(GOMA_TPL_LIBRARIES ${GOMA_TPL_LIBRARIES} ${PETSC_STATIC_LIBRARIES})
endif()

# remove some C++ files from unity build as MACROS conflict with Trilinos
# headers
set_source_files_properties(
  src/dp_ghost.cpp src/sl_epetra_util.cpp src/adapt/omega_h_interface.cpp
  PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)

add_library(goma_user ${GOMA_USER_SOURCES} ${GOMA_USER_INCLUDES})
target_include_directories(goma_user PUBLIC include)
target_include_directories(goma_user SYSTEM PRIVATE ${GOMA_TPL_INCLUDES})
target_compile_definitions(goma_user PUBLIC ${GOMA_COMPILE_DEFINITIONS})
if(ENABLE_WARNINGS)
  target_compile_options(
    goma_user
    PRIVATE
      $<$<AND:$<C_COMPILER_ID:Clang>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter>
      $<$<AND:$<C_COMPILER_ID:GNU>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter
      -Wimplicit-fallthrough>)
endif()

add_library(goma ${GOMA_SOURCES} ${GOMA_INCLUDES})
target_include_directories(goma PUBLIC include)
target_include_directories(goma SYSTEM PRIVATE ${GOMA_TPL_INCLUDES})
target_link_directories(goma PUBLIC ${GOMA_TPL_LIBRARY_DIRS})
target_link_libraries(goma PUBLIC ${GOMA_TPL_LIBRARIES} gds)
target_compile_definitions(goma PUBLIC ${GOMA_COMPILE_DEFINITIONS})
if(ENABLE_WARNINGS)
  target_compile_options(
    goma
    PRIVATE
      $<$<AND:$<C_COMPILER_ID:Clang>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter>
      $<$<AND:$<C_COMPILER_ID:GNU>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter
      -Wimplicit-fallthrough>)
endif()

add_library(gds ${GDS_SOURCES} ${GDS_INCLUDES})
target_include_directories(gds PUBLIC include)
if(ENABLE_WARNINGS)
  target_compile_options(
    gds
    PRIVATE
      $<$<AND:$<C_COMPILER_ID:Clang>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter>
      $<$<AND:$<C_COMPILER_ID:GNU>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter
      -Wimplicit-fallthrough>)
endif()

add_library(goma_util ${GOMA_UTIL_SOURCES} ${GOMA_UTIL_INCLUDES})
target_include_directories(goma_util PUBLIC include)
target_include_directories(goma_util SYSTEM PRIVATE ${GOMA_TPL_INCLUDES})
target_link_libraries(goma_util PUBLIC ${GOMA_TPL_LIBRARIES} gds)
target_link_directories(goma_util PUBLIC ${GOMA_TPL_LIBRARY_DIRS})
target_compile_definitions(goma_util PUBLIC ${GOMA_COMPILE_DEFINITIONS})
if(ENABLE_WARNINGS)
  target_compile_options(
    goma_util
    PRIVATE
      $<$<AND:$<C_COMPILER_ID:Clang>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter>
      $<$<AND:$<C_COMPILER_ID:GNU>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter
      -Wimplicit-fallthrough>)
endif()

add_executable(goma_exe src/main.c)
set_target_properties(goma_exe PROPERTIES OUTPUT_NAME "goma")
target_include_directories(goma_exe PUBLIC include)
target_include_directories(goma_exe SYSTEM PRIVATE ${GOMA_TPL_INCLUDES})
target_link_libraries(goma_exe PUBLIC ${GOMA_TPL_LIBRARIES} goma goma_user
                                      goma_util)
target_compile_definitions(goma_exe PUBLIC ${GOMA_COMPILE_DEFINITIONS})
if(ENABLE_WARNINGS)
  target_compile_options(
    goma_exe
    PRIVATE
      $<$<AND:$<C_COMPILER_ID:Clang>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter>
      $<$<AND:$<C_COMPILER_ID:GNU>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter
      -Wimplicit-fallthrough>)
endif()

add_executable(fix_exe src/fix_main.c)
set_target_properties(fix_exe PROPERTIES OUTPUT_NAME "fix")
target_include_directories(fix_exe PUBLIC include)
target_include_directories(fix_exe SYSTEM PRIVATE ${GOMA_TPL_INCLUDES})
target_link_libraries(fix_exe PUBLIC ${GOMA_TPL_LIBRARIES} goma goma_user
                                     goma_util)
target_compile_definitions(fix_exe PUBLIC ${GOMA_COMPILE_DEFINITIONS})
if(ENABLE_WARNINGS)
  target_compile_options(
    fix_exe
    PRIVATE
      $<$<AND:$<C_COMPILER_ID:Clang>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter>
      $<$<AND:$<C_COMPILER_ID:GNU>,$<OR:$<COMPILE_LANGUAGE:C>,$<COMPILE_LANGUAGE:CXX>>>:-Wall
      -Wextra
      -Wno-unused-parameter
      -Wimplicit-fallthrough>)
endif()

install(
  TARGETS goma_exe fix_exe
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include/goma)

option(ENABLE_TESTING "ENABLE_TESTING" OFF)
if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()
