cmake_minimum_required (VERSION 3.13)
cmake_policy(VERSION 3.13)

project(goma)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

set(PROJECT_NAME goma)
list(APPEND CMAKE_MODULE_PATH "${goma_SOURCE_DIR}/cmake")

set(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})

# Get Trilinos as one entity
FIND_PACKAGE(Trilinos REQUIRED)
MESSAGE(STATUS "Trilinos_DIR = ${Trilinos_DIR}")
#message(STATUS "Trilinos LD Flags = ${Trilinos_EXTRA_LD_FLAGS}")

set(MPI_C_COMPILER ${Trilinos_C_COMPILER})
find_package(MPI REQUIRED)

# Avoid linking to too many libraries by removing duplicates
# reverse so that correct order for linking is preserved
list(REVERSE Trilinos_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
list(REVERSE Trilinos_LIBRARIES)

list(REVERSE Trilinos_TPL_LIBRARIES)
list(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list(REVERSE Trilinos_TPL_LIBRARIES)

set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} ${MPI_C_INCLUDE_PATH})

# Trilinos required packages
foreach(package Amesos AztecOO Epetra)
  LIST(FIND Trilinos_PACKAGE_LIST ${package} trilinos_package_index)
  if(${trilinos_package_index} LESS 0)
  message(FATAL_ERROR "Trilinos not built with ${package}")
  endif()
endforeach()

# Trilinos required TPLs
foreach(package MPI UMFPACK)
  LIST(FIND Trilinos_TPL_LIST ${package} trilinos_package_index)
  if(${trilinos_package_index} LESS 0)
    message(FATAL_ERROR "Trilinos not built with TPL ${package}")
  endif()
endforeach()

# Trilinos extra TPLs
foreach(package SuperLUDist MUMPS)
  LIST(FIND Trilinos_TPL_LIST ${package} trilinos_package_index)
  if(${trilinos_package_index} LESS 0)
    message(WARNING "Trilinos not built with TPL ${package} some features of Goma won't work")
  endif()
endforeach()

# UMFPACK in trilinos
LIST(FIND Trilinos_TPL_LIST UMFPACK trilinos_package_index)
if(${trilinos_package_index} LESS 0)
  message(STATUS "Trilinos not built with UMFPACK TPL")
else()
  message(STATUS "Trilinos built with UMFPACK TPL")
  message(STATUS "\tEnabling UMFPACK in Goma")
  list(APPEND GOMA_COMPILE_DEFINITIONS HAVE_UMFPACK)
endif()

# stratimikos and teko
LIST(FIND Trilinos_PACKAGE_LIST Stratimikos stratimikos_package_index)
if(${stratimikos_package_index} GREATER_EQUAL 0)
  message(STATUS "TRILINOS: Stratimikos found, enabling in Goma")
  list(APPEND GOMA_COMPILE_DEFINITIONS HAVE_STRATIMIKOS)
  LIST(FIND Trilinos_PACKAGE_LIST Teko teko_package_index)
  if(${teko_package_index} GREATER_EQUAL 0)
    message(STATUS "TRILINOS: Teko found, enabling in Goma")
    list(APPEND GOMA_COMPILE_DEFINITIONS HAVE_TEKO)
  endif()
endif()

# Exodus
find_package(SEACASExodus REQUIRED HINTS ${Trilinos_DIR}/../SEACASExodus)
MESSAGE(STATUS "SEACASExodus_DIR = ${SEACASExodus_DIR}")
find_package(SEACASChaco REQUIRED HINTS ${Trilinos_DIR}/../SEACASChaco)
MESSAGE(STATUS "SEACASChaco_DIR = ${SEACASChaco_DIR}")

message("")

set(GOMA_TPL_LIBRARIES "")

set(ARPACK_PREFIX "${Trilinos_DIR}/../../../../ARPACK" CACHE PATH "set ARPACK library location")

find_package(ARPACK)

if(ARPACK_FOUND)
  message(STATUS "Using ARPACK ${ARPACK_LIBRARIES}")
  set(GOMA_TPL_LIBRARIES ${GOMA_TPL_LIBRARIES} ${ARPACK_LIBRARIES})
else()
  message(STATUS "ARPACK library not found")
endif()

set(Sparse_PREFIX "${Trilinos_DIR}/../../../../sparse" CACHE PATH "set sparse library location")
find_package(Sparse)
if (Sparse_FOUND)
  message(STATUS "Using Sparse ${Sparse_LIBRARIES}")
  set(GOMA_TPL_LIBRARIES ${GOMA_TPL_LIBRARIES} ${Sparse_LIBRARIES})
  list(APPEND GOMA_COMPILE_DEFINITIONS HAVE_SPARSE)
else()
  message(STATUS "Sparse library not found")
  set(Sparse_INCLUDES "")
endif()

macro(option_definition OPTION VALUE)
  option(${OPTION} "${OPTION}" ${VALUE})
  if(${OPTION})
    list(APPEND GOMA_COMPILE_DEFINITIONS ${OPTION})
  endif()
endmacro()

option_definition(COMPILER_64BIT ON)
option_definition(ENABLE_AMESOS ON)
option_definition(TRILINOS ON)
option_definition(CHECK_FINITE ON)
option_definition(NO_CHEBYSHEV_PLEASE ON)
option_definition(COUPLED_FILL ON)
option_definition(PARALLEL ON)
option_definition(FP_EXCEPT OFF)
if (ARPACK_FOUND)
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

option(PRINT_STACK_TRACE_ON_EH "PRINT_STACK_TRACE_ON_EH" ON)
if(PRINT_STACK_TRACE_ON_EH)
  list(APPEND GOMA_COMPILE_DEFINITIONS PRINT_STACK_TRACE_ON_EH)
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} -rdynamic)
endif()

option_definition(DISABLE_COLOR_ERROR_PRINT OFF)

set(MDE "27" CACHE STRING "set MDE")
list(APPEND GOMA_COMPILE_DEFINITIONS MDE=${MDE})
set(MAX_EXTERNAL_FIELD "4" CACHE STRING "set MAX_EXTERNAL_FIELD")
list(APPEND GOMA_COMPILE_DEFINITIONS MAX_EXTERNAL_FIELD=${MAX_EXTERNAL_FIELD})
set(MAX_CONC "4" CACHE STRING "set MAX_CONC")
list(APPEND GOMA_COMPILE_DEFINITIONS MAX_CONC=${MAX_CONC})
set(MAX_PROB_VAR "15" CACHE STRING "set MAX_PROB_VAR")
list(APPEND GOMA_COMPILE_DEFINITIONS MAX_PROB_VAR=${MAX_PROB_VAR})

set(GOMA_INCLUDES
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
        include/brkfix/
        include/brk_utils.h
        include/dp_comm.h
        include/dpi.h
        include/dp_map_comm_vec.h
        include/dp_types.h
        include/dp_utils.h
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
        include/mm_fill_stress.h
        include/mm_fill_terms.h
        include/mm_fill_util.h
        include/mm_flux.h
        include/mm_input.h
        include/mm_interface.h
        include/mm_more_utils.h
        include/mm_mp_const.h
        include/mm_mp.h
        include/mm_mp_structs.h
        include/mm_names.h
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
        include/rf_node.h
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
        include/sl_matrix_util.h
        include/sl_rcm.h
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
        include/brkfix/brkfix.h
        include/brkfix/brkfix_types.h
        include/brkfix/brk.h
        include/brkfix/emuck.h
        include/brkfix/exo_utils.h
        include/brkfix/fix.h
        include/brkfix/mk_dm.h
        include/brkfix/nodesc.h
        include/brkfix/ppi.h
        include/brkfix/rd_in.h
        include/brkfix/sam_perea.h
        include/brkfix/utils.h
        include/brkfix/wr_coords.h
        include/brkfix/wr_graph_file.h
        include/gds/gds_vector.h)

set(GOMA_SOURCES
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
        src/brkfix/
        src/brk_utils.c
        src/dp_comm.c
        src/dp_map_comm_vec.c
        src/dp_utils.c
        src/dp_vif.c
        src/el_elm_info.c
        src/el_quality.c
        src/exo_conn.c
        src/loca_bord.c
        src/loca_eigen_c2f.F
        src/loca_eigen_cayley.F
        src/loca_eigenvalue.c
        src/loca_lib.c
        src/loca_util.c
        src/md_ieee.c
        src/md_timer.c
        src/mm_as_alloc.c
        src/mm_augc_util.c
        src/mm_bc.c
        src/mm_bc_conflict.c
        src/mm_chemkin.c
        src/mm_dil_viscosity.c
        src/mm_fill_aux.c
        src/mm_fill.c
        src/mm_fill_common.c
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
        src/mm_fill_util.c
        src/mm_flux.c
        src/mm_input_bc.c
        src/mm_input.c
        src/mm_input_mp.c
        src/mm_input_particles.c
        src/mm_input_util.c
        src/mm_interface.c
        src/mm_matrl.c
        src/mm_more_utils.c
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
        src/sl_amesos_interface.C
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
        src/sl_front_setup.c
        src/sl_leastsquares.C
        src/sl_lu.c
        src/sl_lu_fill.c
        src/sl_lustat.c
        src/sl_ma28.c
        src/sl_matrix_dump.c
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
        src/brkfix/brk_exo_file.c
        src/brkfix/emuck.c
        src/brkfix/fix_exo_file.c
        src/brkfix/mk_dm.c
        src/brkfix/ppi.c
        src/brkfix/rd_in.c
        src/brkfix/sam_perea.c
        src/brkfix/utils.c
        src/brkfix/wr_coords.c
        src/brkfix/wr_graph_file.c)

set(GOMA_USER_INCLUDES
        include/user_bc.h
        include/user_continuation.h
        include/user_mp_gen.h
        include/user_mp.h
        include/user_out_hkm.h
        include/user_post.h
        include/user_pre.h
        include/user_shell.h
        include/usr_print.h
)
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
        src/user_senkin.F
)

set(GOMA_UTIL_INCLUDES
  include/bc/rotate_util.h
  include/mm_eh.h
)

set(GOMA_UTIL_SOURCES
  src/bc/rotate_util.c
  src/mm_eh.c
)

set(GDS_INCLUDES
        include/gds/gds_vector.h)

set(GDS_SOURCES
        src/gds/gds_vector.c)

# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
set(GOMA_WARNINGS "-Wall -Wextra -Wno-unused-parameter" CACHE STRING "set warnings for goma C and CXX")
set(GOMA_DEBUG_FLAGS "-O0 -ggdb3" CACHE STRING "set debug flags for goma C and CXX")
set(GOMA_RELEASE_FLAGS "-O2 -DNDEBUG" CACHE STRING "set release flags for goma C and CXX")

# Select flags.
SET(CMAKE_CXX_FLAGS ${GOMA_WARNINGS})
SET(CMAKE_CXX_FLAGS_RELEASE "${GOMA_RELEASE_FLAGS}")
SET(CMAKE_CXX_FLAGS_DEBUG  "${GOMA_DEBUG_FLAGS}")
SET(CMAKE_C_FLAGS ${GOMA_WARNINGS})
SET(CMAKE_C_FLAGS_RELEASE "${GOMA_RELEASE_FLAGS}")
SET(CMAKE_C_FLAGS_DEBUG  "${GOMA_DEBUG_FLAGS}")

message(STATUS "CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}")
message("\n")
message(STATUS "Goma compile defines:")
foreach(d ${GOMA_COMPILE_DEFINITIONS})
  message(STATUS "${d}")
endforeach()
message("\n")

set(GOMA_TPL_INCLUDES
  ${Sparse_INCLUDES}
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${SEACASExodus_INCLUDE_DIRS}
  ${SEACASChaco_INCLUDE_DIRS}
  ${SEACASExodus_TPL_INCLUDE_DIRS}
  ${SEACASChaco_TPL_INCLUDE_DIRS}
)

set(GOMA_TPL_LIBRARY_DIRS
  ${Trilinos_LIBRARY_DIRS}
  ${Trilinos_TPL_LIBRARY_DIRS}
  ${SEACASExodus_LIBRARY_DIRS}
  ${SEACASExodus_TPL_LIBRARY_DIRS}
  ${SEACASChaco_LIBRARY_DIRS}
  ${SEACASChaco_TPL_LIBRARY_DIRS}
)

set(GOMA_TPL_LIBRARIES
  ${GOMA_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_EXTRA_LD_FLAGS}
  ${SEACASExodus_LIBRARIES}
  ${SEACASExodus_TPL_LIBRARIES}
  ${SEACASChaco_LIBRARIES}
  ${SEACASChaco_TPL_LIBRARIES}
)

add_library(goma_user ${GOMA_USER_SOURCES} ${GOMA_USER_INCLUDES})
target_include_directories(goma_user PUBLIC include)
target_compile_definitions(goma_user PUBLIC ${GOMA_COMPILE_DEFINITIONS})

add_library(goma ${GOMA_SOURCES} ${GOMA_INCLUDES})
target_include_directories(goma PUBLIC include)
target_include_directories(goma PRIVATE include ${GOMA_TPL_INCLUDES})
target_link_directories(goma PUBLIC ${GOMA_TPL_LIBRARY_DIRS})
target_link_libraries(goma ${GOMA_TPL_LIBRARIES} gds)
target_compile_definitions(goma PUBLIC ${GOMA_COMPILE_DEFINITIONS})

add_library(gds ${GDS_SOURCES} ${GDS_INCLUDES})
target_include_directories(gds PUBLIC include)

add_library(goma_util ${GOMA_UTIL_SOURCES} ${GOMA_UTIL_INCLUDES})
target_include_directories(goma_util PUBLIC include)
target_link_libraries(goma_util gds)
target_compile_definitions(goma_util PUBLIC ${GOMA_COMPILE_DEFINITIONS})

add_executable(goma_exe src/main.c)
set_target_properties(goma_exe PROPERTIES OUTPUT_NAME "goma")
target_include_directories(goma_exe PUBLIC include PRIVATE ${GOMA_TPL_INCLUDES})
target_link_libraries(goma_exe PUBLIC ${GOMA_TPL_LIBRARIES} goma goma_user goma_util gds)
target_compile_definitions(goma_exe PUBLIC ${GOMA_COMPILE_DEFINITIONS})

install(TARGETS goma_exe #gds goma_util goma goma_user
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/goma
)

option(ENABLE_TESTING "ENABLE_TESTING" OFF)
if(ENABLE_TESTING)
  add_subdirectory(tests)
endif()
