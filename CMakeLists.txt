#/************************************************************************ *
#* Goma - Multiphysics finite element software                             *
#* Sandia National Laboratories                                            *
#*                                                                         *
#* Copyright (c) 2016 Sandia Corporation.                                  *
#*                                                                         *
#* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
#* the U.S. Government retains certain rights in this software.            *
#*                                                                         *
#* This software is distributed under the GNU General Public License?      *
#\************************************************************************/
#
# CMake file for Goma 6.1
# #######################
#
# If you are a developer on this project: Do not commit any files generated
# by this script to the repository or everything will (might) break.
#

cmake_minimum_required (VERSION 3.0)
cmake_policy(VERSION 3.0)
cmake_policy(SET CMP0053 OLD)

#Party like it's 1977
enable_language(Fortran)


include(${CMAKE_SOURCE_DIR}/ProjectName.cmake)
include(${CMAKE_SOURCE_DIR}/CMakeScripts.cmake)


#It is required that you name the project in the main file in addition to an include
project(${PROJECT_NAME} NONE)



#Run "ccmake ." to change any options after running "cmake ."
option(Debug_Mode "This enables the debug flag in make, building it with Debug enabled." OFF)

#Turn off if you used the "noc++11 script when installing libraries
#option(Goma_ENABLE_CXX11_DEFAULT "Use the CXX11 flag. Required if you used the normal build script" ON)
option(Has_Stratimikos "If you have stratimikos and want support for it enable this. This modifies the definitions" OFF)


set(User_Define "")
set(User_Flags "")

if(Has_Stratimikos)
    set(User_Define "${User_Define} -DHAVE_STRATIMIKOS")
endif()
if(Debug_Mode)
    set(User_Flags "${User_Flags} -Wall")
    set(DCMAKE_BUILD_TYPE "Debug")
endif()


set(User_Define "${User_Define} -Dlinux -DCOMPILER_64BIT -DENABLE_AMESOS -DTRILINOS -DCHECK_FINITE -DNO_CHEBYSHEV_PLEASE -DMDE=27 -DMAX_PROB_VAR=15 -DMAX_EXTERNAL_FIELD=4 -DMAX_CONC=4 -DCOUPLED_FILL -DPARALLEL -DEIGEN_SERIAL" CACHE STRING "Enter or modify the compiler definitions. This overrides any definition related options.")

set(User_Flags "${User_Flags}" CACHE STRING "Enter any compiler flags you want to run on all three compilers. This overrides any cache related options.")
set(User_C_Flags "-w" CACHE STRING "Enter any extra compiler flags you want for all C compilation.")
set(User_CXX_Flags "-w" CACHE STRING "Enter any extra compiler flags you want for all C++ compilation.")
set(User_F_Flags "" CACHE STRING "Enter any extra compiler flags you want for all Fortran compilation.")
set(User_LD_Flags "" CACHE STRING "Enter any extra compiler flags you want when linking files and creating the Goma executable.")



#By default no custom definitions are used. These definitions are in addition to the default definitions
#These customs should not be modified with ccmake or equivalent. Instead, modify "USER_DEFINE" and "USER_FLAGS"


#Finding the locations of Trilinos and MPI by looking for files inside their directories
#This means that this cmake will work no matter what version of Trilinos is used.
find_path(Trilinos_Location NAMES bin/cth_pressure_map CACHE FILEPATH "Specify where you installed the Trilinos libraries.")
find_path(MPI_Location NAMES bin/mpicc CACHE FILEPATH "Specify where you installed the MPI libraries.")

#Converting variable names from ones which are easy for people to read to ones which follow the CMAKE standards.
set(MPI_TOP ${MPI_Location})
set(TRILINOS_TOP ${Trilinos_Location})
set(USER_DEFINE ${User_Define})
set(USER_FLAGS ${User_Flags})
set(USER_C_FLAGS ${User_C_Flags})
set(USER_CXX_FLAGS ${User_C_Flags})
set(USER_F_FLAGS ${User_C_Flags})
set(USER_LD_FLAGS ${User_C_Flags})


#Setting other libraries (can be overridden by user)
set(Sparse_Location "${TRILINOS_TOP}/../sparse" CACHE FILEPATH "Specify where you installed the Sparse libraries.")
set(Arpack_Location "${TRILINOS_TOP}/../ARPACK" CACHE FILEPATH "Specify where you installed the Arpack libraries.")
set(Library_Location "${TRILINOS_TOP}/.." CACHE FILEPATH "Specify where you want the Goma libraries to be temporarilly placed.")



set(SPARSE_LIB "-L${Sparse_Location}/lib -lsparse")
set(ARPACK_LIB "-L${Arpack_Location}/../ARPACK/lib -larpack")

#Converting variable names from ones which are easy for people to read to ones which follow the CMAKE standards.
set(LIB_DIR ${Library_Location})
set(MPI_TOP ${MPI_Location})
set(TRILINOS_TOP ${Trilinos_Location})

set(GOMA_MPI_LIB "-L${MPI_TOP}/lib -lmpi -lmpi_f77")
set(FORTRAN_LIB "-lc -lgcc -lstdc++ -lgfortran")
set(SYS_LIB "-lm -lz")





set(GCC_TOP /usr)
set(PREFIX ${CMAKE_SOURCE_DIR}/build)



set (PROJECT_HEADER_DIR "include")
set (PROJECT_SRC_DIR "src")


include_directories("${PROJECT_SOURCE_DIR}/${PROJECT_SRC_DIR}")
add_subdirectory(${PROJECT_SRC_DIR})



#Now done in SRC

#file(GLOB_RECURSE GOMA_SRC ${CMAKE_SOURCE_DIR}/${PROJECT_SRC_DIR}/*.cpp ${CMAKE_SOURCE_DIR}/${PROJECT_SRC_DIR}/*.C ${CMAKE_SOURCE_DIR}/${PROJECT_SRC_DIR}/*.F ${CMAKE_SOURCE_DIR}/${PROJECT_SRC_DIR}/*.c)
#file(GLOB_RECURSE GOMA_HEADERS ${CMAKE_SOURCE_DIR}/${PROJECT_HEADER_DIR}/*.h ${CMAKE_SOURCE_DIR}/${PROJECT_HEADER_DIR}/*.H)


#make install support
install (TARGETS "${CMAKE_SOURCE_DIR}/src/${PROJECT_NAME}" DESTINATION bin)
install (SCRIPT "${CMAKE_SOURCE_DIR}/install.cmake")


#make uninstall support

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)









#add_subdirectory(src)
#add_subdirectory(include)


#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()


