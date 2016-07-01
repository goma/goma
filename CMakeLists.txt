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
# by this script to the repository or everything will break in the future.
#
# Also reject any pull requests attempting to do the same.
# It's just not worth debugging when everything falls apart 2 years later.
#


#CMake standards for this program:
#Editable (with ccmake etc.) variables, also known as human readable variables follow this format Word_Word_Word
#Every other variable follows this format WORD_WORD_WORD
#All arguments are in all caps (EG REQUIRED)
#Endif() needs empty parenthesis
#4 spaces of indentation for all loops/statements
#Comments should usually go before the code they explain
#All variables are called with ${} and not $()


#This tells cmake to work in the future when cmake decides to change everything.
#Also to not work with really old hardware. I didn't do any version testing
#But I am fairly certain this works with cmake 3.0
cmake_minimum_required (VERSION 3.0)
cmake_policy(VERSION 3.0)

#This tells the labelling system to behave a certain way, singlehandedly eliminating hundreds of errors.
#It makes cmake work with cmake versions newer than 3.23
cmake_policy(SET CMP0053 OLD)

#Party like it's 1977. This is in here for the two fortran files inside Goma.
enable_language(Fortran)

#This tells goma where to compile to, and the prefix of every variable
set(PROJECT_NAME goma)
include(${CMAKE_SOURCE_DIR}/CMakeScripts.cmake)
project(${PROJECT_NAME} NONE)



#Run "ccmake ." to change any options here after running "cmake ."
option(Debug_Mode "This enables the debug flag in make, building it with Debug enabled." OFF)
option(Has_Stratimikos "If you have stratimikos and want support for it enable this. This modifies the definitions" OFF)

#There is no option for the c++11 flag because Trilinos already handles it.
#I am not quite sure why but it compiles just fine without it.


#This adds a flag to the "User_Flags" section for stratimikos support. This flag is the same as the old makefile but I haven't tested it.
if(Has_Stratimikos)
    set(User_Define "${User_Define} -DHAVE_STRATIMIKOS")
endif()
if(Debug_Mode)
#Warn all. Can also do more in the future depending on what you all need. This requires overriding User_C_Flags and User_CXX_Flags
    set(User_Flags "${User_Flags} -Wall")
    set(DCMAKE_BUILD_TYPE "Debug")
endif()


#These are all set here so you can modify them with ccmake. This works because they are set with the "CACHE" flag. The last input in set is the description you will see
set(User_Define "${User_Define} -Dlinux -DCOMPILER_64BIT -DENABLE_AMESOS -DTRILINOS -DCHECK_FINITE -DNO_CHEBYSHEV_PLEASE -DMDE=27 -DMAX_PROB_VAR=15 -DMAX_EXTERNAL_FIELD=4 -DMAX_CONC=4 -DCOUPLED_FILL -DPARALLEL -DEIGEN_SERIAL" CACHE STRING "Enter or modify the compiler definitions. This overrides any definition related options.")

set(User_Flags "${User_Flags}" CACHE STRING "Enter any compiler flags you want to run on all three compilers. This overrides any cache related options.")
#-w suppresses warnings becasue gosh darn it the Trilinos people are engineers not programmers.
#And besides, there is no OSHA for code. Yet.
set(User_C_Flags "-w" CACHE STRING "Enter any extra compiler flags you want for all C compilation.")
set(User_CXX_Flags "-w" CACHE STRING "Enter any extra compiler flags you want for all C++ compilation.")
set(User_F_Flags "" CACHE STRING "Enter any extra compiler flags you want for all Fortran compilation.")
set(User_LD_Flags "" CACHE STRING "Enter any extra compiler flags you want when linking files and creating the Goma executable.")





#Finding the locations of Trilinos and MPI by looking for files inside their directories
#This means that this cmake will work no matter what version of Trilinos is used.
find_path(Trilinos_Location NAMES bin/cth_pressure_map CACHE FILEPATH "Specify where you installed the Trilinos libraries.")
find_path(MPI_Location NAMES bin/mpicc CACHE FILEPATH "Specify where you installed the MPI libraries.")

#Converting variable names from ones which are easy for people to read to ones which follow the CMAKE variable standards.
set(MPI_TOP ${MPI_Location})
set(TRILINOS_TOP ${Trilinos_Location})
set(USER_DEFINE ${User_Define})
set(USER_FLAGS ${User_Flags})
set(USER_C_FLAGS ${User_C_Flags})
set(USER_CXX_FLAGS ${User_C_Flags})
set(USER_F_FLAGS ${User_C_Flags})
set(USER_LD_FLAGS ${User_C_Flags})


#Setting other libraries (can be overridden by user). These all are set off of Trilinos by default (under the assumption you installed all your libraries in the same place).
#This is done because find_path isn't great and it's best to avoid it whenever possible to make the program less likely to fail.
set(Sparse_Location "${TRILINOS_TOP}/../sparse" CACHE FILEPATH "Specify where you installed the Sparse libraries.")
set(Arpack_Location "${TRILINOS_TOP}/../ARPACK" CACHE FILEPATH "Specify where you installed the Arpack libraries.")
set(Library_Location "${TRILINOS_TOP}/.." CACHE FILEPATH "Specify where you want the Goma libraries to be temporarilly placed.")


#Once again back to following standards. Also done so that some flags can be placed on the libraries.
#People aren't going to be seeing these variable names anyway.
set(SPARSE_LIB "-L${Sparse_Location}/lib -lsparse")
set(ARPACK_LIB "-L${Arpack_Location}/../ARPACK/lib -larpack")
set(LIB_DIR ${Library_Location})
set(MPI_TOP ${MPI_Location})
set(TRILINOS_TOP ${Trilinos_Location})

set(GOMA_MPI_LIB "-L${MPI_TOP}/lib -lmpi -lmpi_f77")
set(FORTRAN_LIB "-lc -lgcc -lstdc++ -lgfortran")
set(SYS_LIB "-lm -lz")


#This is set to match the makefile and might not be nessisary but I haven't tested it yet. TODO: Test this.
set(GCC_TOP /usr)


#Unfortunately this is overridden by Trilinos but if it wasn't then all files would go here.
set(PREFIX ${CMAKE_SOURCE_DIR}/build)
#This is set so it can be used by the CMakeLists in src.
set (PROJECT_HEADER_DIR "include")
#This tells the program where to find the Goma source. It is names as such to separate it from ${PROJECT_SOURCE_DIR}
set (GOMA_SOURCE_FILES "src")


#Runs cmake on the file in src/
include_directories("${PROJECT_SOURCE_DIR}/${GOMA_SOURCE_FILES}")
add_subdirectory(${GOMA_SOURCE_FILES})




#make install support
#This passes some data to install.cmake which uses it to add a few lines to your .profile
install (CODE "execute_process(${CMAKE_COMMAND} -DMPI_TOP=${MPI_TOP} -DTRILINOS_TOP=${TRILINOS_TOP} -P ${CMAKE_SOURCE_DIR}/install.cmake)" )



#make uninstall support
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

