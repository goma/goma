###/************************************************************************ *
###* Goma - Multiphysics finite element software                             *
###* Sandia National Laboratories                                            *
###*                                                                         *
###* Copyright (c) 2016 Sandia Corporation.                                  *
###*                                                                         *
###* Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,  *
###* the U.S. Government retains certain rights in this software.            *
###*                                                                         *
###* This software is distributed under the GNU General Public License.      *
###\************************************************************************/
###
### CMake file for Goma 6.1
###########################
###
### If you are a developer on this project: Do not commit any files generated
### by this script to the repository or everything will break in the future.
###
### Also reject any pull requests attempting to do the same.
### It's just not worth debugging when everything falls apart 2 years later.
###


### CMake standards for this program:
### Most variables follows this format ${PROJECT_NAME}_FOO_DATA
### All arguments are in all caps (EG REQUIRED)
### Endif() needs empty parenthesis
### Comments should usually go before the code they explain


### This tells cmake to work in the future when cmake updates.
### 2.8.4 required for fully functional make debug
cmake_minimum_required (VERSION 2.8)
cmake_policy(VERSION 2.8)


### Party like it's 1977. This is in here for the two fortran files inside Goma.
enable_language(Fortran)

### This tells goma where to compile to, and the prefix of every variable
set(PROJECT_NAME goma)
include(${CMAKE_SOURCE_DIR}/cmake/CMakeScripts.cmake)

### Run "ccmake ." to change any options here after running "cmake ."
option(${PROJECT_NAME}_CXX11_STD "Unless you have an older compiler, or you used the NoCXX11 library, keep this enabled" ON)
option(${PROJECT_NAME}_Write_To_Profile "This will automatically add the libraries to your profile, it is recommended to keep this off" OFF)
option(${PROJECT_NAME}_Create_Run_Script "This will create rungoma, a script which runs goma with all nessisary settings" ON)


if(${PROJECT_NAME}_CXX11_STD)
  set(${PROJECT_NAME}_EXTRA_CXX_FLAGS "${${PROJECT_NAME}_EXTRA_CXX_FLAGS} -std=c++11")
endif()

### If you ever want to add additional flags when running debug (other than -g which is automatically added)
if (${CMAKE_BUILD_TYPE} MATCHES DEBUG)
  set(${${PROJECT_NAME}_EXTRA_C_FLAGS} "${${PROJECT_NAME}_EXTRA_C_FLAGS} -Wall")
  set(${${PROJECT_NAME}_EXTRA_CXX_FLAGS} "${${PROJECT_NAME}_EXTRA_CXX_FLAGS} -Wall")
  set(${${PROJECT_NAME}_EXTRA_F_FLAGS} "${${PROJECT_NAME}_EXTRA_F_FLAGS}")

else()
  set(${${PROJECT_NAME}_EXTRA_C_FLAGS} "${${PROJECT_NAME}_EXTRA_C_FLAGS} -Wall")
  set(${${PROJECT_NAME}_EXTRA_CXX_FLAGS} "${${PROJECT_NAME}_EXTRA_CXX_FLAGS} -Wall")
  set(${${PROJECT_NAME}_EXTRA_F_FLAGS} "${${PROJECT_NAME}_EXTRA_F_FLAGS} ")
endif()

### These are all set here so you can modify them with ccmake. This works because they are set with the "CACHE" flag. The last input in set is the description you will see
set(${PROJECT_NAME}_User_Define "-Dlinux"
  "-DCOMPILER_64BIT"
  "-DENABLE_AMESOS"
  "-DTRILINOS"
  "-DCHECK_FINITE"
  "-DNO_CHEBYSHEV_PLEASE"
  "-DMDE=27"
  "-DMAX_PROB_VAR=15"
  "-DMAX_EXTERNAL_FIELD=4"
  "-DMAX_CONC=4"
  "-DCOUPLED_FILL"
  "-DPARALLEL"
  "-DEIGEN_SERIAL" CACHE STRING "Enter or modify the compiler definitions. This is in addition to any options above.")

set(${PROJECT_NAME}_User_Flags "" CACHE STRING "Enter any compiler flags you want to run on all three compilers. This is in addition to any options above.")
set(${PROJECT_NAME}_User_C_Flags "" CACHE STRING "Enter any extra compiler flags you want for all C compilation.")
set(${PROJECT_NAME}_User_CXX_Flags "" CACHE STRING "Enter any extra compiler flags you want for all C++ compilation.")
set(${PROJECT_NAME}_User_F_Flags "" CACHE STRING "Enter any extra compiler flags you want for all Fortran compilation.")
set(${PROJECT_NAME}_User_LD_Flags "" CACHE STRING "Enter any extra compiler flags you want when linking files and creating the Goma executable.")


if("${${PROJECT_NAME}_LIBS}" STREQUAL "")
  find_path(${PROJECT_NAME}_Trilinos_DIR NAMES lib/libamesos.a CACHE FILEPATH "Specify where you installed the Trilinos libraries.")
  find_path(${PROJECT_NAME}_MPI_DIR NAMES bin/mpicc CACHE FILEPATH "Specify where you installed the MPI libraries.")
  set(${PROJECT_NAME}_LIBS "${${PROJECT_NAME}_Trilinos_DIR}/.." CACHE FILEPATH "Specify where you want the Goma libraries to be temporarilly placed.")
  set(${PROJECT_NAME}_Sparse_DIR "${${PROJECT_NAME}_LIBS}/sparse" CACHE FILEPATH "Specify where you installed the Sparse libraries.")
  set(${PROJECT_NAME}_Arpack_DIR "${${PROJECT_NAME}_LIBS}/ARPACK" CACHE FILEPATH "Specify where you installed the Arpack libraries.")
  set(${PROJECT_NAME}_SEACAS_DIR "${${PROJECT_NAME}_LIBS}/seacas" CACHE FILEPATH "If SEACAS is standalone, specify where it is installed.")
else()
  find_path(${PROJECT_NAME}_Trilinos_DIR NAMES lib/libamesos.a PATHS ${${PROJECT_NAME}_LIBS}/* CACHE FILEPATH "Specify where you installed the Trilinos libraries.")
  find_path(${PROJECT_NAME}_MPI_DIR NAMES bin/mpicc PATHS ${${PROJECT_NAME}_LIBS}/* CACHE FILEPATH "Specify where you installed the MPI libraries.")
  set(${PROJECT_NAME}_Sparse_DIR "${${PROJECT_NAME}_LIBS}/sparse" CACHE FILEPATH "Specify where you installed the Sparse libraries.")
  set(${PROJECT_NAME}_Arpack_DIR "${${PROJECT_NAME}_LIBS}/ARPACK" CACHE FILEPATH "Specify where you installed the Arpack libraries.")
  set(${PROJECT_NAME}_SEACAS_DIR "${${PROJECT_NAME}_LIBS}/seacas" CACHE FILEPATH "If SEACAS is standalone, specify where it is installed.")
endif()




### Finding the locations of Trilinos and MPI by looking for files inside their directories
### This means that this cmake will work no matter what version of Trilinos is used.
find_path(${PROJECT_NAME}_Trilinos_DIR NAMES bin/cth_pressure_map CACHE FILEPATH "Specify where you installed the Trilinos libraries.")
find_path(${PROJECT_NAME}_MPI_DIR NAMES bin/mpicc CACHE FILEPATH "Specify where you installed the MPI libraries.")


### Converting variable names from ones which are easy for people to read to ones which follow the CMAKE variable standards.
### Also adding the flags defined by the variables
set(${PROJECT_NAME}_MPI_TOP "${${PROJECT_NAME}_MPI_DIR}")
set(${PROJECT_NAME}_TRILINOS_TOP "${${PROJECT_NAME}_Trilinos_DIR}")
set(${PROJECT_NAME}_USER_DEFINE "${${PROJECT_NAME}_User_Define} ${${PROJECT_NAME}_EXTRA_DEFINE}")
set(${PROJECT_NAME}_USER_FLAGS "${${PROJECT_NAME}_User_Flags} ${${PROJECT_NAME}_EXTRA_FLAGS}")
set(${PROJECT_NAME}_USER_C_FLAGS "${${PROJECT_NAME}_User_C_Flags} ${${PROJECT_NAME}_EXTRA_C_FLAGS}")
set(${PROJECT_NAME}_USER_CXX_FLAGS "${${PROJECT_NAME}_User_CXX_Flags} ${${PROJECT_NAME}_EXTRA_CXX_FLAGS}")
set(${PROJECT_NAME}_USER_F_FLAGS "${${PROJECT_NAME}_User_F_Flags}")
set(${PROJECT_NAME}_USER_LD_FLAGS "${${PROJECT_NAME}_User_LD_Flags}")



### Setting other libraries (can be overridden by user). These all are set off of Trilinos by default (under the assumption you installed all your libraries in the same place).
### This is done because find_path isn't great and it's best to avoid it whenever possible to make the program less likely to fail.

set(${PROJECT_NAME}_Sparse_DIR "${${PROJECT_NAME}_Trilinos_DIR}/../sparse" CACHE FILEPATH "Specify where you installed the Sparse libraries.")
set(${PROJECT_NAME}_Arpack_DIR "${${PROJECT_NAME}_Trilinos_DIR}/../ARPACK" CACHE FILEPATH "Specify where you installed the Arpack libraries.")
set(${PROJECT_NAME}_SEACAS_DIR "${${PROJECT_NAME}_Trilinos_DIR}/../seacas" CACHE FILEPATH "If SEACAS is standalone, specify where it is installed.")
set(${PROJECT_NAME}_LIBS "${${PROJECT_NAME}_Trilinos_DIR}/.." CACHE FILEPATH "Specify where you want the Goma libraries to be temporarilly placed.")


### Once again back to following standards. Also done so that some flags can be placed on the libraries.
### People aren't going to be seeing these variable names anyway.
set(${PROJECT_NAME}_SPARSE_LIB "-L${${PROJECT_NAME}_Sparse_DIR}/lib -lsparse")
set(${PROJECT_NAME}_ARPACK_LIB "-L${${PROJECT_NAME}_Arpack_DIR}/lib -larpack")
set(${PROJECT_NAME}_LIB_DIR ${${PROJECT_NAME}_LIBS})
set(${PROJECT_NAME}_MPI_LIB "-L${${PROJECT_NAME}_MPI_TOP}/lib -lmpi -lmpi_f77")
set(${PROJECT_NAME}_FORTRAN_LIB "-lc -lgcc -lstdc++ -lgfortran")
set(${PROJECT_NAME}_SYS_LIB "-lm -lz")


### I really wish I could combine these two together but you can't. (set (Var EXISTS) doesn't work)
### This will look for a standalone seacas directory.
if(IS_DIRECTORY "${${PROJECT_NAME}_SEACAS_DIR}/")
  set(${PROJECT_NAME}_FOUND_SEACAS ON)
else()
  set(${PROJECT_NAME}_FOUND_SEACAS OFF)
endif()

option(${PROJECT_NAME}_Standalone_SEACAS "This will build Goma using standalone SEACAS libraries, instead of the ones build with Trilinos." ${${PROJECT_NAME}_FOUND_SEACAS})

if(${${PROJECT_NAME}_Standalone_SEACAS})
  set(${PROJECT_NAME}_SEACAS_LIB "-L${${PROJECT_NAME}_SEACAS_DIR}/lib -lexoIIv2c -L${${PROJECT_NAME}_SEACAS_DIR}/lib -lnetcdf -lhdf5_hl -lhdf5 -lchaco")
else()
  set(${PROJECT_NAME}_SEACAS_LIB "")
endif()


### This is set so it can be used by the CMakeLists in src.
set (${PROJECT_NAME}_HEADER_DIR "include")
### This tells the program where to find the Goma source. It is names as such to separate it from ${${PROJECT_NAME}_PROJECT_SOURCE_DIR}
set (${PROJECT_NAME}_SOURCE_FILES "src")




message(STATUS "Adding ${PROJECT_NAME} files to Make")

### Telling cmake where to find all the trilinos cmake files. These import all the libraries used by trilinos automatically
set(Trilinos_PREFIX ${${PROJECT_NAME}_TRILINOS_TOP}/lib/cmake/Trilinos)
set(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})

### make install support
### This passes some data to install.cmake which uses it to add a few lines to your .profile
if (${${PROJECT_NAME}_Write_To_Profile})
  if (${${PROJECT_NAME}_Standalone_SEACAS})
    install (CODE "execute_process(COMMAND ${CMAKE_COMMAND} -DTRILINOS_TOP=${${PROJECT_NAME}_SEACAS_DIR} -DMPI_TOP=${${PROJECT_NAME}_MPI_TOP} -P ${CMAKE_SOURCE_DIR}/cmake/install.cmake)" )
  else()
    install (CODE "execute_process(COMMAND ${CMAKE_COMMAND} -DMPI_TOP=${${PROJECT_NAME}_MPI_TOP} -DTRILINOS_TOP=${${PROJECT_NAME}_TRILINOS_TOP} -P ${CMAKE_SOURCE_DIR}/cmake/install.cmake)" )
  endif()
endif()

if (${${PROJECT_NAME}_Create_Run_Script})
  if (${${PROJECT_NAME}_Standalone_SEACAS})
    if(${CMAKE_BUILD_TYPE} MATCHES DEBUG)
      install (CODE "execute_process(COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/debugscript ${CMAKE_INSTALL_PREFIX} ${${PROJECT_NAME}_MPI_TOP} ${${PROJECT_NAME}_SEACAS_DIR})" )
    else()
      install (CODE "execute_process(COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/runscript ${CMAKE_INSTALL_PREFIX} ${${PROJECT_NAME}_MPI_TOP} ${${PROJECT_NAME}_SEACAS_DIR})" )
    endif()
  else()
    if(${CMAKE_BUILD_TYPE} MATCHES DEBUG)
      install (CODE "execute_process(COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/debugscript ${CMAKE_INSTALL_PREFIX} ${${PROJECT_NAME}_MPI_TOP} ${${PROJECT_NAME}_TRILINOS_TOP})" )
    else()
      install (CODE "execute_process(COMMAND /bin/bash ${CMAKE_SOURCE_DIR}/cmake/runscript ${CMAKE_INSTALL_PREFIX} ${${PROJECT_NAME}_MPI_TOP} ${${PROJECT_NAME}_TRILINOS_TOP})" )
    endif()
  endif()
endif()

### The required flag will cause cmake to fail if it cannot find Trilinos. It's better that cmake breaks first before make breaks because
### Make errors are more confusing and annoying to debug
find_package(Trilinos)

### ID10-t proofing the libraries.
if(NOT ${Trilinos_FOUND})
  ### I have no faith in humanity if they cannot figure this out
  message(FATAL_ERROR "


BEGIN ERROR MESSAGE:
-------------------------------------------------------
Unable to locate Trilinos, please add your libraries to cmake-config.
Or manually specify the trilinos path.


The program will now exit.")
endif()

### Alert the user that everything is going okay.
message("\nFound Trilinos version: ${Trilinos_VERSION}\n")


option(${PROJECT_NAME}_Build_Stratimikos "If you have stratimikos and want support for it enable this. This modifies the definitions" ${Trilinos_Stratimikos_FOUND})
if(${PROJECT_NAME}_Build_Stratimikos)
  set(${PROJECT_NAME}_USER_DEFINE "${${PROJECT_NAME}_USER_DEFINE} -DHAVE_STRATIMIKOS")
endif()

### This uses the MPI compilers instead of the normal ones which allows for multithreaded operations
set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

### This adds the trilinos flags to the Link Directories operation
set(LD_FLAGS "${Trilinos_EXTRA_LD_FLAGS}")


### These are the last called libs, and kept in a separate variable, to keep consistancy with the original makefile
set(${PROJECT_NAME}_ENDING_LIBS "${${PROJECT_NAME}_FORTRAN_LIB} ${${PROJECT_NAME}_SYS_LIB}")


### CMAKE Flags can be specified with CMAKE -DCMAKE_C_FLAGS or by editing the user flags.
set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} ${${PROJECT_NAME}_USER_FLAGS} ${${PROJECT_NAME}_USER_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${${PROJECT_NAME}_USER_FLAGS} ${${PROJECT_NAME}_USER_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} ${${PROJECT_NAME}_USER_FLAGS} ${${PROJECT_NAME}_USER_F_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${${PROJECT_NAME}_LD_FLAGS} ${${PROJECT_NAME}_USER_FLAGS} ${${PROJECT_NAME}_USER_LD_FLAGS}")

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${CMAKE_Fortran_FLAGS}")


### Imports all nessisary headers.
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(SYSTEM ${${PROJECT_NAME}_MPI_TOP}/include)
include_directories(SYSTEM ${${PROJECT_NAME}_TRILINOS_TOP}/include)
if(${${PROJECT_NAME}_Standalone_SEACAS})
  include_directories(SYSTEM "${${PROJECT_NAME}_SEACAS_DIR}/include")
endif()

### Adds all definitions in the old Goma makefile. The equivalent of Make -D${USER_DEFINE} But this variable can also be customized to add any definitions specified by the user.
add_definitions(${${PROJECT_NAME}_USER_DEFINE})

include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

file(GLOB SOURCE_FILES
  "src/*.cpp"
  "src/*.c"
  "src/*.F"
  "src/*.C"
  "src/brkfix/*.c"
  "include/*.h"
  "include/brkfix/*.h")


if(${CMAKE_BUILD_TYPE} MATCHES DEBUG)
  add_executable("${PROJECT_NAME}d" ${SOURCE_FILES})
else()
  add_executable(${PROJECT_NAME} ${SOURCE_FILES})
endif()

if(${CMAKE_BUILD_TYPE} MATCHES DEBUG)
  target_link_libraries("${PROJECT_NAME}d" ${${PROJECT_NAME}_SPARSE_LIB} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${${PROJECT_NAME}_ARPACK_LIB} ${${PROJECT_NAME}_MPI_LIB} ${Trilinos_EXTRA_LD_FLAGS} ${${PROJECT_NAME}_SEACAS_LIB} ${${PROJECT_NAME}_ENDING_LIBS})
else()
  target_link_libraries(${PROJECT_NAME} ${${PROJECT_NAME}_SPARSE_LIB} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${${PROJECT_NAME}_ARPACK_LIB} ${${PROJECT_NAME}_MPI_LIB} ${Trilinos_EXTRA_LD_FLAGS} ${${PROJECT_NAME}_SEACAS_LIB} ${${PROJECT_NAME}_ENDING_LIBS})
endif()

### Install goma to $CMAKE_INSTALL_PREFIX/bin (by default /usr/local/bin). This can only be specified with -DCMAKE_INSTALL_PREFIX="LOCATION"
if(${CMAKE_BUILD_TYPE} MATCHES DEBUG)
  install(TARGETS "${PROJECT_NAME}d" DESTINATION bin)
else()
  install(TARGETS ${PROJECT_NAME} DESTINATION bin)
endif()

### make uninstall support
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake)

### make debug will reconfigure with debugging enabled and then build in the library set by make config
add_custom_target(debug COMMAND "/bin/bash" ${CMAKE_SOURCE_DIR}/cmake-debug-example "-D${PROJECT_NAME}_LIBS\\=${${PROJECT_NAME}_LIB_DIR}" AUTOMAKE $(MAKE))
