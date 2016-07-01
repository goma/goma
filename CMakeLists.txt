# CMAKE file for Goma 6.1
cmake_minimum_required (VERSION 3.0)

#Party like it's 1977
enable_language(Fortran)


include(${CMAKE_SOURCE_DIR}/ProjectName.cmake)
include(${CMAKE_SOURCE_DIR}/CMakeScripts.cmake)


#It is required that you name the project in the main file in addition to an include
project(${PROJECT_NAME} NONE)


#Sets CXX11 on so you can use STD11
#This can be turned off if the person uses the graphical cmake though not recommended
option(Goma_ENABLE_CXX11_DEFAULT "Use the CXX11 flag. Required if you used the normal build script" ON)
option(HAS_STRATIMIKOS "If you have stratimikos and want support for it enable this" OFF)

if(HAS_STRATIMIKOS)
    set(CUSTOM_DEFINE "${CUSTOM_DEFINE} -DHAVE_STRATIMIKOS")
endif()



find_path(TRILINOS_TOP NAMES bin/cth_pressure_map)
find_path(MPI_TOP NAMES bin/mpicc)


set(SPARSE_LIB "-L${TRILINOS_TOP}/../sparse/lib -lsparse")
set(ARPACK_LIB "-L${TRILINOS_TOP}/../ARPACK/lib -larpack")

set(LIB_DIR "${TRILINOS_TOP}/..")

set(GOMA_MPI_LIB "-L/home/eli/goma/lib/openmpi-1.6.4/lib -lmpi -lmpi_f77")
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










#add_subdirectory(src)
#add_subdirectory(include)


#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()


