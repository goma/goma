project(Goma)
message(STATUS "Adding Goma files to Make")

set(TARGETDIR bin)

set(TRILINOS_PREFIX ${TRILINOS_TOP}/lib/cmake/Trilinos)

set(CMAKE_PREFIX_PATH ${TRILINOS_PREFIX} ${CMAKE_PREFIX_PATH})
find_package(Trilinos REQUIRED)

MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("End of Trilinos details\n")
MESSAGE("Switching to Trilinos mode!\n")

set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )


set(STARTING_LIBS "${SPARSE_LIB} ${Trilinos_LIBRARY_DIRS} ${Trilinos_LIBRARIES} ${ARPACK_LIB}")



#These are NOT libraries but rather compiler commands
#They are called libs to keep consistancy with the original makefile
set(ENDING_LIBS "${GOMA_MPI_LIB} ${FORTRAN_LIB} ${SYS_LIB}")

#Suppress warnings because Goma isn't responsible for (most of) the spagetti code in trilinos and I don't want people thinking otherwise
set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS} ${STARTING_LIBS} ${ENDING_LIBS} ${Trilinos_EXTRA_LD_FLAGS} -w")
set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS} -w")
set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")





include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${MPI_TOP}/include)
include_directories(${TRILINOS_TOP}/include)
include_directories(${UMFPACK_TOP}) #UMF please




#These defs are in the old makefile and used to make Goma properly compile
add_definitions(-Dlinux
-DCOMPILER_64BIT
-DENABLE_AMESOS
-DTRILINOS
-DCHECK_FINITE
-DNO_CHEBYSHEV_PLEASE
-DMDE=27
-DMAX_PROB_VAR=15
-DMAX_EXTERNAL_FIELD=4
-DMAX_CONC=4
-DCOUPLED_FILL
-DPARALLEL
-DEIGEN_SERIAL)



get_source_files(SOURCE_FILES)



INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}) 

set_project_info(${PROJECT_NAME} "${CMAKE_SOURCE_DIR}/${TARGETDIR}")



#A quick way to debug these scripts, uncomment the 4 lines below:

get_cmake_property(_variableNames VARIABLES)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
