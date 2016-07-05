#The project name must be set in every CMakeLists files
#This is so you can have subprojects inside projects with a different name
project(goma)
message(STATUS "Adding Goma files to Make")

#Telling cmake where to find all the trilinos cmake files. These import all the libraries used by trilinos automatically
set(TRILINOS_PREFIX ${TRILINOS_TOP}/lib/cmake/Trilinos)
set(CMAKE_PREFIX_PATH ${TRILINOS_PREFIX} ${CMAKE_PREFIX_PATH})

#The required flag will cause cmake to fail if it cannot find Trilinos. It's better that cmake breaks first before make breaks because
#Make errors are more confusing and annoying to debug
find_package(Trilinos REQUIRED)

#Alert the user that everything is going okay.
message("\nFound Trilinos!")
message("   Trilinos_DIR = ${Trilinos_DIR}")
message("   Trilinos_VERSION = ${Trilinos_VERSION}")
message("Switching to Trilinos mode!\n")

#This uses the MPI compilers instead of the normal ones which allows for multithreaded operations
set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

#This adds the trilinos flags to the Link Directories operation
set(LD_FLAGS "-O ${Trilinos_EXTRA_LD_FLAGS}")

#These libraries are used as the link list. STARTING_LIBS and ENDING_LIBS are for consistancy with the original makefile
set(TRILINOS_LIB_COMPLETELIST "${Trilinos_LIBRARY_DIRS} ${Trilinos_LIBRARIES}")
set(STARTING_LIBS "${SPARSE_LIB} ${TRILINOS_LIB_COMPLETELIST} ${ARPACK_LIB}")

#These are called libs, and kept in a separate variable, to keep consistancy with the original makefile
set(ENDING_LIBS "${FORTRAN_LIB} ${SYS_LIB}")


#When debug mode is off, -w is run on C++ and C to suppress all the trilinos warnings that aren't Goma's responsibility.
#There are a lot of warnings and they can confuse the average user

#CMAKE Flags can be specified with CMAKE -DCMAKE_C_FLAGS or by editing the user flags.
set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS} ${USER_FLAGS} ${USER_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${USER_FLAGS} ${USER_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} ${USER_FLAGS} ${USER_F_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${LD_FLAGS} ${USER_FLAGS} ${USER_LD_FLAGS}")

#Imports all nessisary headers.
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${MPI_TOP}/include)
include_directories(${TRILINOS_TOP}/include)

#Adds all definitions in the old Goma makefile. The equivalent of Make -D${USER_DEFINE} But this variable can also be customized to add any definitions specified by the user.
add_definitions(${USER_DEFINE})



get_source_files(SOURCE_FILES)
include_directories(${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})



add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${SPARSE_LIB} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${ARPACK_LIB} ${GOMA_MPI_LIB} ${Trilinos_EXTRA_LD_FLAGS} ${ENDING_LIBS})


#Install goma to $CMAKE_INSTALL_PREFIX/bin (by default /usr/local/bin). This can only be specified with -DCMAKE_INSTALL_PREFIX="LOCATION"
install(TARGETS goma DESTINATION bin)


#make uninstall support done at the CMake above this one.
