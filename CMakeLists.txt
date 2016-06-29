# CMAKE file for Goma 6.1
cmake_minimum_required (VERSION 3.0)
project(Goma)

# Goma Version
set (Goma_VERSION_MAJOR 6)
set (Goma_VERSION_MINOR 1)


set (PROJECT_HEADER_DIR "include")
set (PROJECT_SRC_DIR "src")

#TODO: Add a file to Goma src that can take advantage of the settings set by CMake
#configure_file ()



# Where are the library files you preinstalled.
set(DIR_LIB "lib" CACHE FILEPATH "Please tell me where you installed the libraries")

option(ADD_BASHRC "Use .bash_rc, and not update-alternatives, to add trilinos to your PATH?" OFF)



# Add the src directory and its cmake file
#add_subdirectory (src)
