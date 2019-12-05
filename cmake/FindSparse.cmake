find_path(Sparse_INCLUDES spMatrix.h PATHS ${Sparse_PREFIX}/include ${Sparse_PREFIX}/src)
find_library(Sparse_LIB NAMES sparse PATHS ${Sparse_PREFIX}/lib)

set(Sparse_LIBRARIES "${Sparse_LIB}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sparse DEFAULT_MSG Sparse_LIBRARIES Sparse_INCLUDES)
