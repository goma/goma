find_path(Sparse_INCLUDES spMatrix.h
          PATHS ${Sparse_PREFIX}/sparse $ENV{SPARSEDIR} ${Sparse_PREFIX}
          PATH_SUFFIXES include src)
find_library(
  Sparse_LIBRARIES
  NAMES sparse
  PATHS ${Sparse_PREFIX}/sparse ${Sparse_PREFIX} $ENV{SPARSEDIR} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sparse DEFAULT_MSG Sparse_LIBRARIES
                                  Sparse_INCLUDES)
