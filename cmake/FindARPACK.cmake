find_library(ARPACK_LIB NAMES arpack HINTS ${ARPACK_PREFIX}/lib ${ARPACK_PREFIX})

set(ARPACK_LIBRARIES "${ARPACK_LIB}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK DEFAULT_MSG ARPACK_LIBRARIES)
