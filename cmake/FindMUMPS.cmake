find_library(
  MUMPS_LIBRARIES
  NAMES dmumps mumps_common pord
  HINTS ${MUMPS_PREFIX}/MUMPS/lib ${MUMPS_PREFIX}/MUMPS
        ${MUMPS_PREFIX}/lib ${MUMPS_PREFIX} $ENV{MUMPSDIR} $ENV{MUMPSD_DIR}
        $ENV{MUMPSDIR}/lib $ENV{MUMPS_DIR}/lib ${LIB_INSTALL_DIR}
  PATH_SUFFIXES lib)

find_path(
  MUMPS_INCLUDES
  NAMES dmumps_c.h
  HINTS ${MUMPS_DIR}/include $ENV{MUMPS_DIR}/include ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES mumps)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS DEFAULT_MSG MUMPS_INCLUDES MUMPS_LIBRARIES)
mark_as_advanced(MUMPS_INCLUDES MUMPS_LIBRARIES)
