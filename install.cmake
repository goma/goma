set (IMPORT_NAME "importgomalibraries.sh")

install (TARGETS "${CMAKE_SOURCE_DIR}/src/${PROJECT_NAME}" DESTINATION bin)
install (TARGETS "${CMAKE_SOURCE_DIR}/src/${IMPORT_NAME}" DESTINATION "/etc/profile.d")

