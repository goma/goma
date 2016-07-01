set (IMPORT_NAME "importgomalibraries.sh")

file (COPY "${CMAKE_SOURCE_DIR}/src/${IMPORT_NAME}" DESTINATION "/etc/profile.d")

