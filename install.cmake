set (IMPORT_NAME "importgomalibraries.sh")

file (COPY "${CMAKE_SOURCE_DIR}/src/${IMPORT_NAME}" DESTINATION "/etc/profile.d")
file (APPEND "${CMAKE_SOURCE_DIR}/install_manifest.txt" "
/etc/profile.d/${IMPORT_NAME}")
