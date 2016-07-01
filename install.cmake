set (IMPORT_NAME "importgomalibraries.sh")
set (PROJECT_NAME "goma")

file (COPY "${CMAKE_SOURCE_DIR}/src/${PROJECT_NAME}" DESTINATION "/usr/local/bin/")
file (COPY "${CMAKE_SOURCE_DIR}/src/${IMPORT_NAME}" DESTINATION "/etc/profile.d")
file (WRITE "${CMAKE_SOURCE_DIR}/install_manifest.txt" "/usr/local/bin/${PROJECT_NAME}
/etc/profile.d/${IMPORT_NAME}")
