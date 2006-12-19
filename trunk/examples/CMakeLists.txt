
FILE( GLOB SOURCES ${PROJECT_SOURCE_DIR}/examples/*.cpp )

FOREACH (src_file ${SOURCES} )
   GET_FILENAME_COMPONENT( file_name ${src_file} NAME )
   INSTALL_FILES( "/core/examples/simmath/" FILES ${file_name} )
ENDFOREACH (src_file ${SOURCES})

INSTALL_FILES( "/core/examples/simmath/" FILES ${PROJECT_SOURCE_DIR}/examples/README.txt)
