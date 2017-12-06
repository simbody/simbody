include(FindDoxygen)

if(DOXYGEN_EXECUTABLE-NOTFOUND)
else(DOXYGEN_EXECUTABLE-NOTFOUND)
    set(DOXY_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/Doxyfile")

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in 
          ${DOXY_CONFIG}
          @ONLY )

    add_custom_target(doxygen ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG}) 

    file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/html/")
    install(DIRECTORY "${PROJECT_BINARY_DIR}/html/"
            DESTINATION "${SIMBODY_INSTALL_DOXYGENDIR}"
            )
    # This is just a shortcut to the Doxygen index.html.
    install(FILES "SimbodyAPI.html" DESTINATION "${CMAKE_INSTALL_DOCDIR}")
endif(DOXYGEN_EXECUTABLE-NOTFOUND)

