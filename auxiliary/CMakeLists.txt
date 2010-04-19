SET(AUX_SHARED_TARGET ${SIMBODY_LIBRARY_NAME}_aux)
SET(AUX_STATIC_TARGET ${SIMBODY_LIBRARY_NAME}_aux_static)
SET(AUX_SHARED_TARGET_VN ${SIMBODY_LIBRARY_NAME}_aux${VN})
SET(AUX_STATIC_TARGET_VN ${SIMBODY_LIBRARY_NAME}_aux${VN}_static)

IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(AUX_SHARED_TARGET ${AUX_SHARED_TARGET}_d)
    SET(AUX_STATIC_TARGET ${AUX_STATIC_TARGET}_d)
    SET(AUX_SHARED_TARGET_VN ${AUX_SHARED_TARGET_VN}_d)
    SET(AUX_STATIC_TARGET_VN ${AUX_STATIC_TARGET_VN}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

## Test against the unversioned libraries if they are being built;
## otherwise against the versioned libraries.
IF(BUILD_UNVERSIONED_LIBRARIES)
	SET(TEST_AUX_SHARED_TARGET ${AUX_SHARED_TARGET})
	SET(TEST_AUX_STATIC_TARGET ${AUX_STATIC_TARGET})
ELSE(BUILD_UNVERSIONED_LIBRARIES)
	SET(TEST_AUX_SHARED_TARGET ${AUX_SHARED_TARGET_VN})
	SET(TEST_AUX_STATIC_TARGET ${AUX_STATIC_TARGET_VN})
ENDIF(BUILD_UNVERSIONED_LIBRARIES)


INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)

FILE(GLOB_RECURSE AUX_SOURCE_FILES ${PROJECT_SOURCE_DIR}/auxiliary/src/*.cpp )
FILE(GLOB_RECURSE AUX_INCLUDE_FILES ${PROJECT_SOURCE_DIR}/auxiliary/src/*.h 
${PROJECT_SOURCE_DIR}/auxiliary/include/*.h )


## Find VTK
find_path(vtk_version_include_dir
    NAMES 
        vtkRenderWindow.h
    PATHS
        "${SimTK_SDK}/include"
        "$ENV{PROGRAMFILES}/VTK/include"
        /usr/local/include
        /usr/local/vtk/include
    PATH_SUFFIXES
        vtk-5.6
        vtk-5.4
        vtk-5.2
        vtk-5.0
    DOC "Path to VTK include files, for use in finding VTK install location"
)
if(vtk_version_include_dir)
    mark_as_advanced(vtk_version_include_dir)
    get_filename_component(vtk_version_temp ${vtk_version_include_dir} NAME)
    SET(VTK_VERSION "${vtk_version_temp}" CACHE STRING "VTK version to use")
    # Go up two directory levels to find VTK_DIR
    get_filename_component(vtk_include_dir_temp "${vtk_version_include_dir}" PATH)
    get_filename_component(vtk_dir_temp "${vtk_include_dir_temp}" PATH)
    SET(VTK_DIR "${vtk_dir_temp}" CACHE PATH "VTK installation directory")
    SET(VTK_LIB_DIR "${VTK_DIR}/lib/${VTK_VERSION}" CACHE PATH "Location of VTK libraries")
    SET(VTK_INCLUDE_DIR "${VTK_DIR}/include/${VTK_VERSION}" CACHE PATH "Location of VTK libraries")
else(vtk_version_include_dir)
    mark_as_advanced(CLEAR vtk_version_include_dir)
    message("Unable to find VTK include files")
endif(vtk_version_include_dir)

INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIR})
LINK_DIRECTORIES(${VTK_LIB_DIR})

SET(VTK_LIBS_TO_USE
     vtkHybrid
     vtkFiltering
     vtkGenericFiltering
     vtkVolumeRendering
     vtkRendering
     vtkIO
     vtkGraphics
     vtkCommon
)

# explicitly add Carbon framework as one part of weaning ourselves from 
# accursed include($USE_VTK) fiasco
if(APPLE)
    find_library(CARBON_FRAMEWORK Carbon)
endif(APPLE)

ADD_SUBDIRECTORY(staticTarget)
ADD_SUBDIRECTORY(sharedTarget)

IF (BUILD_TESTING)
    ADD_SUBDIRECTORY(tests)
ENDIF (BUILD_TESTING)

IF (BUILD_EXAMPLES)
    ADD_SUBDIRECTORY(examples)
ENDIF (BUILD_EXAMPLES)

