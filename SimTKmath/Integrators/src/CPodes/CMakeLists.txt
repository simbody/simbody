# visit the sundials directory to gather up source and tests
add_subdirectory(sundials)

# add in local source and headers here
file(GLOB src_files  ./*.cpp)
file(GLOB incl_files ./*.h)

# append to the local scope copy, and then copy up to parent scope
list(APPEND SOURCE_FILES ${src_files})
set(SOURCE_FILES ${SOURCE_FILES} PARENT_SCOPE)

list(APPEND SOURCE_INCLUDE_FILES ${incl_files})
set(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} PARENT_SCOPE)

list(APPEND SOURCE_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/sundials/include)
set(SOURCE_INCLUDE_DIRS ${SOURCE_INCLUDE_DIRS} PARENT_SCOPE)

# just pass these up
set(SOURCE_GROUPS ${SOURCE_GROUPS} PARENT_SCOPE)
set(SOURCE_GROUP_FILES ${SOURCE_GROUP_FILES} PARENT_SCOPE)

# define tests that depend on sources here
if(BUILD_TESTING)
    add_subdirectory(tests)
endif(BUILD_TESTING)
