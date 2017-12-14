# visit the CPodes directory to gather up source and tests
add_subdirectory(src/CPodes)

# add in local source and headers here
file(GLOB src_files  ./src/*.cpp)
file(GLOB incl_files ./src/*.h)

# append to the local scope copy, and then copy up to parent scope
list(APPEND SOURCE_FILES ${src_files})
set(SOURCE_FILES ${SOURCE_FILES} PARENT_SCOPE)

list(APPEND SOURCE_INCLUDE_FILES ${incl_files})
set(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} PARENT_SCOPE)

# We don't have any include dirs to add here but we have to pass
# up anything added in a subdirectory.
set(SOURCE_INCLUDE_DIRS ${SOURCE_INCLUDE_DIRS} PARENT_SCOPE)
set(SOURCE_GROUPS ${SOURCE_GROUPS} PARENT_SCOPE)
set(SOURCE_GROUP_FILES ${SOURCE_GROUP_FILES} PARENT_SCOPE)
