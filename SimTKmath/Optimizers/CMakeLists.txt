# visit the following directories to gather up source
add_subdirectory(src/IpOpt)
add_subdirectory(src/c-cmaes)

# add in local source and headers here
file(GLOB src_files  ./src/*.cpp)
file(GLOB incl_files ./src/*.h)

# append to the local scope copy, and then copy up to parent scope
list(APPEND SOURCE_FILES ${src_files})
set(SOURCE_FILES ${SOURCE_FILES} PARENT_SCOPE)

list(APPEND SOURCE_INCLUDE_FILES ${incl_files})
set(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} PARENT_SCOPE)

# pass up the chain
set(SOURCE_INCLUDE_DIRS ${SOURCE_INCLUDE_DIRS} PARENT_SCOPE)
set(SOURCE_GROUPS ${SOURCE_GROUPS} PARENT_SCOPE)
set(SOURCE_GROUP_FILES ${SOURCE_GROUP_FILES} PARENT_SCOPE)
