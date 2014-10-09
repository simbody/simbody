# add in local source and headers here
FILE(GLOB src_files  ./*.c)
FILE(GLOB incl_files ./*.h)

foreach( src ${src_files} )
    list(APPEND SOURCE_GROUPS "Source Files\\c-cmaes")
    list(APPEND SOURCE_GROUP_FILES ${src})
endforeach()

foreach( incl ${incl_files} )
    list(APPEND SOURCE_GROUPS "Header Files\\c-cmaes")
    list(APPEND SOURCE_GROUP_FILES ${incl})
endforeach()

set(SOURCE_GROUPS ${SOURCE_GROUPS} PARENT_SCOPE)
set(SOURCE_GROUP_FILES ${SOURCE_GROUP_FILES} PARENT_SCOPE)

# append to the local scope copy, and then copy up to parent scope
list(APPEND SOURCE_FILES ${src_files})
set(SOURCE_FILES ${SOURCE_FILES} PARENT_SCOPE)

list(APPEND SOURCE_INCLUDE_FILES ${incl_files})
set(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} PARENT_SCOPE)
