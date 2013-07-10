cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

include (FindPkgConfig)
if (PKG_CONFIG_FOUND)
  pkg_check_modules(GAZEBO gazebo)
endif()

include(FindBoost)
find_package(Boost ${MIN_BOOST_VERSION} REQUIRED program_options system filesystem regex)

find_package(Protobuf REQUIRED)

include_directories(${GAZEBO_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
link_directories(${GAZEBO_LIBRARY_DIRS} ${Boost_LIBRARY_DIRS})

add_executable(listener listener.cc)
add_executable(fakeserver server.hh server.cc)
add_executable(publisher publisher.cc)

message("Gazebo message: ${GAZEBO_LIBRARIES}")

target_link_libraries(listener ${GAZEBO_LIBRARIES} ${Boost_LIBRARIES} ${PROTOBUF_LIBRARIES})
target_link_libraries(fakeserver libgazebo.so ${GAZEBO_LIBRARIES} ${Boost_LIBRARIES} ${PROTOBUF_LIBRARIES})
target_link_libraries(publisher ${GAZEBO_LIBRARIES} ${Boost_LIBRARIES} ${PROTOBUF_LIBRARIES})
