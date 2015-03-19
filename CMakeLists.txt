
cmake_minimum_required(VERSION 2.8)

project( RegisterHrsc )

find_package( OpenCV REQUIRED )

message("dir = ${OpenCV_INCLUDE_DIRS}")
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable( RegisterHrsc RegisterHrsc.cpp )

target_link_libraries( RegisterHrsc ${OpenCV_LIBS} )

