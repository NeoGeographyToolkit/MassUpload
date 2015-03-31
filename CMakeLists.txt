
#===============================================================================
# Project information


project( RegisterHrsc )
cmake_minimum_required(VERSION 2.8)

# Add path for custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}" )



#===============================================================================
# Configuration options

# --- Directory paths ---
set(BASESYSTEM_INSTALL_DIR "" CACHE FILEPATH "Path to AutoTools installation directory")
if(NOT BASESYSTEM_INSTALL_DIR)
  message(SEND_ERROR "You need to set BASESYSTEM_INSTALL_DIR: ${BASESYSTEM_INSTALL_DIR}")
endif()

set(VISIONWORKBENCH_INSTALL_DIR "" CACHE FILEPATH "Path to Vision Workbench installation directory")
if(NOT VISIONWORKBENCH_INSTALL_DIR)
  message(SEND_ERROR "You need to set VISIONWORKBENCH_INSTALL_DIR: ${VISIONWORKBENCH_INSTALL_DIR}")
endif()


message("Base System path      = ${BASESYSTEM_INSTALL_DIR}")
message("Vision Workbench path = ${VISIONWORKBENCH_INSTALL_DIR}")

# Copy the nice input variable names to the messy internal names
set(VISIONWORKBENCH_ROOT   ${VISIONWORKBENCH_INSTALL_DIR})

# --- Custom options ---
option(BUILD_SHARED_LIBS "Produce shared libraries." TRUE)

# --- Fixed options ---
set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON )
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})



#==================================================================================
# Find dependencies

message("--------------------------------------------------------------")
message("Searching for boost...")

set(BOOST_ROOT "${BASESYSTEM_INSTALL_DIR}/include/boost-1_55")
set(BOOST_LIBRARYDIR  "${BASESYSTEM_INSTALL_DIR}/lib")
set(BOOST_INCLUDEDIR "${BASESYSTEM_INSTALL_DIR}/include/boost-1_55/boost")
set(Boost_NO_SYSTEM_PATHS ON) #This does not seem to work
#set(Boost_DEBUG 1) # Uncomment for useful finding boost debug output
set(BOOST_MIN_VERSION "1.55.0")
set(REQUIRED_BOOST_LIBS program_options system filesystem regex date_time thread iostreams)

find_package(Boost ${BOOST_MIN_VERSION} REQUIRED COMPONENTS ${REQUIRED_BOOST_LIBS})

if(Boost_FOUND)
    message("Successfully found Boost version ${Boost_VERSION}!")
    message("---Boost Include Dir = ${Boost_INCLUDE_DIR}")
    #set(Boost_LIBRARY_DIR ${BASESYSTEM_INSTALL_DIR}/lib) # Override this value
    message("---Boost Library Dir = ${Boost_LIBRARY_DIRS}")
    message("---Boost libraries = ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})
endif(Boost_FOUND)
message("--------------------------------------------------------------")



find_package(VisionWorkbench REQUIRED COMPONENTS  core math image camera fileio
                                                  cartography geometry interestPoint
                                                  stereo mosaic )

message("Found vision workbench include directory: ${VISIONWORKBENCH_INCLUDE_DIRS}")


find_package( OpenCV REQUIRED )

#message("dir = ${OpenCV_INCLUDE_DIRS}")




#set(PROTOBUF_FOUND)
#set(PROTOBUF_INCLUDE_DIR ${BASESYSTEM_INSTALL_DIR}/include   )
#set(PROTOBUF_LIBRARY     ${BASESYSTEM_INSTALL_DIR}/lib/libprotobuf.so)
#include_directories("${PROTOBUF_INCLUDE_DIR}")
#include_directories("${QT_INCLUDE_DIRS}")
##include_directories("${NAIF_INCLUDE_DIRS}")

include_directories(${VISIONWORKBENCH_INCLUDE_DIRS})
#include_directories(${StereoPipeline_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
#include_directories(${PROTOBUF_INCLUDE_DIR})
include_directories( ${OpenCV_INCLUDE_DIRS} )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

link_directories(${Boost_LIBRARY_DIRS})



#==================================================================================
# Executables

add_executable( RegisterHrsc RegisterHrsc.cpp )
target_link_libraries( RegisterHrsc ${OpenCV_LIBS} ${VISIONWORKBENCH_LIBRARIES})

add_executable( writeHrscColorPairs writeHrscColorPairs.cpp )
target_link_libraries( writeHrscColorPairs ${OpenCV_LIBS} ${Boost_LIBRARIES})

add_executable( transformHrscImageColor transformHrscImageColor.cpp )
target_link_libraries( transformHrscImageColor ${OpenCV_LIBS} ${Boost_LIBRARIES})

add_executable( hrscMosaic hrscMosaic.cpp )
target_link_libraries( hrscMosaic ${OpenCV_LIBS} ${Boost_LIBRARIES})

add_executable( computeBrightnessCorrection computeBrightnessCorrection.cpp )
target_link_libraries( computeBrightnessCorrection ${OpenCV_LIBS} ${Boost_LIBRARIES})



