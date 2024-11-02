
# -----------------------------------------------------------------------------
# Project definition
# -----------------------------------------------------------------------------
PROJECT(wind_mesh_map)
cmake_minimum_required (VERSION 3.1)
cmake_policy(SET CMP0057 NEW)
# -----------------------------------------------------------------------------
# CPP11
# -----------------------------------------------------------------------------
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED TRUE)
# -----------------------------------------------------------------------------
# DGtal dependencies
# -----------------------------------------------------------------------------
FIND_PACKAGE(DGtal 1.3 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
message(STATUS "DGtal found.")
# -----------------------------------------------------------------------------
# Thread dependencies
# -----------------------------------------------------------------------------
FIND_PACKAGE(Threads)
# -----------------------------------------------------------------------------
# CLI11 et CSV
# -----------------------------------------------------------------------------
include_directories( "${PROJECT_SOURCE_DIR}/ext/" )
# -----------------------------------------------------------------------------
# LIBPOINTMATCHER
# -----------------------------------------------------------------------------
find_package(libpointmatcher REQUIRED)
# -----------------------------------------------------------------------------
# TENSORFLOW_LIB
# -----------------------------------------------------------------------------
#find_library(TENSORFLOW_LIB tensorflow HINT $ENV{HOME}/libtensorflow2/lib)

#include_directories(${CMAKE_CURRENT_BINARY_DIR} ${libpointmatcher_INCLUDE_DIRS})

FIND_PACKAGE(GSL REQUIRED)
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
LINK_DIRECTORIES(${GSL_LIBRARIES})

#add_subdirectory(/volWork/these/source/TLDDC/ext/matplotplusplus/)

#find_package(PCL 1.3 REQUIRED COMPONENTS common segmentation)
#include_directories(${PCL_INCLUDE_DIRS})
#link_directories(${PCL_LIBRARY_DIRS})
#add_definitions(${PCL_DEFINITIONS})

#find_package( OpenCV REQUIRED )#core imgproc highgui imgcodecs
#if ( OpenCV_FOUND )
#  message(STATUS "OpenCV and boost_program_options found.")
#endif( OpenCV_FOUND)

#SET(PCLLib ${PCL_SEGMENTATION_LIBRARIES})


#ADD_EXECUTABLE(b_box Main_bbox ../Common/IOHelper ../Common/CylindricalCoordinateSystem)
#TARGET_LINK_LIBRARIES(b_box ${DGTAL_LIBRARIES} ${PCLLib})

#ADD_EXECUTABLE(segunroll_mesh Main_Unroll_MESH segmentation/SegmentationHelper
#    ../Common/IOHelper ../Common/IOHelper  ../Common/Statistic ../Common/CylindricalPoint ../Common/Rosin ../Common/Morphology
#    ../Common/CylindricalCoordinateSystem
#    ../segmentation/Patch ../ReliefMap/UnrolledMap ../segmentation/DefectSegmentation
#    ../segmentation/SegmentationAbstract ../Centerline/SimpleNormalAccumulator
#     Centerline/CenterlineMultiCore Centerline/Centerline )
#TARGET_LINK_LIBRARIES(segunroll_mesh ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} ${PCLLib} ${CMAKE_THREAD_LIBS_INIT} ${GSL_LIBRARIES})


ADD_EXECUTABLE(main Main
      UnrolledMap)
TARGET_LINK_LIBRARIES(main ${DGTAL_LIBRARIES}) 