cmake_minimum_required(VERSION 3.15)

project(LASSO VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)

#add_subdirectory(WINDOW)
#list(APPEND MYLIBS WINDOW)

add_executable(lasso lasso.cpp)

#find_package(OpenGL REQUIRED)
#find_package(GLUT REQUIRED)
find_package(HDF5 REQUIRED)

#include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

target_link_libraries(lasso ${HDF5_CXX_LIBRARIES}) # ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${MYLIBS})
