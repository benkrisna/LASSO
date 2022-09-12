cmake_minimum_required(VERSION 3.15)

project(LASSO VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)

#add_subdirectory(WINDOW)
#list(APPEND MYLIBS WINDOW)
add_library(grid grid.cpp)
list(APPEND SRC_LIBS grid)

add_executable(lasso lasso.cpp)

#find_package(OpenGL REQUIRED)
#find_package(GLUT REQUIRED)
find_package(HDF5 REQUIRED)
find_package(toml11 REQUIRED)

#include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

target_link_libraries(lasso ${SRC_LIBS} ${HDF5_CXX_LIBRARIES} ${TOML11_LIBRARIES}) # ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${MYLIBS})
