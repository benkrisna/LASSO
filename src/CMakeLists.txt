cmake_minimum_required(VERSION 3.15)

project(LASSO VERSION 1.0)

set(CMAKE_CXX_STANDARD 14)

#add_subdirectory(WINDOW)
#list(APPEND MYLIBS WINDOW)

add_executable(lasso lasso.cpp)

#find_package(OpenGL REQUIRED)
#find_package(GLUT REQUIRED)

#include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})

target_link_libraries(lasso) # ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${MYLIBS})
