cmake_minimum_required(VERSION 3.2)

# Add the path to the custom CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")


#add_executable(kmouse kinect_touch.cpp)
#set(OpenCV_DIR ../build)
#find_package(OpenCV REQUIRED)
# find_package(OpenNI REQUIRED)
# include_directories(${OPENNI_INCLUDE_DIRS})

# include_directories(${OpenCV_INCLUDE_DIRS})
# include_directories(
#   "/usr/include/ni"
#   "src/TUIO"
#   "src/oscpack"
# )

add_executable(kmouse kinect_mouse_mm.c)

find_package(Threads REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(freenect REQUIRED)

if (NOT GLUT_FOUND)
  message(FATAL_ERROR "GLUT not found!")
endif()
if (NOT freenect_FOUND)
  message(FATAL_ERROR "libfreenect not found!")
endif()

message(STATUS "GLUT library: ${GLUT_LIBRARY}")
message(STATUS "libfreenect library: ${FREENECT_LIBRARIES}")

include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${USB_INCLUDE_DIRS} ${FREENECT_INCLUDE_DIRS})
target_link_libraries(kmouse ${FREENECT_LIBRARIES} X11 Xtst ncurses ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} m)

install (TARGETS kmouse
  DESTINATION bin)
