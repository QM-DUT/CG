cmake_minimum_required(VERSION 2.8.12)
project(Assignment4)

find_package(OpenGL REQUIRED)
find_package(GLU REQUIRED)

# Suppress warnings of the deprecation of glut functions on macOS.
if(APPLE)
 add_definitions(-Wno-deprecated-declarations)
endif()

### Output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

### Compilation flags: adapt to your needs ###
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") ### Enable parallel compilation
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

### Add src to the include directories
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")

### Add OpenGL
# set(INCLUDE_DIRS ${OPENGL_INCLUDE_DIR})
# set(LIBRARIES ${OPENGL_LIBRARIES})

#链接动态链接库
#set(LINK_DIR_freetype /usr/local/freetype/lib)
#link_directories(${LINK_DIR_freetype})#添加动态库的路径，相当于-L
#link_libraries(freetype)#已被废弃，和target_link_libraries一个作用 -lfreetype
#包含头文件
include_directories("/usr/local/freetype/include/freetype2/")#相当于-I 

#添加opencv
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

### Include Eigen for linear algebra
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../ext/eigen")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../ext/glm")
### Compile GLFW3 statically
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_INSTALL OFF CACHE BOOL " " FORCE)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../ext/glfw" "glfw")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../ext/glfw/include")
set(LIBRARIES "glfw" ${GLFW_LIBRARIES})

### On windows, you also need glew
if((UNIX AND NOT APPLE) OR WIN32)
  set(GLEW_INSTALL OFF CACHE BOOL " " FORCE)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../ext/glew" "glew")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../ext/glew/include")
  list(APPEND LIBRARIES "glew")
endif()

if(APPLE)
list(APPEND LIBRARIES "-framework OpenGL")
endif()

# ### Compile all the cpp files in src
# #main其实就是路径，现在就是得到编译成最终二进制文件所需的所有cpp路径
# file(GLOB main "${CMAKE_CURRENT_SOURCE_DIR}/src/main*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/extra/main*.cpp")
# file(GLOB helper "${CMAKE_CURRENT_SOURCE_DIR}/src/Helpers*.cpp")
# foreach(mainfile IN LISTS main)
#     # Get file name without directory
#     file(GLOB main_helper ${mainfile} ${helper})
#     #message(${main_helper})
#     get_filename_component(mainname ${mainfile} NAME_WE)
#     add_executable(${mainname} ${main_helper})
#     target_link_libraries(${mainname} ${LIBRARIES} ${OPENGL_LIBRARIES})
# endforeach()

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(${PROJECT_NAME}_bin ${SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${LIBRARIES} ${OPENGL_LIBRARIES} freetype ${OpenCV_LIBS})