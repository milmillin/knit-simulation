cmake_minimum_required(VERSION 3.1)
project(knit-simulator)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("./src/")

# Add project files
if (${CMAKE_VERSION} VERSION_EQUAL "3.11.0"
    OR ${CMAKE_VERSION} VERSION_GREATER "3.11.0")
  message("-- Using GLOB_RECURSE to find source files")
  file(GLOB_RECURSE SRCFILES CONFIGURE_DEPENDS src/*.cpp src/*.h)
else()
  file(GLOB SRCFILES src/*.cpp src/*/*.cpp src/*/*/*.cpp third_party/)
endif()

# spdlog
include_directories("./third_party/spdlog/include/")

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use IMGUI"          ON)

find_package(LIBIGL REQUIRED QUIET)

# GLM
find_package(GLM REQUIRED)
include_directories(${GLM_INCLUDE_DIRS})
link_directories(${GLM_LIBRARY_DIRS})
add_definitions(${GLM_DEFINITIONS})

# Easy Profiler
list(APPEND CMAKE_PREFIX_PATH "C:/Users/milmi/Documents/1_Projects/KnittingSimulation/easy-profiler/lib/cmake/easy_profiler")
find_package(easy_profiler REQUIRED)

# AABB tree
include_directories(third_party/aabbcc/src)
set(SRCFILES ${SRCFILES} third_party/aabbcc/src/AABB.cc)

# Threading
include_directories(third_party/CTPL)

# Linking
add_executable(${PROJECT_NAME} ${SRCFILES})
target_link_libraries(${PROJECT_NAME} easy_profiler igl::core igl::opengl_glfw igl::opengl_glfw_imgui) # Add stdc++fs
