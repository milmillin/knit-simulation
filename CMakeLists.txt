cmake_minimum_required(VERSION 3.1)
project(knit-simulator)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_BUILD_TYPE Debug)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories("./")

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use IMGUI"          ON)

find_package(LIBIGL REQUIRED QUIET)

# glm
find_package(GLM REQUIRED)
include_directories(${GLM_INCLUDE_DIRS})
link_directories(${GLM_LIBRARY_DIRS})
add_definitions(${GLM_DEFINITIONS})

# Easy Profiler
list(APPEND CMAKE_PREFIX_PATH "C:/Users/milmi/Documents/1_Projects/KnittingSimulation/easy-profiler/lib/cmake/easy_profiler")
find_package(easy_profiler REQUIRED)

# Add your project files
file(GLOB SRCFILES *.cpp */*.cpp simulator/*/*.cpp)
add_executable(${PROJECT_NAME} ${SRCFILES})
target_link_libraries(${PROJECT_NAME} easy_profiler igl::core igl::opengl_glfw igl::opengl_glfw_imgui stdc++fs)
