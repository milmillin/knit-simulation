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
include_directories("./third_party/")

# Add project files
if (${CMAKE_VERSION} VERSION_EQUAL "3.11.0"
    OR ${CMAKE_VERSION} VERSION_GREATER "3.11.0")
  message("-- Using GLOB_RECURSE to find source files")
  file(GLOB_RECURSE SRCFILES CONFIGURE_DEPENDS src/*.cpp src/*.h)
else()
  file(GLOB SRCFILES src/*.cpp src/*/*.cpp src/*/*/*.cpp third_party/)
endif()

# Enable warnings
if (MSVC)
  # For Microsoft compilers
  set_source_files_properties(
    ${SRCFILES}
    PROPERTIES
    COMPILE_FLAGS "/W4"
  )
else()
  # For other compilers
  set_source_files_properties(
    ${SRCFILES}
    PROPERTIES
    COMPILE_FLAGS "-Wall -Wextra -pedantic"
  )
endif()

# Link libraries
set(LINK_LIBS igl::core igl::opengl_glfw igl::opengl_glfw_imgui)
if (NOT MSVC)
    list(APPEND LINK_LIBS stdc++fs)
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
set(ENABLE_EASY_PROFILER TRUE)
if (ENABLE_EASY_PROFILER)
    add_definitions(-DUSE_EASY_PROFILER)
    set(EASY_PROFILER_DIR "${PROJECT_SOURCE_DIR}/../easy-profiler") # change this to easy-profiler-dir
    message(STATUS "Using easy-profiler from ${EASY_PROFILER_DIR}")
	list(APPEND CMAKE_PREFIX_PATH "${EASY_PROFILER_DIR}/lib/cmake/easy_profiler")
    list(APPEND LINK_LIBS easy_profiler)
	find_package(easy_profiler REQUIRED)
endif()


# AABB tree
include_directories(third_party/aabbcc/src)
set(SRCFILES ${SRCFILES} third_party/aabbcc/src/AABB.cc)

# Threading
include_directories(third_party/CTPL)

# cxxopts
include_directories(third_party/cxxopts/include)

# Linking
add_executable(${PROJECT_NAME} ${SRCFILES})
target_link_libraries(${PROJECT_NAME} ${LINK_LIBS})

# Copy easy_profiler.dll to binary folder
if (ENABLE_EASY_PROFILER)
if ("${CMAKE_SYSTEM_NAME}" MATCHES "[Ww]indows")
    add_custom_command(TARGET ${PROJECT_NAME}
      POST_BUILD
      COMMAND
        ${CMAKE_COMMAND} -E copy_if_different
        "${EASY_PROFILER_DIR}/bin/easy_profiler.dll"
        $<TARGET_FILE_DIR:${PROJECT_NAME}>
    )
endif()
endif()
