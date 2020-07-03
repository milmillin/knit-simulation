
if(GLM_FOUND)
    return()
endif()

find_path(GLM_INCLUDE_DIR glm/glm.hpp
    HINTS
        ${GLM_DIR}
        ENV GLM_DIR
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/glm
        ${CMAKE_SOURCE_DIR}/../glm
        ${CMAKE_SOURCE_DIR}/../../glm
        /usr
        /usr/local
        /usr/local/glm
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLM
    "\nGLM not found"
    GLM_INCLUDE_DIR)

if (GLM_FOUND)
    set(GLM_INCLUDE_DIRS ${GLM_INCLUDE_DIR})
endif()

mark_as_advanced(GLM_INCLUDE_DIR)

