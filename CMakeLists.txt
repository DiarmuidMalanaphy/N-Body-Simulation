cmake_minimum_required(VERSION 3.5)
project(YourProject)

# Specify the paths to GLFW and GLAD
set(GLFW_INCLUDE_DIR "C:/Users/diarmuid/Downloads/glfw-3.4.bin.WIN64/glfw-3.4.bin.WIN64/include")
set(GLFW_LIBRARY_DIR "C:/Users/diarmuid/Downloads/glfw-3.4.bin.WIN64/glfw-3.4.bin.WIN64/lib-vc2022")
set(GLAD_INCLUDE_DIR "C:/Users/diarmuid/Downloads/glad/include")
set(GLM_INCLUDE_DIR "C:/Users/diarmuid/Downloads/glm-1.0.1-light")
set(GLAD_SRC "C:/Users/diarmuid/Downloads/glad/src/glad.c")
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
# Include GLFW and GLAD headers
include_directories(${GLFW_INCLUDE_DIR} ${GLAD_INCLUDE_DIR} ${GLM_INCLUDE_DIR})




# Find GLFW package
find_library(GLFW_LIBRARY
    NAMES glfw3 glfw
    PATHS ${GLFW_LIBRARY_DIR}
    NO_DEFAULT_PATH
    REQUIRED
)

# Find OpenGL package
find_package(OpenGL REQUIRED)

# Add your executable and include GLAD source
add_executable(
    YourExecutable
    code/ForceInteraction.cpp
    code/Bound.cpp
    code/SphereRender.cpp
    code/Particle.cpp
    code/Node.cpp
    code/BarnesHut.cpp
    
    code/main.cpp
    ${GLAD_SRC}
    )
# add_executable(YourExecutable main.cpp ${GLAD_SRC})

# Include OpenGL headers
include_directories(${OPENGL_INCLUDE_DIR})

# Link GLFW and OpenGL libraries
target_link_libraries(YourExecutable ${GLFW_LIBRARY} ${OPENGL_LIBRARIES})

# Specify C++11 standard
set_target_properties(YourExecutable PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

# Set the runtime library to Multi-threaded DLL (/MD)
if(MSVC)
    set_property(TARGET YourExecutable PROPERTY
        MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
endif()

# Add GLFW DLL to the output directory (if using shared library)
if(WIN32)
    add_custom_command(TARGET YourExecutable POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${GLFW_LIBRARY_DIR}/glfw3.dll"
        $<TARGET_FILE_DIR:YourExecutable>
    )
endif()


