cmake_minimum_required(VERSION 3.5)
project(loop_weaver)

# Set the C++ standard to C++17
#set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto" )
#add_definitions(-DSAVE_STATUS_REMOVE)

set(IMGUI_DIR "/Users/134750/Desktop/imgui/" CACHE INTERNAL "")
set(TRACEFIELD_DIR = "../TracingFieldLib/include" CACHE INTERNAL "")
set(GLDRAW_DIR = "../TracingFieldLib/include/GLFunctions/" CACHE INTERNAL "")
set(EIGENPATH  = "../TracingFieldLib/external_lib/eigen-3.4.0/" CACHE INTERNAL "")


#set(SANITIZE_FLAGS "-fsanitize=address -fsanitize=undefined")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${SANITIZE_FLAGS} -g -O1")
#set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${SANITIZE_FLAGS}")

set(IGLLIB  = "../TracingFieldLib/external_lib/libigl/include/" CACHE INTERNAL "")

include_directories(${TRACEFIELD_DIR})
include_directories(${IMGUI_DIR})
include_directories(${IMGUI_DIR}/backends)
include_directories(/opt/homebrew/include)
include_directories(${EIGENPATH})

#set(GUROBI_PATH  = "/Library/gurobi951/macos_universal2/" CACHE INTERNAL "")
#set(GUROBI_COMPILER     = gurobi_c++)
#set(GUROBI_LIB          = gurobi95)
#set(IGLLIB  = "../TracingFieldLib/external_lib/libigl/include/" CACHE INTERNAL "")

#include_directories(${GUROBI_PATH})
#include_directories(${GUROBI_PATH}/include)

#if(USE_DIRECTIONAL)
#    include_directories(${DIRECTIONALPATH})
#endif()

include_directories(${IGLLIB})

# Find GLEW (installed via Homebrew)
find_library(GLEW_LIB NAMES GLEW glew PATHS /opt/homebrew/lib /usr/local/lib)
find_path(GLEW_INCLUDE_DIR GL/glew.h PATHS /opt/homebrew/include /usr/local/include)
if(NOT GLEW_LIB OR NOT GLEW_INCLUDE_DIR)
  message(FATAL_ERROR "GLEW not found")
endif()

add_executable(${PROJECT_NAME} main.cpp
               ${IMGUI_DIR}/imgui.cpp
               ${IMGUI_DIR}/imgui_demo.cpp
               ${IMGUI_DIR}/imgui_draw.cpp
               ${IMGUI_DIR}/imgui_tables.cpp
               ${IMGUI_DIR}/imgui_widgets.cpp
               ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
               ${IMGUI_DIR}/backends/imgui_impl_opengl2.cpp
               )

target_include_directories(${PROJECT_NAME} PRIVATE ${GLEW_INCLUDE_DIR})


target_link_libraries(${PROJECT_NAME}
    ${GLEW_LIB} 
    "-framework OpenGL"
    "-framework Cocoa"
    "-framework IOKit"
    "-framework CoreVideo"
    "-L/opt/homebrew/lib"
    "-lglfw"
)

#target_link_libraries(${PROJECT_NAME}
#"/Library/gurobi951/macos_universal2/lib/libgurobi_c++.a"
#"/Library/gurobi951/macos_universal2/lib/libgurobi95.dylib"
#)