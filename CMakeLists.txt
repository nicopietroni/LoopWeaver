cmake_minimum_required(VERSION 3.5)
project(loop_weaver)

# ==============================
# C++ Standard
# ==============================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ==============================
# ===== USER SWITCHES (EDIT ONLY THESE) =====
# ==============================
set(BUILD_MODE "Release")     # Debug | Release
set(ENABLE_SANITIZER OFF)      # ON | OFF

# ==============================
# Apply build mode
# ==============================
set(CMAKE_BUILD_TYPE ${BUILD_MODE} CACHE STRING "" FORCE)

message(STATUS "Build mode: ${CMAKE_BUILD_TYPE}")
message(STATUS "Sanitizer: ${ENABLE_SANITIZER}")

# ==============================
# Base flags
# ==============================
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

# ==============================
# Sanitizer flags (can be used in BOTH Debug and Release)
# ==============================
if(ENABLE_SANITIZER)
  add_compile_options(
    -g
    -O1
    -fno-omit-frame-pointer
    -fsanitize=address,undefined
    -fsanitize-address-use-after-scope
    -fno-optimize-sibling-calls
  )

  add_link_options(
    -fsanitize=address,undefined
  )
endif()

# ==============================
# Paths
# ==============================
set(IMGUI_DIR "/Users/134750/Desktop/imgui/" CACHE INTERNAL "")
set(TRACEFIELD_DIR "../TracingFieldLib/include" CACHE INTERNAL "")
set(GLDRAW_DIR "../TracingFieldLib/include/GLFunctions/" CACHE INTERNAL "")
set(EIGENPATH "../TracingFieldLib/external_lib/eigen-3.4.0/" CACHE INTERNAL "")
set(IGLLIB "../TracingFieldLib/external_lib/libigl/include/" CACHE INTERNAL "")

include_directories(
    ${TRACEFIELD_DIR}
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
    /opt/homebrew/include
    ${EIGENPATH}
    /Users/134750/Desktop/00_GeoProjects/CurveSurfacingLib/Curve-Surfacing-Core/include/curve_surfacing
    ${IGLLIB}
)

#include(libigl)
#igl_include(restricted triangle)

# ==============================
# Find GLEW
# ==============================
find_library(GLEW_LIB NAMES GLEW glew PATHS /opt/homebrew/lib /usr/local/lib)
find_path(GLEW_INCLUDE_DIR GL/glew.h PATHS /opt/homebrew/include /usr/local/include)

if(NOT GLEW_LIB OR NOT GLEW_INCLUDE_DIR)
  message(FATAL_ERROR "GLEW not found")
endif()

# ==============================
# Executable
# ==============================
add_executable(${PROJECT_NAME}
    main.cpp
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_demo.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl2.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE ${GLEW_INCLUDE_DIR})

# ==============================
# Static Libraries
# ==============================
target_link_libraries(${PROJECT_NAME}
   /Users/134750/Desktop/00_GeoProjects/CurveSurfacingLib/Curve-Surfacing-Core/build/libcurvesurfacing-core.a
   /Users/134750/Desktop/00_GeoProjects/TracingFieldLib/external_lib/TetGen-1.4.3/libtet.a
   /Users/134750/Desktop/00_GeoProjects/CurveSurfacingLib/Curve-Surfacing-Core/build/_deps/triangle-build/libtriangle.a
)

# ==============================
# System + OpenGL
# ==============================
target_link_libraries(${PROJECT_NAME}
    ${GLEW_LIB}
    "-framework OpenGL"
    "-framework Cocoa"
    "-framework IOKit"
    "-framework CoreVideo"
    "-L/opt/homebrew/lib"
    "-lglfw"
)
# ==============================