cmake_minimum_required(VERSION 3.5)
project(loop_weaver)

# ==============================
# C++ Standard
# ==============================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ==============================
# Build mode
# ==============================
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type")
endif()
option(ENABLE_SANITIZER "Enable sanitizers" OFF)

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
set(IMGUI_DIR "../tools/imgui")
set(TRACEFIELD_DIR "../CurveSurfacingLib/TracingFieldLib/include")
set(EIGENPATH "../CurveSurfacingLib/TracingFieldLib/external_lib/eigen-3.4.0")
set(IGLLIB "../libigl/include")
set(CURVE_SURFACING_CORE_DIR "${CMAKE_SOURCE_DIR}/../Refactor/CurveSurfacingLib-V2/Curve-Surfacing-Core")

set(LUT_PATH "../../../include/Lut.txt")
add_definitions(-DLUT_PATH=\"${LUT_PATH}\")

# ==============================
# Dependencies
# ==============================
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)

# Use local Eigen instead of system package
if(NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/${EIGENPATH}"
    )
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

# ==============================
# Includes
# ==============================
target_include_directories(${PROJECT_NAME} PRIVATE
    ${TRACEFIELD_DIR}
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
    ${CURVE_SURFACING_CORE_DIR}/include
    ${CURVE_SURFACING_CORE_DIR}/include/curve_surfacing
    ${CURVE_SURFACING_CORE_DIR}/include/curve_surfacing/CycleDiscovery/tetgen
    ${CURVE_SURFACING_CORE_DIR}/build/_deps/triangle-src
    ${CURVE_SURFACING_CORE_DIR}/build/_deps/triangle-build
    ${CURVE_SURFACING_CORE_DIR}/external/Directional/include
    ${IGLLIB}
)

# Ensure Eigen headers are available (must come before other includes)
target_include_directories(${PROJECT_NAME} BEFORE PRIVATE ${EIGENPATH})

# ==============================
# Link libraries
# no libcurvesurfacing-core.a — headers self-include their .cpp files
# ==============================
target_link_libraries(${PROJECT_NAME}
    ${CURVE_SURFACING_CORE_DIR}/build/libtetgen.a
    ${CURVE_SURFACING_CORE_DIR}/build/_deps/triangle-build/libtriangle.a
    OpenGL::GL
    OpenGL::GLU
    GLEW::GLEW
    glfw
    Eigen3::Eigen
)