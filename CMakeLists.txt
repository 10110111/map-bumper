cmake_minimum_required(VERSION 3.10.0)
project(equirect2cubemap VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type -Wall -Wextra -fno-strict-aliasing")
endif()

include(CheckIncludeFileCXX)
check_include_file_cxx("glm/glm.hpp" HAVE_GLM)
if(NOT HAVE_GLM)
    unset(HAVE_GLM CACHE)
    message(FATAL_ERROR "GLM was not found")
endif()

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
find_package(Qt6Core 6.3 REQUIRED)
find_package(Qt6Gui 6.3 REQUIRED)

add_executable(equirect2cubemap equirect2cubemap.cpp)
target_link_libraries(equirect2cubemap Qt6::Core Qt6::Gui)

add_executable(equirect2hips equirect2hips.cpp healpix.cpp)
target_link_libraries(equirect2hips Qt6::Core Qt6::Gui)

add_executable(wms-tiles-to-hips wms-tiles-to-hips.cpp healpix.cpp)
target_link_libraries(wms-tiles-to-hips Qt6::Core Qt6::Gui)

add_executable(heightmap2normalmap heightmap2normalmap.cpp)
target_link_libraries(heightmap2normalmap Qt6::Core Qt6::Gui)

add_executable(heightmap2horizonmap heightmap2horizonmap.cpp)
target_link_libraries(heightmap2horizonmap Qt6::Core Qt6::Gui)

add_executable(lroc-raw-to-srgb lroc-raw-to-srgb.cpp)
target_link_libraries(lroc-raw-to-srgb Qt6::Core Qt6::Gui)

add_executable(lroc-emp-to-srgb lroc-emp-to-srgb.cpp)
target_link_libraries(lroc-emp-to-srgb Qt6::Core Qt6::Gui)

add_executable(wms-to-srgb wms-to-srgb.cpp)
target_link_libraries(wms-to-srgb Qt6::Core Qt6::Gui)

add_executable(cgi-moon-kit-color-to-srgb cgi-moon-kit-color-to-srgb.cpp)
target_link_libraries(cgi-moon-kit-color-to-srgb Qt6::Core Qt6::Gui)

add_executable(moonkit-poles-recolor moonkit-poles-recolor.cpp)
target_link_libraries(moonkit-poles-recolor Qt6::Core Qt6::Gui)

add_executable(combine-lroc-and-moonkit combine-lroc-and-moonkit.cpp)
target_link_libraries(combine-lroc-and-moonkit Qt6::Core Qt6::Gui)

add_executable(sldem-img2png sldem-img2png.cpp)
target_link_libraries(sldem-img2png Qt6::Core Qt6::Gui)
target_compile_options(sldem-img2png PRIVATE -fno-strict-aliasing)
