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

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
find_package(Qt6Core 6.3 REQUIRED)
find_package(Qt6Gui 6.3 REQUIRED)

find_package(TIFF REQUIRED)
check_include_file_cxx("geotiff/geotiff.h" HAVE_GEOTIFF)
if(NOT HAVE_GEOTIFF)
    unset(HAVE_GEOTIFF CACHE)
    message(FATAL_ERROR "geotiff include header was not found")
endif()
find_library(GEOTIFF_LIB geotiff)
if(NOT GEOTIFF_LIB)
    unset(GEOTIFF_LIB CACHE)
    message(FATAL_ERROR "geotiff library was not found")
endif()
find_library(PROJ_LIB proj)
if(NOT PROJ_LIB)
    unset(PROJ_LIB CACHE)
    message(FATAL_ERROR "proj library was not found")
endif()

include(FindPkgConfig)
pkg_check_modules(CFITSIO cfitsio)
if(CFITSIO_FOUND)
    include_directories(${CFITSIO_INCLUDE_DIRS})
else()
    message(WARNING "cfitsio library wasn't found, so sun-img-to-map won't be built")
endif()

add_executable(equirect2cubemap equirect2cubemap.cpp)
target_link_libraries(equirect2cubemap Qt6::Core Qt6::Gui)

add_executable(equirect2hips equirect2hips.cpp healpix.cpp hips.cpp timing.cpp)
target_link_libraries(equirect2hips Qt6::Core Qt6::Gui)

add_executable(complete-hips-tiles complete-hips-tiles.cpp healpix.cpp hips.cpp timing.cpp)
target_link_libraries(complete-hips-tiles Qt6::Core Qt6::Gui)

add_executable(lroc-to-hips lroc-to-hips.cpp healpix.cpp hips.cpp timing.cpp)
target_link_libraries(lroc-to-hips Qt6::Core Qt6::Gui)

add_executable(gld100-to-hips-normalmap gld100-to-hips-normalmap.cpp healpix.cpp hips.cpp timing.cpp)
target_link_libraries(gld100-to-hips-normalmap Qt6::Core Qt6::Gui)

add_executable(sldem-to-hips-normalmap sldem-to-hips-normalmap.cpp healpix.cpp hips.cpp timing.cpp)
target_link_libraries(sldem-to-hips-normalmap Qt6::Core Qt6::Gui TIFF::TIFF ${GEOTIFF_LIB} ${PROJ_LIB})

add_executable(cubemap2hips-horizonmap cubemap2hips-horizonmap.cpp healpix.cpp hips.cpp timing.cpp)
target_link_libraries(cubemap2hips-horizonmap Qt6::Core Qt6::Gui)

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

add_executable(ldem-img2png ldem-img2png.cpp)
target_link_libraries(ldem-img2png Qt6::Core Qt6::Gui)

if(CFITSIO_FOUND)
    add_executable(sun-img-to-map sun-img-to-map.cpp)
    target_link_libraries(sun-img-to-map TIFF::TIFF ${CFITSIO_LIBRARIES})

    add_executable(fits2tiff fits2tiff.cpp)
    target_link_libraries(fits2tiff TIFF::TIFF ${CFITSIO_LIBRARIES})
endif()
