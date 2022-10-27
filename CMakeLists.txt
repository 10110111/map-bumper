cmake_minimum_required(VERSION 3.10.0)
project(ApertureDiffraction VERSION 0.0.1 LANGUAGES CXX)

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
find_package(Qt5 5.9 REQUIRED Core Gui)

add_executable(equirect2cubemap equirect2cubemap.cpp)
target_link_libraries(equirect2cubemap Qt5::Core Qt5::Gui)

add_executable(equirect2normalmap equirect2normalmap.cpp)
target_link_libraries(equirect2normalmap Qt5::Core Qt5::Gui)

add_executable(equirect2horizonmap equirect2horizonmap.cpp)
target_link_libraries(equirect2horizonmap Qt5::Core Qt5::Gui)
