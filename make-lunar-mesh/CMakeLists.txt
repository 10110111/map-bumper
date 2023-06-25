cmake_minimum_required(VERSION 3.10.0)
project(make-lunar-mesh VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type -Wall -Wextra -fno-strict-aliasing")
endif()

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
find_package(Qt6Core 6.3 REQUIRED)
find_package(Qt6Gui 6.3 REQUIRED)

add_executable(make-lunar-mesh make-lunar-mesh.cpp)
target_link_libraries(make-lunar-mesh Qt6::Core Qt6::Gui)

add_executable(obj2bin obj2bin.cpp)
target_link_libraries(obj2bin)
