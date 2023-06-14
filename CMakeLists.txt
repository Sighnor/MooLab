cmake_minimum_required(VERSION 3.10)
project(Engine)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(RAYLIB_VERSION 4.2.0)
find_package(raylib ${RAYLIB_VERSION} QUIET) # QUIET or REQUIRED
if (NOT raylib_FOUND) # If there's none, fetch and build raylib
  include(FetchContent)
  FetchContent_Declare(
    raylib
    URL https://github.com/raysan5/raylib/archive/refs/tags/${RAYLIB_VERSION}.tar.gz
  )
  FetchContent_GetProperties(raylib)
  if (NOT raylib_POPULATED) # Have we downloaded raylib yet?
    set(FETCHCONTENT_QUIET NO)
    FetchContent_Populate(raylib)
    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) # don't build the supplied examples
    add_subdirectory(${raylib_SOURCE_DIR} ${raylib_BINARY_DIR})
  endif()
endif()

find_package( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(OpenCV REQUIRED)

set(CMAKE_CXX_STANDARD 17)

add_definitions(-D_DEFAULT_SOURCE)
add_definitions(-DRAYLIB_BUILD_MODE=RELEASE)
add_definitions(-DPLATFORM_DESKTOP)

add_executable(Engine main.cpp array.hpp camera.hpp character.hpp controller.hpp database.hpp global.hpp light.hpp material.hpp mat.hpp mesh.hpp model.hpp motion.hpp nnet.hpp quat.hpp renderer.hpp robot.hpp shader.hpp texture.hpp vec.hpp)
target_link_libraries(Engine ${OpenCV_LIBRARIES})
target_link_libraries(Engine raylib)
