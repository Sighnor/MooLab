cmake_minimum_required(VERSION 3.10)
project(Engine)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(OpenCV REQUIRED)

set(CMAKE_CXX_STANDARD 17)

ADD_EXECUTABLE(Engine main.cpp array.hpp camera.hpp character.hpp controller.hpp database.hpp global.hpp light.hpp material.hpp mat.hpp mesh.hpp model.hpp motion.hpp nnet.hpp quat.hpp renderer.hpp robot.hpp shader.hpp texture.hpp vec.hpp)
target_link_libraries(Engine ${OpenCV_LIBRARIES})
