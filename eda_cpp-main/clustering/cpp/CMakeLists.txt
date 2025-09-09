cmake_minimum_required(VERSION 3.10)
project(clustering LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)


add_executable(test
    test.cpp
    src/utils.cpp
    src/matrix.cpp
    src/cluster.cpp
    src/simsearch.cpp
)