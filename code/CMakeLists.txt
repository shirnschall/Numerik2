cmake_minimum_required(VERSION 3.15)
project(code)

set(CMAKE_CXX_STANDARD 14)
set ( CMAKE_C_FLAGS "pthread")
find_package( Threads )

add_executable(code main.cpp Linag/vector.h Linag/densematrix.h Linag/sparsematrix.h)

target_link_libraries(code ${CMAKE_THREAD_LIBS_INIT} )