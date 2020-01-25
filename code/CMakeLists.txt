cmake_minimum_required(VERSION 3.15)
project(code)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

add_executable(code main.cpp Linag/vector.h Linag/densematrix.h Linag/sparsematrix.h "eigen.cpp")

target_link_libraries(code Threads::Threads)