cmake_minimum_required(VERSION 3.18)
project(InflatorConsumer LANGUAGES CXX)

find_package(inflator CONFIG REQUIRED)
add_executable(consumer ../../src/api/test_api.cpp)
target_link_libraries(consumer PRIVATE inflator::inflator)
enable_testing()
add_test(NAME installed_inflator COMMAND consumer)
