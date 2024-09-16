cmake_minimum_required(VERSION 3.15)

project(Gravitacek2 VERSION 0.1)

add_executable(Gravitacek2 main.cpp)

add_library(setup STATIC ${PROJECT_SOURCE_DIR}/src/gravitacek2/setup.cpp)

target_include_directories(setup PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(Gravitacek2 PUBLIC setup)