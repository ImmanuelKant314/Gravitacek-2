cmake_minimum_required(VERSION 3.15)

project(Gravitacek2 VERSION 0.1)

add_executable(Gravitacek2 main.cpp)
add_subdirectory(${PROJECT_SOURCE_DIR}/external/googletest)
add_subdirectory(${PROJECT_SOURCE_DIR}/tests)

add_library(setup STATIC ${PROJECT_SOURCE_DIR}/src/gravitacek2/setup.cpp)
add_library(ode STATIC ${PROJECT_SOURCE_DIR}/src/gravitacek2/odesolver/ode.cpp)

target_include_directories(setup PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ode PRIVATE ${PROJECT_SOURCE_DIR}/include/gravitacek2/odesolver)

target_link_libraries(Gravitacek2 PUBLIC setup)
target_link_libraries(ode PUBLIC setup)

enable_testing()