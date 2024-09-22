cmake_minimum_required(VERSION 3.15)

project(Gravitacek2 VERSION 0.1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(Gravitacek2 main.cpp)
add_subdirectory(${PROJECT_SOURCE_DIR}/external/googletest)
add_subdirectory(${PROJECT_SOURCE_DIR}/tests)

set(ODESOLVER_DIR ${PROJECT_SOURCE_DIR}/src/gravitacek2/odesolver)

add_library(setup STATIC ${PROJECT_SOURCE_DIR}/src/gravitacek2/setup.cpp)
add_library(odesolver 
            STATIC 
            ${ODESOLVER_DIR}/ode.cpp
            ${ODESOLVER_DIR}/stepper.cpp
            ${ODESOLVER_DIR}/stepcontroller.cpp
            ${ODESOLVER_DIR}/stepper_types/rk4.cpp
            ${ODESOLVER_DIR}/stepcontroller_types/standardstepcontroller.cpp)

target_include_directories(setup PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(odesolver PRIVATE ${ODESOLVER_DIR})

target_link_libraries(Gravitacek2 PUBLIC setup)
target_link_libraries(odesolver PUBLIC setup)

enable_testing()