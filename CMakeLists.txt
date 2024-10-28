# minimum version for CMake
cmake_minimum_required(VERSION 3.15)

# project name and version
project(Gravitacek2 VERSION 0.1)

# minimum C++ standard required
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# add executable 
add_executable(Gravitacek2 main.cpp)

# add subdirectories
add_subdirectory(${PROJECT_SOURCE_DIR}/external/googletest)
add_subdirectory(${PROJECT_SOURCE_DIR}/tests)

# directory shortcuts
set(INTEGRATOR_DIR ${PROJECT_SOURCE_DIR}/src/gravitacek2/integrator)
set(GEOMOTION_DIR ${PROJECT_SOURCE_DIR}/src/gravitacek2/geomotion)
set(MYMATH_DIR ${PROJECT_SOURCE_DIR}/src/gravitacek2)
set(INTF_DIR ${PROJECT_SOURCE_DIR}/src/interface)

# add libraries
add_library(setup STATIC ${PROJECT_SOURCE_DIR}/src/gravitacek2/setup.cpp)
add_library(integrator
            STATIC 
            ${INTEGRATOR_DIR}/odesystem.cpp
            ${INTEGRATOR_DIR}/event.cpp
            ${INTEGRATOR_DIR}/stepperbase.cpp
            ${INTEGRATOR_DIR}/stepcontrollerbase.cpp
            ${INTEGRATOR_DIR}/integrator.cpp
            ${INTEGRATOR_DIR}/steppers/rk4.cpp
            ${INTEGRATOR_DIR}/steppers/dopr853.cpp
            ${INTEGRATOR_DIR}/stepcontrollers/standardstepcontroller.cpp
            ${INTEGRATOR_DIR}/stepcontrollers/stepcontrollernr.cpp
            ${INTEGRATOR_DIR}/odesystems.cpp)
add_library(geomotion
            STATIC 
            ${GEOMOTION_DIR}/geomotion.cpp
            ${GEOMOTION_DIR}/weyl.cpp
            ${GEOMOTION_DIR}/spacetimes/schwarzschild.cpp
            ${GEOMOTION_DIR}/spacetimes/weylschwarzschild.cpp
            ${GEOMOTION_DIR}/spacetimes/bachweylring.cpp)
add_library(mymath
            STATIC
            ${MYMATH_DIR}/mymath.cpp)
add_library(interface
            STATIC
            ${INTF_DIR}/interface.cpp)


# include directiories
target_include_directories(setup PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(integrator PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(geomotion PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(mymath PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(interface PUBLIC ${PROJECT_SOURCE_DIR}/include)

# link libraries
target_link_libraries(Gravitacek2 PUBLIC setup interface)
target_link_libraries(integrator PUBLIC setup)
target_link_libraries(geomotion PUBLIC setup integrator mymath)
target_link_libraries(mymath PUBLIC setup)
target_link_libraries(interface PUBLIC setup geomotion)

# testing 
enable_testing()