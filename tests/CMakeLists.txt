# minimum version for CMake
cmake_minimum_required(VERSION 3.15)

include(GoogleTest)
include(CTest)

# add executable 
add_executable(test_odesolver test_odesolver.cpp)
add_executable(test_integrator test_integrator.cpp)
add_executable(test_stepcontrollertypes test_stepcontrollertypes.cpp)
add_executable(test_steppertypes test_steppertypes.cpp)

# link libraries
target_link_libraries(test_odesolver PRIVATE gtest gtest_main odesolver setup)
target_link_libraries(test_integrator PRIVATE gtest gtest_main odesolver setup)
target_link_libraries(test_stepcontrollertypes PRIVATE gtest gtest_main odesolver setup)
target_link_libraries(test_steppertypes PRIVATE gtest gtest_main odesolver setup)

# add tests
add_test(
    NAME    test_odesolver 
    COMMAND test_odesolver
)

add_test(
    NAME    test_integrator
    COMMAND test_integrator
)

add_test(
    NAME    test_stepcontrollertypes
    COMMAND test_stepcontrollertypes
)

add_test(
    NAME    test_steppertypes
    COMMAND test_steppertypes
)

gtest_discover_tests(test_odesolver)
gtest_discover_tests(test_integrator)
gtest_discover_tests(test_stepcontrollertypes)
gtest_discover_tests(test_steppertypes)