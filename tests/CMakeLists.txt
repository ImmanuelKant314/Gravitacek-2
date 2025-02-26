# minimum version for CMake
cmake_minimum_required(VERSION 3.15)

include(GoogleTest)
include(CTest)

# move folder
# file(COPY ${PROJECT_SOURCE_DIR}/tests/test-WeylSpacetime DESTINATION ${CMAKE_BINARY_DIR}/tests)

# add executable 
add_executable(test_integratorcomponents test_integratorcomponents.cpp)
# add_executable(test_integrator test_integrator.cpp)
add_executable(test_stepcontrollertypes test_stepcontrollertypes.cpp)
add_executable(test_steppertypes test_steppertypes.cpp)
add_executable(test_spacetimes test_spacetimes.cpp)
add_executable(test_weylspacetimes test_weylspacetimes.cpp)
add_executable(test_majumdarpapapetrouwspacetimes test_majumdarpapapetrouwspacetimes.cpp)
add_executable(test_mymath test_mymath.cpp)
add_executable(test_chaos test_chaos.cpp)

# link libraries
target_link_libraries(test_integratorcomponents PRIVATE gtest gtest_main integrator setup)
# target_link_libraries(test_integrator PRIVATE gtest gtest_main integrator setup)
target_link_libraries(test_stepcontrollertypes PRIVATE gtest gtest_main integrator setup)
target_link_libraries(test_steppertypes PRIVATE gtest gtest_main integrator setup)
target_link_libraries(test_spacetimes PRIVATE gtest gtest_main geomotion setup)
target_link_libraries(test_weylspacetimes PRIVATE gtest gtest_main geomotion setup)
target_link_libraries(test_majumdarpapapetrouwspacetimes PRIVATE gtest gtest_main geomotion setup)
target_link_libraries(test_mymath PRIVATE gtest gtest_main setup mymath)
target_link_libraries(test_chaos PRIVATE gtest gtest_main setup chaos)

# add tests
add_test(
    NAME    test_integratorcomponents
    COMMAND test_integratorcomponents
)

# add_test(
#     NAME    test_integrator
#     COMMAND test_integrator
# )

add_test(
    NAME    test_stepcontrollertypes
    COMMAND test_stepcontrollertypes
)

add_test(
    NAME    test_steppertypes
    COMMAND test_steppertypes
)

add_test(
    NAME    test_spacetimes
    COMMAND test_spacetimes
)

add_test(
    NAME    test_weylspacetimes
    COMMAND test_weylspacetimes
)

add_test(
    NAME    test_majumdarpapapetrouwspacetimes
    COMMAND test_majumdarpapapetrouwspacetimes
)

add_test(
    NAME    test_mymath
    COMMAND test_mymath
)

add_test(
    NAME    test_chaos
    COMMAND test_chaos
)

# discover tests
gtest_discover_tests(test_integratorcomponents)
# gtest_discover_tests(test_integrator)
gtest_discover_tests(test_stepcontrollertypes)
gtest_discover_tests(test_steppertypes)
gtest_discover_tests(test_spacetimes)
gtest_discover_tests(test_weylspacetimes)
gtest_discover_tests(test_majumdarpapapetrouwspacetimes)
gtest_discover_tests(test_mymath)
gtest_discover_tests(test_chaos)
