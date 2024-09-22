# minimum version for CMake
cmake_minimum_required(VERSION 3.15)

# add executable 
add_executable(test_odesolver test_odesolver.cpp)

# link libraries
target_link_libraries(test_odesolver PRIVATE gtest gtest_main odesolver setup)

# add tests
add_test(
    NAME    test_odesolver 
    COMMAND test_odesolver
)

include(GoogleTest)
gtest_discover_tests(test_odesolver)