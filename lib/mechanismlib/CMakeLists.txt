cmake_minimum_required(VERSION 3.13)

# Set CMAKE variables that can be called by the user to build specific sets
option(WITH_TESTS "Run tests after compilation" ON)

# set the project name
project(mechanismlibrary)
# It is a static library when applied with bindings, as the bindings are themselves a shared library
# this skips a linking step and keeps the python cmech library contained in a single file
add_library(cmechanism STATIC "common.h")
#target_compile_definitions(cmechanismlib PRIVATE EIGEN_RUNTIME_NO_MALLOC)

add_subdirectory(Dynamics)
add_subdirectory(Control)
add_subdirectory(Evolution)
add_subdirectory(Solver)

find_package(Eigen3 REQUIRED)

set_target_properties(cmechanism PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_link_libraries(cmechanism PRIVATE Eigen3::Eigen)

# If possible try to build with OpenMP enabled to take advantage of some multithreading functions in Eigen
#find_package(OpenMP)
#if(OpenMP_CXX_FOUND)
#    target_link_libraries(cmechanism PRIVATE OpenMP::OpenMP_CXX)
#endif()

if(WITH_TESTS)
    # Enable the testing suite, currently just a mechanism simulation trial
    add_executable(mechanismtest main.cpp)
    target_link_libraries(mechanismtest PUBLIC cmechanism)
endif()

