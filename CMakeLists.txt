cmake_minimum_required(VERSION 3.13)
## Build the python bindings

project(pymechlib)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	# using Clang
	set(CMAKE_CXX_FLAGS "-march=native -O3 -fcxx-exceptions")
	message(STATUS ">>> Building with Clang")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	# Using GCC
    set(CMAKE_CXX_FLAGS "-march=native -O3")
	message(STATUS ">>> Building with GCC")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# using Visual Studio C++
  	set(CMAKE_CXX_FLAGS "/GL /O2 /EHsc /fp:precise")
	message(STATUS ">>> Building with MSVC")
endif()

set(SOURCE_DIR "src/pymechlib")
add_subdirectory(lib/mechanismlib)

find_package(pybind11 REQUIRED)

pybind11_add_module(cmechlib ${SOURCE_DIR}/bindings.cpp ${SOURCE_DIR}/bindings.h ${SOURCE_DIR}/threader.cpp ${SOURCE_DIR}/threader.h)

target_include_directories(cmechlib PRIVATE ${SOURCE_DIR} "lib/mechanismlib")
target_link_libraries(cmechlib PRIVATE cmechanism)
