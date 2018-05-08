cmake_minimum_required(VERSION 3.0.0)
project(hankelhacker VERSION 0.0.0)
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_VERBOSE_MAKEFILE on)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules")
set (CMAKE_CXX_STANDARD 17)

# Globals
set(GLOBAL_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/install)

# for gcc
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${GLOBAL_OUTPUT_PATH})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_OUTPUT_PATH})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${GLOBAL_OUTPUT_PATH})

# for other compilers
foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${GLOBAL_OUTPUT_PATH})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${GLOBAL_OUTPUT_PATH})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${GLOBAL_OUTPUT_PATH})
endforeach(OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES)

# enable OpenMP
option (USE_OPENMP "Use OpenMP" ON)
if(USE_OPENMP)
    #find_package(OpenMP)
    #if(OPENMP_FOUND)
    #    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    #    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    #endif()
    set(OPENMP_FLAG "-fopenmp=libomp")
    check_cxx_compiler_flag(-fopenmp HAS_OPENMP_C_FLAG)
    check_cxx_compiler_flag(-fopenmp HAS_OPENMP_CXX_FLAG)
    if (HAS_OPENMP_C_FLAG AND HAS_OPENMP_CXX_FLAG)
        message(STATUS "OpenMP enabled")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPENMP_FLAG}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMP_FLAG}")
    endif (HAS_OPENMP_C_FLAG AND HAS_OPENMP_CXX_FLAG)
endif(USE_OPENMP)

include(CTest)
enable_testing()

find_package(GMP REQUIRED)
find_package(GSL REQUIRED)
file(GLOB SOURCES src/*.cpp)
file(GLOB HEADERS src/*.hpp)
add_executable(hankelhacker ${SOURCES})
target_link_libraries(hankelhacker gmp gmpxx gsl gslcblas)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

set(CONFIGURED_ONCE TRUE CACHE INTERNAL "Flag showing whether CMake has already configured once.")