cmake_minimum_required(VERSION 3.11)
project(CodingExercises)

if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND git rev-parse HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH_SHORT
    OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  set(GIT_BRANCH "")
  set(GIT_COMMIT_HASH "")
  set(GIT_COMMIT_HASH_SHORT "")
endif()

message(STATUS "On Git Branch: ${GIT_BRANCH} (${GIT_COMMIT_HASH})")

set(DEFAULT_BUILD_TYPE "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${DEFAULT_BUILD_TYPE}"
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui, ccmake
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()
set(CMAKE_TRY_COMPILE_CONFIGURATION ${CMAKE_BUILD_TYPE})
set(CMAKE_CONFIGURATION_TYPES
    "Debug;Release"
    CACHE STRING "" FORCE)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMake")
message(STATUS "CMake module path: ${CMAKE_MODULE_PATH}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

# Add math defines for Windows
if(WIN32)
    add_definitions(-D_USE_MATH_DEFINES)
endif()

set(MSVC_COMPILER_FLAGS
    $<$<COMPILE_LANGUAGE:CXX>:/MP> # Use as many cores as available for
                                   # compilation.
    $<$<COMPILE_LANGUAGE:CXX>:/W4> # /Wall # Baseline reasonable warnings.
    $<$<COMPILE_LANGUAGE:CXX>:/permissive-> # Standards conformance mode for
                                            # MSVC compiler.
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Release>:/Oi>>
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Release>:/O2>>
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Release>:/Ot>>
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Release>:/GL>>
    $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Release>:/GF>>)

set(CLANG_COMPILER_FLAGS
    -g
    -W
    -Wall
    -Wextra # Reasonable and standard.
    -Wshadow # Warn the user if a variable declaration shadows one from a parent
             # context.
    -Wnon-virtual-dtor # Warn the user if a class with virtual functions has a
                       # non-virtual destructor. This helps catch hard to track
                       # down memory errors.
    -Wcast-align # Warn for potential performance problem casts.
    -Wunused # Warn on anything being unused.
    -Wunused-function
    -Wunused-variable
    -Wunused-parameter
    -Woverloaded-virtual # Warn if you overload (not override) a virtual
                         # function.
    $<$<CONFIG:Debug>:-O0;-ggdb;-fno-elide-constructors>)

set(GCC_COMPILER_FLAGS
    ${CLANG_COMPILER_FLAGS}
)

if(MSVC)
  set(COMPILER_FLAGS ${MSVC_COMPILER_FLAGS})
elseif(CMAKE_CXX_COMPILER_ID MATCHES ".*Clang")
  set(COMPILER_FLAGS ${CLANG_COMPILER_FLAGS})
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(COMPILER_FLAGS ${GCC_COMPILER_FLAGS})
else()
  message(
    AUTHOR_WARNING
      "No compiler flags set for '${CMAKE_CXX_COMPILER_ID}' compiler.")
endif()

file(WRITE "${PROJECT_BINARY_DIR}/.gitignore" "*")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
message(STATUS "CMake module path: ${CMAKE_MODULE_PATH}")

add_subdirectory(thirdparty)

# Blade Element Theory Project
set(BET_SOURCES
    src/main.cpp
    src/constants.hpp
    src/bladeelement.hpp
    src/blade.hpp
    src/rotor.hpp
    src/vtksettings.hpp
)

add_executable(BET-VPM ${BET_SOURCES})
target_compile_options(BET-VPM PRIVATE ${COMPILER_FLAGS})
target_link_libraries(BET-VPM PRIVATE Eigen3::Eigen MatplotlibCpp::Matplotlib)
target_include_directories(BET-VPM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Other executables

add_executable(denemeKernel src/deneme.cpp)
target_compile_options(denemeKernel PRIVATE ${COMPILER_FLAGS}) 
target_link_libraries(denemeKernel PRIVATE Eigen3::Eigen MatplotlibCpp::Matplotlib)

add_executable(DENEME src/denemeCode.cpp)
target_compile_options(DENEME PRIVATE ${COMPILER_FLAGS}) 
target_link_libraries(DENEME PRIVATE Eigen3::Eigen MatplotlibCpp::Matplotlib)