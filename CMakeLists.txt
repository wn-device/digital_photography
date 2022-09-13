# Use of SOURCE_DIR target property begins at cmake 3.7
if(APPLE)
    cmake_minimum_required(VERSION 3.7)
else()
    cmake_minimum_required(VERSION 3.5)
endif()

# Compiler id for Apple Clang
if(APPLE)
    cmake_policy(SET CMP0025 NEW)
    cmake_policy(SET CMP0037 NEW)
endif()

project (RawTherapee)
# The default target is Debug:
if(CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING "One of: None Debug Release RelWithDebInfo MinSizeRel"
              FORCE)
endif()

string(TOUPPER ${CMAKE_BUILD_TYPE} UPPER_CMAKE_BUILD_TYPE)

# Set required C and C++ standards and check GCC version:
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION
                                            VERSION_LESS "4.9")
    message(
        FATAL_ERROR
            "Building RawTherapee requires using GCC version 4.9 or higher!")
endif()

# Warning for GCC vectorization issues, which causes problems #5749 and #6384:
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND ((CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "10.0" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "10.2") OR (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "11.0")))
    message(STATUS "WARNING: gcc ${CMAKE_CXX_COMPILER_VERSION} is known to miscompile RawTherapee when using -ftree-loop-vectorize, forcing the option to be off")    
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-tree-loop-vectorize")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-tree-loop-vectorize")
endif()


message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

# Cache name suffix examples: "" = ~/.config/RawTherapee, "5" =
# ~/.config/RawTherapee-5, "_testing" = ~/.config/RawTherapee_testing Use "" for
# stable releases and "5-dev" for anything else.
set(CACHE_NAME_SUFFIX "" CACHE STRING "RawTherapee's cache folder suffix")
find_package(PkgConfig REQUIRED)
pkg_check_modules(GLIB2 REQUIRED glib-2.0>=2.48)
pkg_check_modules(GLIBMM REQUIRED glibmm-2.4>=2.48)
set(SOURCES
    main.cc
    hilite_recon.cc
    fread_dump.cc
)

include_directories(${GLIBMM_INCLUDE_DIRS})
link_directories(. ${GLIBMM_LIBRARY_DIRS})
add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${GLIBMM_LIBRARY_DIRS} ${PROJECT_SOURCE_DIR}) 
target_link_libraries(${PROJECT_NAME}  PUBLIC  ${GLIBMM_LIBRARY_DIRS} ${PROJECT_SOURCE_DIR})