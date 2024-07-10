cmake_minimum_required(VERSION 3.12...3.27)

project(fastprimes
  DESCRIPTION "Fast small primes generator"
  LANGUAGES CXX
)

include(AICxxProject)

#==============================================================================
# OPTIONS

# Set a default value for OptionFastPrimesCompression.
set(OptionFastPrimesCompression "6" CACHE STRING "Compession level of Primes (3 - 6). Default: 6")

message(DEBUG "OptionFastPrimesCompression is ${OptionFastPrimesCompression}")
set(CW_FASTPRIMES_COMPRESSION ${OptionFastPrimesCompression})

#==============================================================================
# PLATFORM SPECIFIC CHECKS
#

# Specify fastprimes specific configure file.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/config.h
               @ONLY)

#==============================================================================
# BUILD OBJECT LIBRARY
#

# This project is an OBJECT-library, used by other git submodules and the main project.
add_library(fastprimes_ObjLib OBJECT)

# Tell cwds/sys.h that we have a config.h.
target_compile_definitions(fastprimes_ObjLib
  PUBLIC
    HAVE_FASTPRIMES_CONFIG_H
)

# The list of source files.
target_sources(fastprimes_ObjLib
  PRIVATE
    "Primes.cxx"
    "Primes.h"
)

# Required include search-paths.
get_target_property(CWDS_INTERFACE_INCLUDE_DIRECTORIES AICxx::cwds INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(fastprimes_ObjLib
  PUBLIC
    "${CWDS_INTERFACE_INCLUDE_DIRECTORIES}" # For sys.h and debug.h.
    "${top_objdir}"                         # For fastprimes/config.h.
)

# Require support for C++20.
target_compile_features(fastprimes_ObjLib
  PUBLIC cxx_std_20
)

# Set link dependencies.
target_link_libraries(fastprimes_ObjLib
  PUBLIC
    AICxx::utils
)

# Create an ALIAS target.
add_library(AICxx::fastprimes ALIAS fastprimes_ObjLib)

# Prepend this object library to the list.
set(AICXX_OBJECTS_LIST AICxx::fastprimes ${AICXX_OBJECTS_LIST} CACHE INTERNAL "List of OBJECT libaries that this project uses.")
