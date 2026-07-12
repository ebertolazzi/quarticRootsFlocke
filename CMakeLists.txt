############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita degli Studi di Trento                                    #
#      email: enrico.bertolazzi@unitn.it                                   #
#                                                                          #
############################################################################

# --------------------------------------------------------------------------
# Minimum CMake version
# --------------------------------------------------------------------------
# This project uses a modern target-based CMake layout, generator expressions,
# install/export rules, and policy handling.  Version 3.16 is kept as the lower
# bound to preserve compatibility with older systems while still supporting the
# features used below.
cmake_minimum_required(VERSION 3.16)

# --------------------------------------------------------------------------
# CMake policy selection
# --------------------------------------------------------------------------
# Policies are set explicitly so the build behaves predictably across CMake
# releases.  Each policy is guarded because older CMake versions do not know
# about policies introduced later.
#
# CMP0177 NEW, when available, enables the newer behavior for install()
# destination normalization introduced by recent CMake versions.
# CMP0135 NEW, when available, makes timestamp handling for extracted archives
# more reproducible.
# CMP0169 OLD is kept for compatibility with legacy FetchContent behavior used
# by related projects in this code base.
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.31)
  cmake_policy(SET CMP0177 NEW)
endif()
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()
if(POLICY CMP0169)
  cmake_policy(SET CMP0169 OLD)
endif()

# --------------------------------------------------------------------------
# Apple FAT/universal build setup
# --------------------------------------------------------------------------
# On Apple platforms a single library can contain slices for several CPU
# architectures, for example arm64 and x86_64.  This is controlled by
# CMAKE_OSX_ARCHITECTURES.
#
# The variable must be set before project(), otherwise CMake may initialize
# compilers and targets before the requested architectures are known.
#
# Default behavior:
#   - ON  on Apple hosts or Apple cross-compilation targets.
#   - OFF on Linux/Windows and other platforms, where there is no equivalent
#     single-file universal library format.
#
# User overrides:
#   -DUTILS_BUILD_FAT=OFF
#   -DUTILS_FAT_ARCHITECTURES=arm64
#   -DCMAKE_OSX_ARCHITECTURES=arm64\;x86_64
set(UTILS_BUILD_FAT_DEFAULT OFF)
if(CMAKE_HOST_APPLE OR CMAKE_SYSTEM_NAME MATCHES "Darwin|iOS|tvOS|watchOS|visionOS")
  set(UTILS_BUILD_FAT_DEFAULT ON)
endif()
option(UTILS_BUILD_FAT "Build Apple universal/fat binaries when supported" ${UTILS_BUILD_FAT_DEFAULT})
set(UTILS_FAT_ARCHITECTURES "arm64;x86_64" CACHE STRING "Apple architectures to build when UTILS_BUILD_FAT is ON")

if(UTILS_BUILD_FAT AND NOT CMAKE_OSX_ARCHITECTURES)
  set(
    CMAKE_OSX_ARCHITECTURES
    "${UTILS_FAT_ARCHITECTURES}"
    CACHE STRING "Build architectures for Apple platforms"
    FORCE
  )
endif()

# --------------------------------------------------------------------------
# User-configurable build options
# --------------------------------------------------------------------------
# These cache options define the public configuration surface of the project.
# They can be changed from the command line with -D<name>=<value>, from ccmake,
# from cmake-gui, or by a parent project embedding this project with
# add_subdirectory().
#
# UTILS_GIT_VERSION controls whether Git tags are used for the project version.
# UTILS_ENABLE_TESTS controls the test executable section near the end.
# UTILS_BUILD_SHARED adds a shared-library target next to the always-built
# static library.
# UTILS_INSTALL_LOCAL enables the local install layout used by the historical
# project structure.
# UTILS_TWEAK_VERSION_IS_TIMESTAMP changes the fallback version suffix from 0
# to a compact timestamp.
# UTILS_ENABLE_ARBITRARY_PRECISION enables the optional Boost.Multiprecision
# backend when the resolver module can find a usable Boost include tree.
option(UTILS_GIT_VERSION "Get version from GIT command" ON)
option(UTILS_ENABLE_TESTS "Enable compilation of testing targets" ON)
option(UTILS_BUILD_SHARED "Enable shared option compilation" OFF)
option(UTILS_INSTALL_LOCAL "Enable installation on local directory" ON)
option(UTILS_TWEAK_VERSION_IS_TIMESTAMP "Sets tweak version to timestamp" OFF)
option(UTILS_ENABLE_ARBITRARY_PRECISION "Enable Boost-based arbitrary precision HQ solvers when Boost is available" OFF)

# --------------------------------------------------------------------------
# Compatibility namespace and install target accumulator
# --------------------------------------------------------------------------
# The historical cmake_utils-based projects exported targets below a configurable
# namespace.  Keep that namespace cache-visible so external scripts can preserve
# the previous target names, for example utils::quarticRootsFlocke_Static.
#
# TARGETS_TO_BE_INSTALLED is populated as targets are created and is consumed by
# the single install(TARGETS ...) call near the bottom of the file.
set(UTILS_NAMESPACE "utils" CACHE STRING "Sets the namespace used for exports")
set(TARGETS_TO_BE_INSTALLED)

# --------------------------------------------------------------------------
# Global CMake behavior
# --------------------------------------------------------------------------
# These settings affect the whole directory scope.  They intentionally mirror
# the conventions used by the previous build infrastructure: quiet install
# output, compile_commands.json generation, non-verbose Makefiles by default,
# strict C++ standard selection, and an install RPATH suitable for libraries
# loaded from the same installed directory.
set(CMAKE_INSTALL_MESSAGE NEVER)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_VERBOSE_MAKEFILE FALSE)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_INSTALL_RPATH $ORIGIN)
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)

# On Windows, export all symbols automatically for shared libraries.  This keeps
# the project usable without manually annotating every public API with
# __declspec(dllexport).  On Unix-like systems this variable has no effect and is
# therefore not enabled.
if(NOT UNIX)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Make mathematical constants such as M_PI visible on MSVC-compatible standard
# library implementations.  The definition is harmless on compilers where those
# constants are already available or not controlled by this macro.
add_compile_definitions(_USE_MATH_DEFINES)

# --------------------------------------------------------------------------
# Project version resolution
# --------------------------------------------------------------------------
# The preferred version source is the most recent Git tag.  If Git is disabled,
# unavailable, or the source tree is not a Git checkout, the build falls back to
# UTILS_CURRENT_VERSION plus a tweak component.  The tweak component is normally
# 0, but can be replaced by a timestamp for local snapshot builds.
set(
  UTILS_CURRENT_VERSION "1.0"
  CACHE STRING "Default project version, if not taken from git"
)

# Look for Git quietly.  Failure is not an error because exported archives and
# source distributions often do not contain a .git directory.
find_package(Git QUIET)
if(UTILS_GIT_VERSION AND Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE UTILS_PROJECT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
endif()

if(NOT UTILS_PROJECT_VERSION)
  set(UTILS_BUILD_TIMESTAMP "0")
  if(UTILS_TWEAK_VERSION_IS_TIMESTAMP)
    string(TIMESTAMP UTILS_BUILD_TIMESTAMP "%y%m%d%H")
  endif()
  set(UTILS_PROJECT_VERSION "${UTILS_CURRENT_VERSION}.${UTILS_BUILD_TIMESTAMP}")
endif()

if(UTILS_PROJECT_VERSION MATCHES "^v")
  string(SUBSTRING "${UTILS_PROJECT_VERSION}" 1 -1 UTILS_PROJECT_VERSION)
endif()

# --------------------------------------------------------------------------
# Project declaration
# --------------------------------------------------------------------------
# project() enables the requested languages, initializes compiler variables,
# computes PROJECT_* metadata, and makes PROJECT_IS_TOP_LEVEL available.  The
# Apple universal-build settings above must already be present before this call.
project(
  quarticRootsFlocke
  VERSION ${UTILS_PROJECT_VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/quarticRootsFlocke/"
  LANGUAGES CXX
)

# --------------------------------------------------------------------------
# Standard CMake modules
# --------------------------------------------------------------------------
# GNUInstallDirs defines portable install directory variables such as
# CMAKE_INSTALL_INCLUDEDIR, CMAKE_INSTALL_LIBDIR, and CMAKE_INSTALL_BINDIR.
# CTest provides BUILD_TESTING plus the enable_testing()/add_test() machinery.
# ResolveBoostMultiprecision.cmake is a project-local helper that locates the
# optional header-only Boost.Multiprecision backend.
include(GNUInstallDirs)
include(CTest)
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/ResolveBoostMultiprecision.cmake")

# --------------------------------------------------------------------------
# Local install prefix
# --------------------------------------------------------------------------
# When this repository is configured as the top-level project, default installs
# are redirected into the source tree instead of /usr/local.  This preserves the
# historical layout where headers and libraries are staged under ./lib.
#
# If a parent project embeds this repository, or if the user explicitly provides
# CMAKE_INSTALL_PREFIX, the prefix is not overwritten.
if(PROJECT_IS_TOP_LEVEL AND (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR CMAKE_INSTALL_PREFIX STREQUAL "/usr/local"))
  set(
    CMAKE_INSTALL_PREFIX
    "${CMAKE_CURRENT_SOURCE_DIR}"
    CACHE PATH "Install path prefix, prepended onto install directories."
    FORCE
  )
endif()

# For standalone builds, force the install subdirectories to the local project
# layout and create the directories eagerly.  This ensures that post-build copy
# and symlink operations have a predictable destination tree.
if(PROJECT_IS_TOP_LEVEL)
  set(CMAKE_INSTALL_INCLUDEDIR "lib/include" CACHE PATH "Install include directory" FORCE)
  set(CMAKE_INSTALL_LIBDIR "lib/lib" CACHE PATH "Install library directory" FORCE)
  set(CMAKE_INSTALL_BINDIR "lib/bin" CACHE PATH "Install binary directory" FORCE)
  file(MAKE_DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/lib"
    "${CMAKE_CURRENT_SOURCE_DIR}/lib/bin"
  )
endif()

# --------------------------------------------------------------------------
# Default build type
# --------------------------------------------------------------------------
# Single-configuration generators such as Unix Makefiles and Ninja use
# CMAKE_BUILD_TYPE.  Multi-configuration generators such as Visual Studio and
# Xcode use CMAKE_CONFIGURATION_TYPES and ignore CMAKE_BUILD_TYPE.  Default to
# Release only when neither mechanism has already been selected.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# --------------------------------------------------------------------------
# Warning and diagnostic policy
# --------------------------------------------------------------------------
# COMMON_CXX_FLAGS contains warnings shared by GCC and Clang.  They favor strict
# diagnostics for numerical code: conversions, shadowing, virtual destructors,
# formatting issues, and suspicious control flow are all reported.
set(COMMON_CXX_FLAGS
  "-pedantic"
  "-Wall"
  "-Wextra"
  "-Wshadow"
  "-Wnon-virtual-dtor"
  "-Wold-style-cast"
  "-Wcast-align"
  "-Wunused"
  "-Woverloaded-virtual"
  "-Wpedantic"
  "-Wconversion"
  "-Wsign-conversion"
  "-Wdouble-promotion"
  "-Wformat=2"
  "-Wimplicit-fallthrough"
  "-Wreturn-local-addr"
)

# GCC exposes several additional diagnostics that are either unavailable or have
# different names on Clang.  Keep them guarded to avoid passing unsupported flags
# to other compilers.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND COMMON_CXX_FLAGS
    "-Wmisleading-indentation"
    "-Wduplicated-cond"
    "-Wduplicated-branches"
    "-Wlogical-op"
    "-Wnull-dereference"
    "-Wuseless-cast"
  )
endif()

string(JOIN " " COMMON_CXX_FLAGS_STR ${COMMON_CXX_FLAGS})

# MSVC warning configuration.  These flags approximate the strictness of the
# GCC/Clang warning set while using MSVC's numeric warning identifiers and the
# /external options to reduce noise from third-party or system headers.
set(MSVC_COMPILE_FLAGS
  "/permissive"
  "/W4"
  "/w14242"
  "/w14254"
  "/w14263"
  "/w14265"
  "/w14287"
  "/we4289"
  "/w14296"
  "/w14311"
  "/w14545"
  "/w14546"
  "/w14547"
  "/w14549"
  "/w14555"
  "/w14619"
  "/w14640"
  "/w14826"
  "/w14905"
  "/w14906"
  "/w14928"
  "/external:W0"
  "/external:anglebrackets"
  "/external:templates-compat"
)

string(JOIN " " MSVC_COMPILE_FLAGS_STR ${MSVC_COMPILE_FLAGS})

# Visual Studio-specific preprocessor definitions and baseline flags.  They are
# separated from MSVC_COMPILE_FLAGS because they are applied to the global C and
# C++ configuration strings for the selected build type.
set(VS_DEFINE "/D_WINDOWS /D_SCL_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY /DYAML_DECLARE_STATIC /DPCRE_STATIC")
set(VS_FLAGS "/nologo /std:c++17 /utf-8 /GS /Gy /GF /EHsc /bigobj")

# --------------------------------------------------------------------------
# Compiler-specific global flags
# --------------------------------------------------------------------------
# This project still uses global CMAKE_<LANG>_FLAGS variables rather than only
# target_compile_options().  The settings below select the C++ dialect,
# optimization/debug behavior, position-independent code, and compiler-specific
# standard library options.
#
# The target-level code later also requests cxx_std_20.  The explicit -std flags
# remain here for compatibility with the existing build conventions.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "-std=c++20 ${COMMON_CXX_FLAGS_STR}")
  set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -O2 -g ")
  set(CMAKE_CXX_FLAGS_DEBUG "-fPIC -O0 -ggdb ")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++20 -stdlib=libc++ ${COMMON_CXX_FLAGS_STR}")
  set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -Os -g -fcolor-diagnostics ")
  set(CMAKE_CXX_FLAGS_DEBUG "-fPIC -O0 -gfull -fcolor-diagnostics ")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  set(CMAKE_CXX_FLAGS " ${MSVC_COMPILE_FLAGS_STR} ${VS_DEFINE} /sdl")
  set(CMAKE_CXX_FLAGS_RELEASE " ${VS_FLAGS} /O2 /MD /Zo ")
  set(CMAKE_CXX_FLAGS_DEBUG " ${VS_FLAGS} /Od /Ob0 /MDd /guard:cf /Zi /DEBUG:FULL ")
else()
  message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
endif()

# Mirror the C++ release/debug flag choice into the C language.  The project is
# primarily C++, but keeping C flags coherent avoids surprises if a future source
# or dependency introduces C compilation in the same directory scope.
set(CMAKE_C_FLAGS "")
set(CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
set(CMAKE_C_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --------------------------------------------------------------------------
# Optional local customization hook
# --------------------------------------------------------------------------
# A developer can place CMakeLists-customize.txt one directory above this source
# tree to inject site-specific flags or paths.  The file is intentionally not
# required and should not be needed for normal builds.
if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/../CMakeLists-customize.txt")
  include("${CMAKE_CURRENT_LIST_DIR}/../CMakeLists-customize.txt")
endif()

# --------------------------------------------------------------------------
# Platform suffix used in artifact names
# --------------------------------------------------------------------------
# The real library file name contains a platform suffix, for example
# libquarticRootsFlocke_osx_static.a or libquarticRootsFlocke_linux_static.a.
# A generic unsuffixed compatibility alias is created later.
#
# On Windows, distinguish MinGW from the native toolchain and encode pointer
# size in the suffix.  On unknown platforms, fall back to a lowercase CMake
# system name.
if(APPLE)
  set(UTILS_PLATFORM_SUFFIX osx)
elseif(WIN32)
  if(MINGW)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
      set(UTILS_PLATFORM_SUFFIX mingw64)
    else()
      set(UTILS_PLATFORM_SUFFIX mingw32)
    endif()
  elseif(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(UTILS_PLATFORM_SUFFIX win64)
  else()
    set(UTILS_PLATFORM_SUFFIX win32)
  endif()
elseif(UNIX)
  set(UTILS_PLATFORM_SUFFIX linux)
else()
  string(TOLOWER "${CMAKE_SYSTEM_NAME}" UTILS_PLATFORM_SUFFIX)
endif()

# --------------------------------------------------------------------------
# Source tree discovery
# --------------------------------------------------------------------------
# SOURCE_DIR is the root of the public headers and implementation files.  The
# project uses globbing with CONFIGURE_DEPENDS to keep the file list updated when
# sources are added or removed and CMake is re-run by the generator.
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Collect implementation units.  Only .cc files are compiled into the library;
# headers are listed separately so IDE generators can display them in the target.
file(
  GLOB SOURCES
  CONFIGURE_DEPENDS
  "${SOURCE_DIR}/*.cc"
)

# Collect public headers recursively.  They are installed preserving their
# relative path below src/, and they are marked HEADER_FILE_ONLY so IDEs show
# them without trying to compile them.
file(
  GLOB_RECURSE PUBLIC_HEADERS
  CONFIGURE_DEPENDS
  "${SOURCE_DIR}/*.h"
  "${SOURCE_DIR}/*.hh"
  "${SOURCE_DIR}/*.hpp"
  "${SOURCE_DIR}/*.hxx"
)

# Mark headers as non-compilable sources.  This is mainly useful for IDE project
# generators, where adding headers to a target improves navigation without
# changing the build graph.
foreach(H IN LISTS PUBLIC_HEADERS)
  set_source_files_properties("${H}" PROPERTIES HEADER_FILE_ONLY TRUE)
endforeach()

# The source directory is both the build-tree public include root and the private
# include root used by the project implementation.
set(INCLUDE_DIRS "${SOURCE_DIR}")

# --------------------------------------------------------------------------
# Optional Boost.Multiprecision backend
# --------------------------------------------------------------------------
# The default build is double-only.  When UTILS_ENABLE_ARBITRARY_PRECISION is ON,
# a project-local resolver tries to find Boost.Multiprecision.  The public macro
# POLYNOMIAL_ROOTS_HAS_MULTIPRECISION is then set per target so client code and
# tests can compile the correct code path.
set(QRF_ENABLE_ARBITRARY_PRECISION OFF)
set(QRF_BOOST_INCLUDE_DIR "")
set(QRF_BOOST_SOURCE "disabled")
if(UTILS_ENABLE_ARBITRARY_PRECISION)
  qrf_resolve_boost_multiprecision(QRF_BOOST_INCLUDE_DIR QRF_BOOST_SOURCE)
  if(QRF_BOOST_INCLUDE_DIR)
    set(QRF_ENABLE_ARBITRARY_PRECISION ON)
  else()
    message(WARNING "UTILS_ENABLE_ARBITRARY_PRECISION=ON but Boost multiprecision was not found. Falling back to double-only build.")
  endif()
endif()

# --------------------------------------------------------------------------
# qrf_configure_target(target_name)
# --------------------------------------------------------------------------
# Applies the common target-level configuration to each library variant.  Keeping
# this in one function prevents the static and shared targets from diverging.
#
# Responsibilities:
#   - request C++20 at target level;
#   - add the optional Boost include directory when multiprecision is enabled;
#   - publish POLYNOMIAL_ROOTS_HAS_MULTIPRECISION to consumers;
#   - expose the source tree in build mode and the install include directory in
#     install mode;
#   - set version metadata, Xcode header search paths, output directories, and
#     position-independent code.
function(qrf_configure_target target_name)
  # Require C++20 for consumers as well as for this target.  Public scope is
  # appropriate because installed headers may expose C++20 constructs.
  target_compile_features(${target_name} PUBLIC cxx_std_20)

  # Boost.Multiprecision is header-only for this use case.  Marking the include
  # directory SYSTEM suppresses third-party warnings while compiling this target.
  if(QRF_ENABLE_ARBITRARY_PRECISION)
    target_include_directories(${target_name} SYSTEM PRIVATE "${QRF_BOOST_INCLUDE_DIR}")
  endif()
  # Publish the multiprecision state as a numeric preprocessor macro.  The
  # generator expression expands to 0 or 1 for each build configuration.
  target_compile_definitions(
    ${target_name}
    PUBLIC
      POLYNOMIAL_ROOTS_HAS_MULTIPRECISION=$<BOOL:${QRF_ENABLE_ARBITRARY_PRECISION}>
  )
  # Use generator expressions to expose different include roots in the build
  # tree and in the install tree.  PRIVATE includes are only needed while
  # compiling this target.
  target_include_directories(
    ${target_name}
    PUBLIC
      "$<BUILD_INTERFACE:${SOURCE_DIR}>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    PRIVATE
      ${INCLUDE_DIRS}
  )

  set_target_properties(
    ${target_name}
    PROPERTIES
      VERSION ${PROJECT_VERSION}
      XCODE_ATTRIBUTE_HEADER_SEARCH_PATHS "${SOURCE_DIR}"
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib/lib"
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib/lib"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib/bin"
      POSITION_INDEPENDENT_CODE ON
  )
endfunction()

# --------------------------------------------------------------------------
# qrf_add_platform_library_alias(target_name, base_name, library_kind)
# --------------------------------------------------------------------------
# Gives the real library a platform-specific file name and creates an unsuffixed
# compatibility alias next to it.
#
# Example for a static macOS build:
#   real file : libquarticRootsFlocke_osx_static.a
#   alias     : libquarticRootsFlocke.a -> libquarticRootsFlocke_osx_static.a
#
# On Windows the alias is a copied file because symlink behavior depends on user
# privileges and filesystem settings.  On Unix-like systems the alias is a
# symbolic link.
#
# The same aliasing logic is repeated at install time so the installed tree has
# the same compatibility layout as the build tree.
function(qrf_add_platform_library_alias target_name base_name library_kind)
  # The real library keeps the platform/static suffix, e.g.
  #   libquarticRootsFlocke_osx_static.a
  # and the unsuffixed name is kept only as a compatibility alias, e.g.
  #   libquarticRootsFlocke.a -> libquarticRootsFlocke_osx_static.a
  if(library_kind STREQUAL "STATIC")
    set(output_name "${base_name}_${UTILS_PLATFORM_SUFFIX}_static")
    set(alias_filename "${CMAKE_STATIC_LIBRARY_PREFIX}${base_name}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if(WIN32)
      set(output_name "lib${output_name}")
      set(real_filename "${output_name}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    else()
      set(real_filename "${CMAKE_STATIC_LIBRARY_PREFIX}${output_name}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()
  else()
    set(output_name "${base_name}_${UTILS_PLATFORM_SUFFIX}")
    set(alias_filename "${CMAKE_SHARED_LIBRARY_PREFIX}${base_name}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    if(WIN32)
      set(output_name "lib${output_name}")
      set(real_filename "${output_name}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    else()
      set(real_filename "${CMAKE_SHARED_LIBRARY_PREFIX}${output_name}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
  endif()

  # Rename the physical artifact.  This affects the file written to disk, not
  # the logical CMake target name used by target_link_libraries().
  set_target_properties(${target_name} PROPERTIES OUTPUT_NAME "${output_name}")

  if(WIN32)
    # Windows compatibility alias: copy the produced artifact to the generic
    # file name after every successful build.
    add_custom_command(
      TARGET ${target_name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              "$<TARGET_FILE:${target_name}>"
              "$<TARGET_FILE_DIR:${target_name}>/${alias_filename}"
      VERBATIM
    )

    if(UTILS_INSTALL_LOCAL)
      string(CONCAT _alias_install_code
        "file(MAKE_DIRECTORY \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")\n"
        "execute_process(\n"
        "  COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different\n"
        "          \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${real_filename}\"\n"
        "          \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${alias_filename}\"\n"
        ")\n"
      )
      install(CODE "${_alias_install_code}")
    endif()
  else()
    # Unix compatibility alias: recreate the symbolic link after every build
    # so it always points to the newest platform-suffixed artifact.
    add_custom_command(
      TARGET ${target_name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E rm -f
              "$<TARGET_FILE_DIR:${target_name}>/${alias_filename}"
      COMMAND ${CMAKE_COMMAND} -E create_symlink
              "$<TARGET_FILE_NAME:${target_name}>"
              "$<TARGET_FILE_DIR:${target_name}>/${alias_filename}"
      VERBATIM
    )

    if(UTILS_INSTALL_LOCAL)
      string(CONCAT _alias_install_code
        "file(MAKE_DIRECTORY \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")\n"
        "execute_process(\n"
        "  COMMAND \"${CMAKE_COMMAND}\" -E rm -f\n"
        "          \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${alias_filename}\"\n"
        ")\n"
        "execute_process(\n"
        "  COMMAND \"${CMAKE_COMMAND}\" -E create_symlink\n"
        "          \"${real_filename}\"\n"
        "          \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${alias_filename}\"\n"
        ")\n"
      )
      install(CODE "${_alias_install_code}")
    endif()
  endif()
endfunction()

# --------------------------------------------------------------------------
# qrf_final_messages()
# --------------------------------------------------------------------------
# Prints a verbose configuration summary at the end of configure.  This is
# intentionally more detailed than a typical CMake project because it records
# compiler IDs, flags, platform data, discovered sources, install paths, and the
# optional Boost backend state.  It is useful for reproducing build issues across
# Linux, macOS, Windows, and cross-compilation environments.
function(qrf_final_messages)
  message(STATUS "C compiler                  = ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER}")
  message(STATUS "C++ compiler                = ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER}")
  message(STATUS "CMAKE_CXX_FLAGS_RELEASE     = ${CMAKE_CXX_FLAGS_RELEASE}")
  message(STATUS "CMAKE_CXX_FLAGS             = ${CMAKE_CXX_FLAGS}")
  message(STATUS "CMAKE_C_FLAGS_RELEASE       = ${CMAKE_C_FLAGS_RELEASE}")
  message(STATUS "CMAKE_C_FLAGS               = ${CMAKE_C_FLAGS}")
  message(STATUS "CMAKE_SYSTEM_NAME           = ${CMAKE_SYSTEM_NAME}")
  message(STATUS "CMAKE_SYSTEM_PROCESSOR      = ${CMAKE_SYSTEM_PROCESSOR}")
  if(APPLE)
    message(STATUS "CMAKE_OSX_SYSROOT           = ${CMAKE_OSX_SYSROOT}")
    message(STATUS "CMAKE_OSX_ARCHITECTURES     = ${CMAKE_OSX_ARCHITECTURES}")
    message(STATUS "CMAKE_OSX_DEPLOYMENT_TARGET = ${CMAKE_OSX_DEPLOYMENT_TARGET}")
  endif()
  message(STATUS "SOURCES = ${SOURCES}")
  message(STATUS "HEADERS = ${PUBLIC_HEADERS}")
  message(STATUS "TARGETS_TO_BE_INSTALLED  = ${TARGETS_TO_BE_INSTALLED}")
  message(STATUS "CMAKE_HOME_DIRECTORY     = ${CMAKE_HOME_DIRECTORY}")
  message(STATUS "CMAKE_CURRENT_SOURCE_DIR = ${CMAKE_CURRENT_SOURCE_DIR}")
  message(STATUS "CMAKE_INSTALL_PREFIX     = ${CMAKE_INSTALL_PREFIX}")
  message(STATUS "UTILS_GIT_VERSION   = ${UTILS_GIT_VERSION}")
  message(STATUS "UTILS_ENABLE_TESTS  = ${UTILS_ENABLE_TESTS}")
  message(STATUS "UTILS_BUILD_SHARED  = ${UTILS_BUILD_SHARED}")
  message(STATUS "UTILS_ENABLE_ARBITRARY_PRECISION = ${UTILS_ENABLE_ARBITRARY_PRECISION}")
  message(STATUS "QRF_ENABLE_ARBITRARY_PRECISION   = ${QRF_ENABLE_ARBITRARY_PRECISION}")
  message(STATUS "UTILS_INSTALL_LOCAL = ${UTILS_INSTALL_LOCAL}")
  message(STATUS "Headers install     = ${CMAKE_INSTALL_INCLUDEDIR}")
  message(STATUS "Libraries install   = ${CMAKE_INSTALL_LIBDIR}")
  message(STATUS "boost backend       = ${QRF_BOOST_SOURCE}")
endfunction()

# --------------------------------------------------------------------------
# Static library target
# --------------------------------------------------------------------------
# The static library is always built.  Its internal target name keeps the legacy
# namespace prefix, and an imported-style ALIAS target exposes the conventional
# utils::<project>_Static name to users and tests.
add_library(${UTILS_NAMESPACE}_${PROJECT_NAME}_Static STATIC ${SOURCES} ${PUBLIC_HEADERS})
add_library(${UTILS_NAMESPACE}::${PROJECT_NAME}_Static ALIAS ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static)
set_target_properties(${UTILS_NAMESPACE}_${PROJECT_NAME}_Static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
qrf_configure_target(${UTILS_NAMESPACE}_${PROJECT_NAME}_Static)
qrf_add_platform_library_alias(${UTILS_NAMESPACE}_${PROJECT_NAME}_Static ${PROJECT_NAME} STATIC)
list(APPEND TARGETS_TO_BE_INSTALLED ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static)

# --------------------------------------------------------------------------
# Optional shared library target
# --------------------------------------------------------------------------
# When UTILS_BUILD_SHARED is enabled, build a shared library from the same source
# and header lists.  It receives the same target configuration and naming/alias
# behavior as the static target.
if(UTILS_BUILD_SHARED)
  add_library(${UTILS_NAMESPACE}_${PROJECT_NAME} SHARED ${SOURCES} ${PUBLIC_HEADERS})
  add_library(${UTILS_NAMESPACE}::${PROJECT_NAME} ALIAS ${UTILS_NAMESPACE}_${PROJECT_NAME})
  set_target_properties(${UTILS_NAMESPACE}_${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  qrf_configure_target(${UTILS_NAMESPACE}_${PROJECT_NAME})
  qrf_add_platform_library_alias(${UTILS_NAMESPACE}_${PROJECT_NAME} ${PROJECT_NAME} SHARED)
  list(APPEND TARGETS_TO_BE_INSTALLED ${UTILS_NAMESPACE}_${PROJECT_NAME})
endif()

# --------------------------------------------------------------------------
# Test executables
# --------------------------------------------------------------------------
# Tests are optional and are built from every .cc file in src_tests/.  Each test
# links against the static library so it exercises the same archive that is
# normally installed for downstream projects.
if(UTILS_ENABLE_TESTS)
  enable_testing()

  # Base libraries needed by tests.  CMAKE_DL_LIBS is empty on platforms that do
  # not need an explicit dynamic-loader library and expands to the proper system
  # library where required.
  set(LIBS ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static ${CMAKE_DL_LIBS})

  # Linux and other non-Apple Unix systems may require an explicit pthread target.
  # macOS handles this differently, so the thread dependency is only added here.
  if(UNIX AND NOT APPLE)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    list(APPEND LIBS Threads::Threads)
  endif()

  # Discover test source files.  CONFIGURE_DEPENDS asks CMake-supported generators
  # to re-run configuration when files are added to or removed from src_tests/.
  file(
    GLOB TEST_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/*.cc"
  )

  if(TEST_SOURCES)
    add_custom_target("${PROJECT_NAME}_all_tests" ALL)
  endif()

  # Create one executable per test source.  The executable name is derived from the
  # source file basename and each executable is registered as a CTest test.
  foreach(TEST_SOURCE IN LISTS TEST_SOURCES)
    get_filename_component(TEST_NAME "${TEST_SOURCE}" NAME_WE)
    add_executable(${TEST_NAME} "${TEST_SOURCE}")
    set_target_properties(
      ${TEST_NAME}
      PROPERTIES
        SUFFIX ".exe"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    )
    if(QRF_ENABLE_ARBITRARY_PRECISION)
      target_include_directories(${TEST_NAME} SYSTEM PRIVATE "${QRF_BOOST_INCLUDE_DIR}")
    endif()
    target_link_libraries(${TEST_NAME} PRIVATE ${LIBS})
    add_test(NAME "${TEST_NAME}" COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/bin/${TEST_NAME}.exe")
    add_dependencies("${PROJECT_NAME}_all_tests" ${TEST_NAME})
  endforeach()
endif()

# --------------------------------------------------------------------------
# Install compiled targets
# --------------------------------------------------------------------------
# Install every library target collected in TARGETS_TO_BE_INSTALLED.  Static
# archives, shared libraries, runtime binaries, and public include metadata are
# routed through GNUInstallDirs variables, which may have been overridden above
# for the local ./lib layout.
install(
  TARGETS ${TARGETS_TO_BE_INSTALLED}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install headers one by one while preserving their relative directory below
# src/.  This avoids installing private build artifacts while keeping the public
# include layout identical between the source tree and the install tree.
foreach(HEADER IN LISTS PUBLIC_HEADERS)
  file(RELATIVE_PATH HEADER_REL "${SOURCE_DIR}" "${HEADER}")
  get_filename_component(HEADER_DIR "${HEADER_REL}" DIRECTORY)
  if(HEADER_DIR STREQUAL "")
    install(FILES "${HEADER}" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  else()
    install(FILES "${HEADER}" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${HEADER_DIR}")
  endif()
endforeach()

# --------------------------------------------------------------------------
# Packaging with CPack
# --------------------------------------------------------------------------
# Select a native-ish package generator per platform.  Windows uses NSIS,
# Apple uses productbuild, and Unix-like systems use DEB.  The package metadata
# is intentionally simple and follows the legacy utils-<project> naming scheme.
set(CPACK_PACKAGE_NAME "${UTILS_NAMESPACE}-${PROJECT_NAME}")
if(WIN32)
  set(CPACK_GENERATOR "NSIS64")
  set(CPACK_SET_DESTDIR OFF)
  set(CPACK_NSIS_MODIFY_PATH ON)
  set(CPACK_PACKAGE_INSTALL_DIRECTORY "CPack Component Example")
elseif(APPLE)
  set(CPACK_GENERATOR "productbuild")
  set(CPACK_SET_DESTDIR ON)
elseif(UNIX)
  set(CPACK_GENERATOR "DEB")
  set(CPACK_SET_DESTDIR ON)
endif()
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_NAME} - CPack Component Installation")
set(CPACK_PACKAGE_VENDOR "Enrico Bertolazzi")
include(CPack)

# --------------------------------------------------------------------------
# Final configuration summary
# --------------------------------------------------------------------------
# These messages are printed at configure time and provide a compact overview
# before qrf_final_messages() emits the full diagnostic dump.
message("${PROJECT_NAME} current version = v${PROJECT_VERSION}")
message(STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Platform suffix: ${UTILS_PLATFORM_SUFFIX}")
qrf_final_messages()
