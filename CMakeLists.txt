############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      <mailto:enrico.bertolazzi@unitn.it>                                 #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_INSTALL_MESSAGE NEVER)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

if(ENABLE_CTEST)
  enable_testing()
endif()

# get version from the git tag
execute_process(COMMAND git describe --tags --abbrev=0 OUTPUT_VARIABLE VERSION)
string(REGEX REPLACE "\n$" "" VERSION "${VERSION}")

project(
  Quartic
  VERSION ${VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/quarticRootsFlocke/"
)

# NOTE: Predictable target name
set(STATIC_PROJECT "${PROJECT_NAME}Static")
set(DYNAMIC_PROJECT "${PROJECT_NAME}")

# FIXME: This should refer to functions that acts on targets
include(./CMakeLists-cflags.txt)
include(./CMakeUtilities.txt)

message(STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}")
set(CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
set(CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG})

# Evaluating the suffix to append to target built by this CMake file
ebertolazzi_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if (BUILD_SHARED)
  ebertolazzi_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

# NOTE: Glob shall be avoided. Headers and sources shall be declared 
# manually in order to rebuild cmake-cache automatically.
# Ref: https://cmake.org/cmake/help/v3.15/command/file.html?highlight=glob#filesystem
set(SOURCES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(SOURCES
  "${SOURCES_DIR}/PolynomialRoots-1-Quadratic.cc"
  "${SOURCES_DIR}/PolynomialRoots-2-Cubic.cc"
  "${SOURCES_DIR}/PolynomialRoots-3-Quartic.cc"
  "${SOURCES_DIR}/PolynomialRoots-Jenkins-Traub.cc"
  "${SOURCES_DIR}/PolynomialRoots-Utils.cc")
file(GLOB HEADERS "${SOURCES_DIR}/*.hh")

add_library(${STATIC_PROJECT} STATIC ${SOURCES})
target_include_directories(
  ${STATIC_PROJECT} PUBLIC 
  "$<BUILD_INTERFACE:${SOURCES_DIR}>"
  "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>"
)
set_target_properties(
  ${STATIC_PROJECT} PROPERTIES 
  OUTPUT_NAME "${PROJECT_NAME}${ARTIFACTS_STATIC_SUFFIX}"
  PUBLIC_HEADER "${HEADERS}"
  POSITION_INDEPENDENT_CODE ON
)

if(BUILD_SHARED)
  add_library(${DYNAMIC_PROJECT} SHARED ${SOURCES})
  target_include_directories(${DYNAMIC_PROJECT} PUBLIC "${SOURCES_DIR}")
  set_target_properties(
    ${DYNAMIC_PROJECT} PROPERTIES 
    OUTPUT_NAME "${PROJECT_NAME}${ARTIFACTS_STATIC_SUFFIX}"
    PUBLIC_HEADER "${HEADERS}"
    POSITION_INDEPENDENT_CODE ON
  )
  if(WIN32)
    set_target_properties(
      ${DYNAMIC_PROJECT} PROPERTIES 
      WINDOWS_EXPORT_ALL_SYMBOLS TRUE
    )
  endif()
endif()

if(BUILD_EXECUTABLE OR ENABLE_CTEST)
  set(TEST_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin")
  file(MAKE_DIRECTORY "${TEST_OUTPUT_PATH}")
  set(EXECUTABLES check_1_quadratic  check_2_cubic check_3_quartic)
  foreach(EXE ${EXECUTABLES})
    add_executable(${EXE} "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${EXE}.cc")
    target_link_libraries(${EXE} PRIVATE ${STATIC_PROJECT})
    set_target_properties(
      ${EXE} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${TEST_OUTPUT_PATH}"
    )
    add_test(
      NAME ${EXE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/bin/${EXE}"
    )
  endforeach()
endif()

# INSTALL TARGETS

# FIXME: CMakeLists should not override CMAKE_INSTALL_PREFIX
# ref: https://gitlab.kitware.com/cmake/cmake/-/issues/19322
if(INSTALL_HERE)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/lib")
  set(PROJECT_PUBLIC_HEADER_DESTINATION "${CMAKE_INSTALL_PREFIX}/include")
else()
  set(PROJECT_PUBLIC_HEADER_DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}")
endif()

set(INSTALL_TARGETS ${STATIC_PROJECT})
if(BUILD_SHARED)
  set(INSTALL_TARGETS ${INSTALL_TARGETS} ${DYNAMIC_PROJECT})
endif()
install(
  TARGETS ${INSTALL_TARGETS}
  RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
  LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/dll"
  PUBLIC_HEADER DESTINATION "${PROJECT_PUBLIC_HEADER_DESTINATION}"
)

message(STATUS "Using ${SSE_FLAGS} extensions")
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
message(STATUS "CMAKE_INSTALL_PREFIX        = ${CMAKE_INSTALL_PREFIX}")
message(STATUS "BUILD_EXECUTABLE            = ${BUILD_EXECUTABLE}")
