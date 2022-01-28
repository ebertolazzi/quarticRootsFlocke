############################################################################
#                                                                          #
#  file: CMakeLists.txt                                                    #
#  Copyright (C) 2016                                                      #
#      Enrico Bertolazzi                                                   #
#      Dipartimento di Ingegneria Industriale                              #
#      Universita` degli Studi di Trento                                   #
#      email: enrico.bertolazzi@ing.unitn.it                               #
#                                                                          #
############################################################################

cmake_minimum_required(VERSION 3.14)

project( quarticRootsFlocke )
include(./cmake/CMakeLists-common.txt)

project(
  ${PROJECT_NAME}
  VERSION ${EB_PROJECT_VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/quarticRootsFlocke/"
)

include(./cmake/CMakeLists-cflags.txt)
include(./cmake/CMakeLists-utilities.txt)

# Evaluating the suffix to append to target built by this CMake file
ebertolazzi_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if( EB_BUILD_SHARED )
  ebertolazzi_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )
message( STATUS "BASE SUFFIXES = ${ARTIFACTS_STATIC_SUFFIX} ${ARTIFACTS_DYNAMIC_SUFFIX}" )

if(APPLE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLUA_USE_MACOSX")
elseif(LINUX OR UNIX)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLUA_USE_LINUX")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
endif()

set(LUALIB Lua)

#   ___
#  / __| ___ _  _ _ _ __ ___ ___
#  \__ \/ _ \ || | '_/ _/ -_|_-<
#  |___/\___/\_,_|_| \__\___/__/
#
set( SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" )

set( SOURCES )
file(
  GLOB S
  ${SOURCE_DIR}/*.cc
)
foreach(F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach(F ${S})

#   _  _             _
#  | || |___ __ _ __| |___ _ _ ___
#  | __ / -_) _` / _` / -_) '_(_-<
#  |_||_\___\__,_\__,_\___|_| /__/
#
set( PUBLIC_HEADERS )
file(
  GLOB_RECURSE S
  ${SOURCE_DIR}/*.h*
  ${SOURCE_DIR}/Eigen/*
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF "${SOURCE_DIR}" "${F}" )
  list( APPEND PUBLIC_HEADERS ${RF} )
endforeach(F ${S})

list(FILTER PUBLIC_HEADERS EXCLUDE REGEX "\.DS")

# extra include directory for compilation
include_directories( ./src ./src/Utils )

#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#
# define target
ebertolazzi_setup_target(
  ${PROJECT_NAME}
  TRUE # build static
  "${SOURCE_DIR}"
  "${SOURCES}"
  "${SOURCE_DIR}"
  "${PUBLIC_HEADERS}"
)
if ( EB_BUILD_SHARED )
  ebertolazzi_setup_target(
    ${PROJECT_NAME}
    FALSE # build dynamic
    "${SOURCE_DIR}"
    "${SOURCES}"
    "${SOURCE_DIR}"
    "${PUBLIC_HEADERS}"
  )
endif()

#   _____       _
#  |_   _|__ __| |_ ___
#    | |/ -_|_-<  _(_-<
#    |_|\___/__/\__/__/
#
if ( EB_ENABLE_TESTS )

  enable_testing()

  if(APPLE OR UNIX)
    find_package( ZLIB REQUIRED )
  else()
    if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
      find_library( ZLIB_LIBRARIES zlibstaticd HINTS ./ThirdParties/zlib/lib )
    else()
      find_library( ZLIB_LIBRARIES zlibstatic HINTS ./ThirdParties/zlib/lib )
    endif()
    set( ZLIB_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParties/zlib/include ")
  endif()

  set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin )
  set( LIBS ${EB_NAMESPACE}_${PROJECT_NAME}_Static ${ZLIB_LIBRARIES} ${CMAKE_DL_LIBS} )

  if( UNIX AND NOT APPLE )
    set( THREADS_PREFER_PTHREAD_FLAG ON )
    find_package( Threads REQUIRED )
    set( LIBS ${LIBS} Threads::Threads )
  endif()

  set( EXELISTCPP )
  file(
    GLOB S
    ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/*.cc
  )
  foreach(F ${S})
    file( RELATIVE_PATH RF "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/" "${F}" )
    get_filename_component( RFF ${RF} NAME_WLE )
    list( APPEND EXELISTCPP ${RFF} )
  endforeach(F ${S})

  add_custom_target( "${PROJECT_NAME}_all_tests" ALL )

  foreach( S ${EXELISTCPP} )
    add_executable( ${S} ${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc )
    target_link_libraries( ${S} ${LIBS} )
    add_test( NAME "${S}" COMMAND ./bin/${S} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( "${PROJECT_NAME}_all_tests" ${S} )
  endforeach()

endif()

#   ___         _        _ _
#  |_ _|_ _  __| |_ __ _| | |
#   | || ' \(_-<  _/ _` | | |
#  |___|_||_/__/\__\__,_|_|_|
#
install(
  TARGETS ${TARGETS_TO_BE_INSTALLED}
  # EXPORT (da definire)
  LIBRARY DESTINATION ${INSTALL_PREFIX}/dll
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${INSTALL_PREFIX}/bin
  PUBLIC_HEADER DESTINATION ${PROJECT_PUBLIC_HEADER_DESTINATION}
)

#   ___         _
#  | _ \__ _ __| |____ _ __ _ ___
#  |  _/ _` / _| / / _` / _` / -_)
#  |_| \__,_\__|_\_\__,_\__, \___|
#                       |___/
#

include( ./cmake/CMakeLists-cpack.txt )

ebertolazzi_final_messages()
