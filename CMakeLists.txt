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

cmake_minimum_required(VERSION 3.16)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.31)
  cmake_policy(SET CMP0177 NEW)
endif()
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()
if(POLICY CMP0169)
  cmake_policy(SET CMP0169 OLD)
endif()

option(UTILS_GIT_VERSION "Get version from GIT command" ON)
option(UTILS_ENABLE_TESTS "Enable compilation of testing targets" ON)
option(UTILS_BUILD_SHARED "Enable shared option compilation" OFF)
option(UTILS_INSTALL_LOCAL "Enable installation on local directory" ON)
option(UTILS_TWEAK_VERSION_IS_TIMESTAMP "Sets tweak version to timestamp" OFF)

set(UTILS_NAMESPACE "utils" CACHE STRING "Sets the namespace used for exports")
set(TARGETS_TO_BE_INSTALLED)

set(CMAKE_INSTALL_MESSAGE NEVER)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_VERBOSE_MAKEFILE FALSE)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_INSTALL_RPATH $ORIGIN)
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)

if(NOT UNIX)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

add_compile_definitions(_USE_MATH_DEFINES)

set(
  UTILS_CURRENT_VERSION "1.0"
  CACHE STRING "Default project version, if not taken from git"
)

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

project(
  quarticRootsFlocke
  VERSION ${UTILS_PROJECT_VERSION}
  HOMEPAGE_URL "https://ebertolazzi.github.io/quarticRootsFlocke/"
  LANGUAGES CXX
)

include(GNUInstallDirs)
include(CTest)

if(PROJECT_IS_TOP_LEVEL AND (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR CMAKE_INSTALL_PREFIX STREQUAL "/usr/local"))
  set(
    CMAKE_INSTALL_PREFIX
    "${CMAKE_CURRENT_SOURCE_DIR}"
    CACHE PATH "Install path prefix, prepended onto install directories."
    FORCE
  )
endif()

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

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

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

set(VS_DEFINE "/D_WINDOWS /D_SCL_SECURE_NO_WARNINGS /DHAVE_STRING_H /DNO_GETTIMEOFDAY /DYAML_DECLARE_STATIC /DPCRE_STATIC")
set(VS_FLAGS "/nologo /std:c++17 /utf-8 /GS /Gy /GF /EHsc /bigobj")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "-std=c++17 ${COMMON_CXX_FLAGS_STR}")
  set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -O2 -g ")
  set(CMAKE_CXX_FLAGS_DEBUG "-fPIC -O0 -ggdb ")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++17 -stdlib=libc++ ${COMMON_CXX_FLAGS_STR}")
  set(CMAKE_CXX_FLAGS_RELEASE "-fPIC -Os -g -fcolor-diagnostics ")
  set(CMAKE_CXX_FLAGS_DEBUG "-fPIC -O0 -gfull -fcolor-diagnostics ")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  set(CMAKE_CXX_FLAGS " ${MSVC_COMPILE_FLAGS_STR} ${VS_DEFINE} /sdl")
  set(CMAKE_CXX_FLAGS_RELEASE " ${VS_FLAGS} /O2 /MD /Zo ")
  set(CMAKE_CXX_FLAGS_DEBUG " ${VS_FLAGS} /Od /Ob0 /MDd /guard:cf /Zi /DEBUG:FULL ")
else()
  message(FATAL_ERROR "Unsupported compiler ${CMAKE_CXX_COMPILER_ID}")
endif()

set(CMAKE_C_FLAGS "")
set(CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
set(CMAKE_C_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})

if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/../CMakeLists-customize.txt")
  include("${CMAKE_CURRENT_LIST_DIR}/../CMakeLists-customize.txt")
endif()

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

set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

file(
  GLOB SOURCES
  CONFIGURE_DEPENDS
  "${SOURCE_DIR}/*.cc"
)

file(
  GLOB_RECURSE PUBLIC_HEADERS
  CONFIGURE_DEPENDS
  "${SOURCE_DIR}/*.h"
  "${SOURCE_DIR}/*.hh"
  "${SOURCE_DIR}/*.hpp"
)

foreach(H IN LISTS PUBLIC_HEADERS)
  set_source_files_properties("${H}" PROPERTIES HEADER_FILE_ONLY TRUE)
endforeach()

set(INCLUDE_DIRS "${SOURCE_DIR}")

function(qrf_configure_target target_name)
  target_compile_features(${target_name} PUBLIC cxx_std_17)
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

function(qrf_add_platform_library_alias target_name base_name library_kind)
  if(library_kind STREQUAL "STATIC")
    set(target_filename "${CMAKE_STATIC_LIBRARY_PREFIX}${base_name}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if(WIN32)
      set(alias_filename "lib${base_name}_${UTILS_PLATFORM_SUFFIX}_static${CMAKE_STATIC_LIBRARY_SUFFIX}")
    else()
      set(alias_filename "${CMAKE_STATIC_LIBRARY_PREFIX}${base_name}_${UTILS_PLATFORM_SUFFIX}_static${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()
  else()
    set(target_filename "${CMAKE_SHARED_LIBRARY_PREFIX}${base_name}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    if(WIN32)
      set(alias_filename "lib${base_name}_${UTILS_PLATFORM_SUFFIX}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    else()
      set(alias_filename "${CMAKE_SHARED_LIBRARY_PREFIX}${base_name}_${UTILS_PLATFORM_SUFFIX}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    endif()
  endif()

  if(WIN32)
    add_custom_command(
      TARGET ${target_name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              "$<TARGET_FILE:${target_name}>"
              "$<TARGET_FILE_DIR:${target_name}>/${alias_filename}"
      VERBATIM
    )

    if(UTILS_INSTALL_LOCAL)
      install(CODE
        "file(MAKE_DIRECTORY \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")\n"
        "execute_process(\n"
        "  COMMAND \"${CMAKE_COMMAND}\" -E copy_if_different\n"
        "          \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${target_filename}\"\n"
        "          \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${alias_filename}\"\n"
        ")"
      )
    endif()
  else()
    add_custom_command(
      TARGET ${target_name} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E rm -f
              "$<TARGET_FILE_DIR:${target_name}>/${alias_filename}"
      COMMAND ${CMAKE_COMMAND} -E create_symlink
              "${target_filename}"
              "$<TARGET_FILE_DIR:${target_name}>/${alias_filename}"
      VERBATIM
    )

    if(UTILS_INSTALL_LOCAL)
      install(CODE
        "file(MAKE_DIRECTORY \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}\")\n"
        "execute_process(\n"
        "  COMMAND \"${CMAKE_COMMAND}\" -E rm -f\n"
        "          \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${alias_filename}\"\n"
        ")\n"
        "execute_process(\n"
        "  COMMAND \"${CMAKE_COMMAND}\" -E create_symlink\n"
        "          \"${target_filename}\"\n"
        "          \"\$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${alias_filename}\"\n"
        ")"
      )
    endif()
  endif()
endfunction()

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
  message(STATUS "UTILS_INSTALL_LOCAL = ${UTILS_INSTALL_LOCAL}")
  message(STATUS "Headers install     = ${CMAKE_INSTALL_INCLUDEDIR}")
  message(STATUS "Libraries install   = ${CMAKE_INSTALL_LIBDIR}")
endfunction()

add_library(${UTILS_NAMESPACE}_${PROJECT_NAME}_Static STATIC ${SOURCES} ${PUBLIC_HEADERS})
add_library(${UTILS_NAMESPACE}::${PROJECT_NAME}_Static ALIAS ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static)
set_target_properties(${UTILS_NAMESPACE}_${PROJECT_NAME}_Static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
qrf_configure_target(${UTILS_NAMESPACE}_${PROJECT_NAME}_Static)
qrf_add_platform_library_alias(${UTILS_NAMESPACE}_${PROJECT_NAME}_Static ${PROJECT_NAME} STATIC)
list(APPEND TARGETS_TO_BE_INSTALLED ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static)

if(UTILS_BUILD_SHARED)
  add_library(${UTILS_NAMESPACE}_${PROJECT_NAME} SHARED ${SOURCES} ${PUBLIC_HEADERS})
  add_library(${UTILS_NAMESPACE}::${PROJECT_NAME} ALIAS ${UTILS_NAMESPACE}_${PROJECT_NAME})
  set_target_properties(${UTILS_NAMESPACE}_${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  qrf_configure_target(${UTILS_NAMESPACE}_${PROJECT_NAME})
  qrf_add_platform_library_alias(${UTILS_NAMESPACE}_${PROJECT_NAME} ${PROJECT_NAME} SHARED)
  list(APPEND TARGETS_TO_BE_INSTALLED ${UTILS_NAMESPACE}_${PROJECT_NAME})
endif()

if(UTILS_ENABLE_TESTS)
  enable_testing()

  set(LIBS ${UTILS_NAMESPACE}_${PROJECT_NAME}_Static ${CMAKE_DL_LIBS})

  if(UNIX AND NOT APPLE)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    list(APPEND LIBS Threads::Threads)
  endif()

  file(
    GLOB TEST_SOURCES
    CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/*.cc"
  )

  if(TEST_SOURCES)
    add_custom_target("${PROJECT_NAME}_all_tests" ALL)
  endif()

  foreach(TEST_SOURCE IN LISTS TEST_SOURCES)
    get_filename_component(TEST_NAME "${TEST_SOURCE}" NAME_WE)
    add_executable(${TEST_NAME} "${TEST_SOURCE}")
    set_target_properties(
      ${TEST_NAME}
      PROPERTIES
        SUFFIX ".exe"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
    )
    target_link_libraries(${TEST_NAME} PRIVATE ${LIBS})
    add_test(NAME "${TEST_NAME}" COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/bin/${TEST_NAME}.exe")
    add_dependencies("${PROJECT_NAME}_all_tests" ${TEST_NAME})
  endforeach()
endif()

install(
  TARGETS ${TARGETS_TO_BE_INSTALLED}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

foreach(HEADER IN LISTS PUBLIC_HEADERS)
  file(RELATIVE_PATH HEADER_REL "${SOURCE_DIR}" "${HEADER}")
  get_filename_component(HEADER_DIR "${HEADER_REL}" DIRECTORY)
  if(HEADER_DIR STREQUAL "")
    install(FILES "${HEADER}" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
  else()
    install(FILES "${HEADER}" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${HEADER_DIR}")
  endif()
endforeach()

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

message("${PROJECT_NAME} current version = v${PROJECT_VERSION}")
message(STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "Platform suffix: ${UTILS_PLATFORM_SUFFIX}")
qrf_final_messages()
