include_guard(GLOBAL)

#
# Locate the Boost.Multiprecision headers.
#
# Resolution order:
#   1. an explicit local checkout, QRF_BOOST_LOCAL_DIR (default ../boost);
#   2. the standard CMake package, find_package(Boost) - this is what works on
#      Linux distributions, in vcpkg/conan toolchains and with any Boost install
#      shipping BoostConfig.cmake;
#   3. find_path(), which honours BOOST_ROOT / BOOST_INCLUDEDIR and the usual
#      system include directories;
#   4. a handful of well known Homebrew/local prefixes, as a last resort.
#
# Only headers are needed: Boost.Multiprecision and the Boost.Math special
# functions used here are header only.
#

set(
  QRF_BOOST_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../boost"
  CACHE PATH "Path to a local Boost checkout"
)

function(qrf_resolve_boost_multiprecision out_include_dir out_source)
  set(_probe "boost/multiprecision/cpp_bin_float.hpp")

  # 1 - local sibling checkout
  get_filename_component(_local_abs "${QRF_BOOST_LOCAL_DIR}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  if(EXISTS "${_local_abs}/${_probe}")
    set(${out_include_dir} "${_local_abs}" PARENT_SCOPE)
    set(${out_source} "local: ${_local_abs}" PARENT_SCOPE)
    return()
  endif()

  # 2 - the standard CMake package
  find_package(Boost QUIET)
  if(Boost_FOUND AND Boost_INCLUDE_DIRS)
    foreach(_dir IN LISTS Boost_INCLUDE_DIRS)
      if(EXISTS "${_dir}/${_probe}")
        set(${out_include_dir} "${_dir}" PARENT_SCOPE)
        set(${out_source} "find_package(Boost): ${_dir}" PARENT_SCOPE)
        return()
      endif()
    endforeach()
  endif()

  # 3 - BOOST_ROOT / BOOST_INCLUDEDIR and the system include directories
  find_path(
    QRF_BOOST_INCLUDE_PATH
    "${_probe}"
    HINTS
      ENV BOOST_ROOT
      ENV BOOST_INCLUDEDIR
    PATH_SUFFIXES
      include
  )
  if(QRF_BOOST_INCLUDE_PATH)
    set(${out_include_dir} "${QRF_BOOST_INCLUDE_PATH}" PARENT_SCOPE)
    set(${out_source} "find_path: ${QRF_BOOST_INCLUDE_PATH}" PARENT_SCOPE)
    return()
  endif()

  # 4 - well known prefixes, kept for convenience on macOS
  foreach(_prefix
    "/opt/homebrew/include"
    "/opt/homebrew/opt/boost/include"
    "/usr/local/include"
    "/usr/local/opt/boost/include"
  )
    if(EXISTS "${_prefix}/${_probe}")
      set(${out_include_dir} "${_prefix}" PARENT_SCOPE)
      set(${out_source} "system: ${_prefix}" PARENT_SCOPE)
      return()
    endif()
  endforeach()

  set(${out_include_dir} "" PARENT_SCOPE)
  set(${out_source} "not found" PARENT_SCOPE)
endfunction()
