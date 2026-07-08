include_guard(GLOBAL)

set(
  QRF_BOOST_LOCAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../boost"
  CACHE PATH "Path to a local Boost checkout"
)

function(qrf_resolve_boost_multiprecision out_include_dir out_source)
  get_filename_component(_local_abs "${QRF_BOOST_LOCAL_DIR}" ABSOLUTE BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  set(_local_probe "${_local_abs}/boost/multiprecision/cpp_bin_float.hpp")
  if(EXISTS "${_local_probe}")
    set(${out_include_dir} "${_local_abs}" PARENT_SCOPE)
    set(${out_source} "local: ${_local_abs}" PARENT_SCOPE)
    return()
  endif()

  set(_brew_probe "/opt/homebrew/include/boost/multiprecision/cpp_bin_float.hpp")
  if(EXISTS "${_brew_probe}")
    set(${out_include_dir} "/opt/homebrew/include" PARENT_SCOPE)
    set(${out_source} "system: /opt/homebrew/include" PARENT_SCOPE)
    return()
  endif()

  set(_brew_opt_probe "/opt/homebrew/opt/boost/include/boost/multiprecision/cpp_bin_float.hpp")
  if(EXISTS "${_brew_opt_probe}")
    set(${out_include_dir} "/opt/homebrew/opt/boost/include" PARENT_SCOPE)
    set(${out_source} "system: /opt/homebrew/opt/boost/include" PARENT_SCOPE)
    return()
  endif()

  set(_usr_local_probe "/usr/local/include/boost/multiprecision/cpp_bin_float.hpp")
  if(EXISTS "${_usr_local_probe}")
    set(${out_include_dir} "/usr/local/include" PARENT_SCOPE)
    set(${out_source} "system: /usr/local/include" PARENT_SCOPE)
    return()
  endif()

  set(_usr_local_opt_probe "/usr/local/opt/boost/include/boost/multiprecision/cpp_bin_float.hpp")
  if(EXISTS "${_usr_local_opt_probe}")
    set(${out_include_dir} "/usr/local/opt/boost/include" PARENT_SCOPE)
    set(${out_source} "system: /usr/local/opt/boost/include" PARENT_SCOPE)
    return()
  endif()

  find_path(
    _boost_include_dir
    boost/multiprecision/cpp_bin_float.hpp
    HINTS
      ENV BOOST_ROOT
      ENV BOOST_INCLUDEDIR
    PATH_SUFFIXES
      include
  )
  if(_boost_include_dir)
    set(${out_include_dir} "${_boost_include_dir}" PARENT_SCOPE)
    set(${out_source} "find_path: ${_boost_include_dir}" PARENT_SCOPE)
    return()
  endif()

  set(${out_include_dir} "" PARENT_SCOPE)
  set(${out_source} "not found" PARENT_SCOPE)
endfunction()
