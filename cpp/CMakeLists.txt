cmake_minimum_required(VERSION 3.18)

project(ow3d_directional_generator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(MF12_ENABLE_FFTW "Enable FFTW backend when a compatible FFTW3 install is available" OFF)
option(MF12_ENABLE_OPENMP "Enable OpenMP for hotspot parallelization when available" OFF)
set(MF12_FFTW_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/fftw" CACHE PATH "Optional local FFTW installation root")

if(MF12_ENABLE_FFTW)
  find_package(FFTW3 CONFIG QUIET)
  if(NOT FFTW3_FOUND AND EXISTS "${MF12_FFTW_ROOT}")
    find_package(FFTW3 CONFIG QUIET PATHS "${MF12_FFTW_ROOT}" NO_DEFAULT_PATH)
  endif()
  if(NOT FFTW3_FOUND)
    find_path(
      MF12_FFTW_INCLUDE_DIR
      NAMES fftw3.h
      PATHS
        "${MF12_FFTW_ROOT}/include"
    )
    find_library(
      MF12_FFTW_LIBRARY
      NAMES fftw3 libfftw3 libfftw3-3
      PATHS
        "${MF12_FFTW_ROOT}/lib"
    )
    find_file(
      MF12_FFTW_DLL
      NAMES libfftw3-3.dll fftw3-3.dll
      PATHS
        "${MF12_FFTW_ROOT}/bin"
    )
    if(MF12_FFTW_INCLUDE_DIR AND MF12_FFTW_LIBRARY)
      add_library(MF12::fftw3 UNKNOWN IMPORTED)
      set_target_properties(
        MF12::fftw3
        PROPERTIES
          IMPORTED_LOCATION "${MF12_FFTW_LIBRARY}"
          INTERFACE_INCLUDE_DIRECTORIES "${MF12_FFTW_INCLUDE_DIR}"
      )
      if(MF12_FFTW_DLL)
        set_target_properties(MF12::fftw3 PROPERTIES IMPORTED_RUNTIME_ARTIFACTS "${MF12_FFTW_DLL}")
      endif()
      set(FFTW3_FOUND TRUE)
    endif()
  endif()
endif()

if(MF12_ENABLE_OPENMP)
  find_package(OpenMP QUIET)
endif()

add_library(
  ow3d_directional_lib
  src/config.cpp
  src/four_phase.cpp
  src/mf12_case_io.cpp
  src/mf12_compare.cpp
  src/mf12_spectral.cpp
  src/mf12_runner.cpp
  src/ow3d_export.cpp
  src/spectrum.cpp
  src/verification.cpp
  src/visualization.cpp
)

target_include_directories(
  ow3d_directional_lib
  PUBLIC
    include
    third_party
)

if(MF12_ENABLE_FFTW AND FFTW3_FOUND)
  if(TARGET FFTW3::fftw3)
    target_link_libraries(ow3d_directional_lib PUBLIC FFTW3::fftw3)
  else()
    target_link_libraries(ow3d_directional_lib PUBLIC MF12::fftw3)
    target_include_directories(ow3d_directional_lib PRIVATE "${MF12_FFTW_INCLUDE_DIR}")
  endif()
  target_compile_definitions(ow3d_directional_lib PRIVATE MF12_HAVE_FFTW=1)
endif()

if(MF12_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(ow3d_directional_lib PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(ow3d_directional_lib PRIVATE MF12_HAVE_OPENMP=1)
endif()

add_executable(ow3d_directional_generator src/main.cpp)
target_link_libraries(ow3d_directional_generator PRIVATE ow3d_directional_lib)

if(MF12_ENABLE_FFTW AND MF12_FFTW_DLL)
  add_custom_command(
    TARGET ow3d_directional_generator POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${MF12_FFTW_DLL}"
    "$<TARGET_FILE_DIR:ow3d_directional_generator>"
  )
endif()

message(STATUS "MF12_ENABLE_FFTW=${MF12_ENABLE_FFTW}")
message(STATUS "MF12_ENABLE_OPENMP=${MF12_ENABLE_OPENMP}")
if(MF12_ENABLE_FFTW)
  if(FFTW3_FOUND)
    message(STATUS "FFTW detected: enabled")
  else()
    message(WARNING "FFTW requested but not found. Build will continue without FFTW.")
  endif()
endif()
if(MF12_ENABLE_OPENMP)
  if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP detected: enabled")
  else()
    message(WARNING "OpenMP requested but not found. Build will continue without OpenMP.")
  endif()
endif()

enable_testing()
add_executable(ow3d_directional_tests tests/test_main.cpp)
target_link_libraries(ow3d_directional_tests PRIVATE ow3d_directional_lib)
add_test(NAME ow3d_directional_tests COMMAND $<TARGET_FILE:ow3d_directional_tests>)
