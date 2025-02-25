cmake_minimum_required(VERSION 3.19)
project(QUANTREG_CPP LANGUAGES CXX)

cmake_policy(SET CMP0074 NEW)

## Build compile_commands.json for later use
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

include(utils)
check_for_in_source_build()
release_as_default_build_type()

# option(QUANTREG_CPP_SHARED_LIBS "Build quantreg-cpp shared library" ${BUILD_SHARED_LIBS})
if (DEFINED QUANTREG_CPP_SHARED_LIBS)
  set(BUILD_SHARED_LIBS "${QUANTREG_CPP_SHARED_LIBS}")
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-Wall -Wextra)
endif()

if (NOT BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  ### armadillo without wrapper => target_compile_definitions(TARGET PRIVATE ARMA_DONT_USE_WRAPPER)
  set(BLA_VENDOR OpenBLAS)
  set(BLA_STATIC ON)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  find_package(BLAS REQUIRED)
endif ()

### libhdf5
find_package(HDF5 1.14.0 REQUIRED CXX)

#### The following is always required to get the include dir of Armadillo!
if(DEFINED ENV{Armadillo_ROOT})
  set(Armadillo_ROOT "$ENV{Armadillo_ROOT}")
endif()
if(DEFINED Armadillo_ROOT)
  find_package(Armadillo 10.5.3 REQUIRED NO_MODULE PATHS ${Armadillo_ROOT} NO_DEFAULT_PATH)
else()
  find_package(Armadillo 10.5.3 REQUIRED)
endif()

## VERSION
find_package(Git)
if(Git_FOUND)
  include(GetGitRevisionDescription)
  get_git_head_revision(REFSPEC COMMITHASH)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=8 --dirty
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE res
    OUTPUT_VARIABLE VERSION
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(NOT res EQUAL 0)
    message("${BoldCyan}Git found but failed to retrieve usefull release information…${ColourReset}")
    set(CMAKE_MODULE_PATH
      ${PROJECT_SOURCE_DIR}
      ${CMAKE_MODULE_PATH})
    include(.release OPTIONAL RESULT_VARIABLE res)
    if(res STREQUAL "NOTFOUND")
      message(FATAL_ERROR "${BoldRed}Unable to get version information from git or .release.cmake!!!${ColourReset}")
    endif()
  else()
    string(REGEX MATCH "v([0-9]*)\.([0-9]*)\.([0-9]*)(.*)" _ ${VERSION})

    set(QUANTREG_CPP_VERSION "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}${CMAKE_MATCH_4}")
    set(QUANTREG_CPP_VERSION_MAJOR ${CMAKE_MATCH_1})
    set(QUANTREG_CPP_VERSION_MINOR ${CMAKE_MATCH_2})
    set(QUANTREG_CPP_VERSION_PATCH ${CMAKE_MATCH_3})
  endif()
else()
  message("${BoldCyan}Git not found, trying to retrieve release information from .release.cmake…${ColourReset}")
  set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}
    ${CMAKE_MODULE_PATH})
  include(.release OPTIONAL RESULT_VARIABLE res)
  if(res STREQUAL "NOTFOUND")
    message(FATAL_ERROR "${BoldRed}Unable to get version information from git or .release.cmake!!!${ColourReset}")
  endif()

endif()
message("quantreg-cpp version: ${BoldMagenta}${QUANTREG_CPP_VERSION}${ColourReset}")

add_library(quantreg-cpp src/quantreg-cpp.cpp)
add_library(quantreg-cpp::quantreg-cpp ALIAS quantreg-cpp)

set_target_properties(quantreg-cpp PROPERTIES VERSION ${QUANTREG_CPP_VERSION} SOVERSION ${QUANTREG_CPP_VERSION_MAJOR})
if (NOT BUILD_SHARED_LIBS)
  set_property(TARGET quantreg-cpp PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()

target_compile_definitions(quantreg-cpp
  PUBLIC "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:QUANTREG_CPP_STATIC_DEFINE>"
  PUBLIC ARMA_USE_HDF5)

target_include_directories(quantreg-cpp
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  # PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>"
  PUBLIC ${ARMADILLO_INCLUDE_DIRS})

#target_compile_features(quantreg-cpp PUBLIC cxx_std_17)

if (BUILD_SHARED_LIBS)
  target_link_libraries(quantreg-cpp PUBLIC ${ARMADILLO_LIBRARIES})
else ()
  target_compile_definitions(quantreg-cpp PUBLIC ARMA_DONT_USE_WRAPPER)
  target_link_libraries(quantreg-cpp PUBLIC BLAS::BLAS)
endif ()

## Generate the export header for quantreg and attach it to the target
include(GenerateExportHeader)
generate_export_header(quantreg-cpp EXPORT_FILE_NAME include/quantreg-cpp_export.h)

## Include the install rules if the user wanted them (included by default when top-level)
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}" is_top_level)
option(QUANTREG_CPP_INCLUDE_PACKAGING "Include packaging rules for quantreg-cpp" "${is_top_level}")
if (QUANTREG_CPP_INCLUDE_PACKAGING)
  add_subdirectory(packaging)
endif ()
