cmake_minimum_required(VERSION 3.24)
project(quantreg-cpp-sample CXX)

## Configure CMAKE_MODULE_PATH to import utilities and resolve some local find_package
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake ${CMAKE_MODULE_PATH})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
enable_testing()

if (DEFINED QUANTREG_CPP_SHARED_LIBS)
  set(BUILD_SHARED_LIBS "${QUANTREG_CPP_SHARED_LIBS}")
endif ()

message("BUILD_SHARED_LIBS: ${BUILD_SHARED_LIBS}")
message("QUANTREG_CPP_SHARED_LIBS: ${QUANTREG_CPP_SHARED_LIBS}")

find_package(quantreg-cpp REQUIRED ${example_quantreg-cpp_components})


get_target_property(SHARED_COMPONENT quantreg-cpp::quantreg-cpp IMPORTED_SONAME_RELEASE)
set(AWKWARD_CONDITION_TO_GET_THE_CORRECT_BUILD OFF)
if (example_quantreg-cpp_components STREQUAL "static")
  set(AWKWARD_CONDITION_TO_GET_THE_CORRECT_BUILD OFF)
elseif ((DEFINED BUILD_SHARED_LIBS AND BUILD_SHARED_LIBS) OR SHARED_COMPONENT)
  set(AWKWARD_CONDITION_TO_GET_THE_CORRECT_BUILD ON)
endif ()

if (AWKWARD_CONDITION_TO_GET_THE_CORRECT_BUILD)
  set(BUILD_SHARED_LIBS ON)
  #### The following is always required to get the include dir of Armadillo!
  if(DEFINED ENV{Armadillo_ROOT})
    set(Armadillo_ROOT "$ENV{Armadillo_ROOT}")
  endif()
  if(DEFINED Armadillo_ROOT)
    find_package(Armadillo 12.0.1 REQUIRED NO_MODULE PATHS ${Armadillo_ROOT} NO_DEFAULT_PATH)
  else()
    find_package(Armadillo 12.0.1 REQUIRED)
  endif()
else ()
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(ZLIB_USE_STATIC_LIBS "ON")
  set(BLA_VENDOR OpenBLAS)
  set(BLA_STATIC ON)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
  find_package(ZLIB REQUIRED)
endif ()
find_package(BLAS REQUIRED)
find_package(HDF5 1.14.0 REQUIRED)

add_executable(quantreg-cpp-sample quantreg-cpp-sample.cpp)

if (AWKWARD_CONDITION_TO_GET_THE_CORRECT_BUILD)
  target_link_libraries(quantreg-cpp-sample
    PUBLIC quantreg-cpp::quantreg-cpp
    ${HDF5_LIBRARIES} ${ARMADILLO_LIBRARIES})
else ()
  ### armadillo without wrapper
  target_compile_definitions(quantreg-cpp-sample PUBLIC ARMA_DONT_USE_WRAPPER)
  target_link_libraries(quantreg-cpp-sample
    PUBLIC quantreg-cpp::quantreg-cpp
    ${HDF5_LIBRARIES} ZLIB::ZLIB ${BLAS_LIBRARIES} gfortran quadmath
  )
endif ()

file(COPY ${PROJECT_SOURCE_DIR}/../../data DESTINATION ${PROJECT_BINARY_DIR})

add_test(NAME qr_test COMMAND quantreg-cpp-sample)
# set_tests_properties(qr_test PROPERTIES
#   PASS_REGULAR_EXPRESSION "[matrix size: 365x6]"
#   ENVIRONMENT "PATH=$<TARGET_FILE_DIR:quantreg-cpp::quantreg-cpp>")

set_tests_properties(qr_test PROPERTIES
  PASS_REGULAR_EXPRESSION "   81\.4822    0\.5602"
  ENVIRONMENT "PATH=$<TARGET_FILE_DIR:quantreg-cpp::quantreg-cpp>")
