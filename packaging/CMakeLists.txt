include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if (NOT DEFINED QUANTREG_CPP_INSTALL_CMAKEDIR)
  set(QUANTREG_CPP_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/quantreg-cpp"
    CACHE STRING "Path to quantreg-cpp CMake files")
endif ()

install(TARGETS quantreg-cpp EXPORT quantreg-cpp_Targets
  RUNTIME COMPONENT quantreg-cpp_Runtime
  LIBRARY COMPONENT quantreg-cpp_Runtime
  NAMELINK_COMPONENT quantreg-cpp_Development
  ARCHIVE COMPONENT quantreg-cpp_Development
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(DIRECTORY "${QUANTREG_CPP_SOURCE_DIR}/include/" "${QUANTREG_CPP_BINARY_DIR}/include/"
  TYPE INCLUDE
  COMPONENT quantreg-cpp_Development)

if (BUILD_SHARED_LIBS)
  set(type shared)
else ()
  set(type static)
endif ()

install(EXPORT quantreg-cpp_Targets
  DESTINATION "${QUANTREG_CPP_INSTALL_CMAKEDIR}"
  NAMESPACE quantreg-cpp::
  FILE quantreg-cpp-${type}-targets.cmake
  COMPONENT quantreg-cpp_Development)

write_basic_package_version_file(quantreg-cpp-config-version.cmake
  VERSION ${QUANTREG_CPP_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/quantreg-cpp-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/quantreg-cpp-config-version.cmake"
  DESTINATION "${QUANTREG_CPP_INSTALL_CMAKEDIR}"
  COMPONENT quantreg-cpp_Development)

# TODO: add additional CPack variables here
include(CPack)
