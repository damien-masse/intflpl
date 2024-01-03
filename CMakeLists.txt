# ==================================================================
#  Codac - cmake configuration file
# ==================================================================

  cmake_minimum_required(VERSION 3.0.2)

  list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/scripts/CMakeModules/)
#  include(version_from_git)
#  version_from_git() # Obtains the version number from Git tags

  # To use a specific compiler:
  #set(CMAKE_C_COMPILER "gcc-7")
  #set(CMAKE_CXX_COMPILER "/usr/bin/g++-7")
set(VERSION 0.1)

 project(intflpl VERSION ${VERSION} LANGUAGES CXX)
#
#  if(NOT VERSION_ID)
#    set(PROJECT_VERSION_FULL ${PROJECT_VERSION})
#  else()
#    set(PROJECT_VERSION_FULL "${PROJECT_VERSION}-${VERSION_ID}")
#  endif()
#  message(STATUS "Full project version is ${PROJECT_VERSION_FULL}")

  set(PROJECT_DESCRIPTION
      "Intflpl is a small library for handling polyhedra using intervals and floating-point values.")
  set(PROJECT_LONG_DESCRIPTION
      "${PROJECT_DESCRIPTION}")
  message(STATUS "Configuring build for ${PROJECT_NAME} ${PROJECT_VERSION}")

  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
    message(STATUS "Configuring ${PROJECT_NAME} in DEBUG mode as none was specified.")
    add_compile_options(-O3)
  endif()

################################################################################
# Options for directories
################################################################################
  
  # Install directories
  set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE PATH "C++ header files (include)")
  set(CMAKE_INSTALL_LIBDIR "lib" CACHE PATH "object code libraries (lib)")
  set(CMAKE_INSTALL_BINDIR "bin" CACHE PATH "user executables (bin)")
  set(CMAKE_INSTALL_PKGCONFIG "share/pkgconfig" CACHE PATH "pkg files (share/pkgconfig)")
  set(CMAKE_INSTALL_CMAKE "share/intflpl/cmake" CACHE PATH "cmake files (share/intflpl/cmake)")


################################################################################
# Compilation configuration
################################################################################

#  # Check that the compiler supports c++11
#  include(CheckCXXCompilerFlag)
#  check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORTS_CXX17)
#
#  if(COMPILER_SUPPORTS_CXX17)
#    add_compile_options(-std=c++17)
#  else()
#    message(FATAL_ERROR "Codac needs a compiler with C++17 support")
#  endif()

#  if(WIN32)
#    # We need this for strdup under Windows (see issue #287 of ibex-lib repo)
#    add_definitions(-U__STRICT_ANSI__)
#  endif()

  #if(NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
  #endif()


################################################################################
# Looking for IBEX
################################################################################

  find_package(IBEX REQUIRED)
  ibex_init_common() # IBEX should have installed this function
  message(STATUS "Found IBEX version ${IBEX_VERSION}")



################################################################################
# Compile sources
################################################################################

  add_subdirectory(src) # C++ sources
#add_subdirectory(doc) # documentation (Doxygen + Sphinx manual)

################################################################################
# Tests
################################################################################
  
  # Note: place this tests block before the add_subdirectory(python),
  # otherwise python tests will not be taken into account. 
  option(BUILD_TESTS "Build test" TRUE)
  if(BUILD_TESTS)
    add_custom_target(check
                      COMMAND flpl-tests
                      DEPENDS intflpl COMMENT "Running the tests")

    add_subdirectory(tests)
  endif()

#  option(TEST_EXAMPLES "Testing examples" OFF)
#  add_subdirectory(examples) # examples are tested as integration tests


################################################################################
# Archives and packages
################################################################################

  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Maintainer <damien.masse@univ-brest.fr>")
  # todo: finish deb package

