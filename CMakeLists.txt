cmake_minimum_required(VERSION 3.3 FATAL_ERROR)  # 3.3: EP_add(DEPENDS interfacelib
# project_VERSION* variables undefined, handled in primary EP

project(pylibefp
        LANGUAGES C CXX)
set(pylibefp_AUTHORS      "Lori A. Burns")
set(pylibefp_DESCRIPTION  "Python binding to I. Kaliman's libefp for Effective Fragment Potentials")
set(pylibefp_URL          "https://github.com/loriab/pylibefp")
set(pylibefp_LICENSE      "BSD 3-clause")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

############################  Options: Build How?  #############################
include(psi4OptionsTools)
option_with_flags(ENABLE_XHOST "Enables processor-specific optimization" ON
                  "-xHost" "-march=native")
option_with_default(CMAKE_BUILD_TYPE "Build type (Release or Debug)" Release)
option_with_default(CMAKE_INSTALL_LIBDIR "Directory to which libraries installed" lib)
option_with_default(PYMOD_INSTALL_LIBDIR "Location within CMAKE_INSTALL_LIBDIR to which python modules are installed" /)
option_with_default(ENABLE_GENERIC "Enables mostly static linking of system libraries for shared library" OFF)
option_with_default(pylibefp_CXX_STANDARD "Specify C++ standard for core pylibefp" 11)

########################  Process & Validate Options  ##########################
include(GNUInstallDirs)
include(autocmake_safeguards)
include(custom_color_messages)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr/local/pylibefp" CACHE PATH "Install path" FORCE)
endif()
message(STATUS "pylibefp install: ${CMAKE_INSTALL_PREFIX}")

set(Python_ADDITIONAL_VERSIONS 3.7 3.6 3.5)  # adjust with CMake minimum FindPythonInterp
find_package(PythonLibsNew 2.7 REQUIRED)
message(STATUS "${Cyan}Found Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE} (found version ${PYTHON_VERSION_STRING})")

set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
add_subdirectory(external/upstream)

# external projects manage their own OpenMP and c++YY flags, so only add to CXX_FLAGS for pylibefp-core
include(custom_cxxstandard)
include(custom_static_library)

################################  Main Project  ################################
include(ExternalProject)
ExternalProject_Add(pylibefp-core
   DEPENDS #lapack_external
           libefp_external
           pybind11_external
   SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pylibefp
   CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
              -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
              -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
              -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}
              -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
              -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
              -DCMAKE_CXX_STANDARD=${pylibefp_CXX_STANDARD}
              -DCMAKE_CXX_STANDARD_REQUIRED=ON
              -DCMAKE_CXX_EXTENSIONS=OFF
              -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
              -DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}
              -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
              -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
              -DPYMOD_INSTALL_LIBDIR=${PYMOD_INSTALL_LIBDIR}
              -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}
              -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR}
              -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
              -DPYLIBEFP_ROOT=${CMAKE_CURRENT_SOURCE_DIR}
              -DTargetLAPACK_DIR=${TargetLAPACK_DIR}
              -Dlibefp_DIR=${libefp_DIR}
              -Dpybind11_DIR=${pybind11_DIR}
              -DLIBC_INTERJECT=${LIBC_INTERJECT}
              -DDESTDIR=${CMAKE_BINARY_DIR}/stage
   CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
   BUILD_ALWAYS 1
   INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
                   DESTDIR=${CMAKE_BINARY_DIR}/stage)

add_subdirectory(doc)
add_subdirectory(tests)

# <<<  Install  >>>

install(DIRECTORY ${STAGED_INSTALL_PREFIX}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        USE_SOURCE_PERMISSIONS
        PATTERN "*.pyc" EXCLUDE)

