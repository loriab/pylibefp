cmake_minimum_required(VERSION 3.0)
# project_VERSION* variables populated from project(... VERSION x.x.x) string
cmake_policy(SET CMP0048 NEW)

project(pylibefp
        VERSION 0.1.0
        LANGUAGES CXX C)
set(pylibefp_AUTHORS      "Lori A. Burns")
set(pylibefp_DESCRIPTION  "Python binding to I. Kaliman's libefp for Effective Fragment Potentials")
set(pylibefp_URL          "https://github.com/loriab/pylibefp")
set(pylibefp_LICENSE      "BSD 3-clause")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

################################### Options ####################################
include(psi4OptionsTools)
option_with_default(CMAKE_BUILD_TYPE "Build type" Release)
option_with_print(ENABLE_GENERIC "Enable mostly static linking in shared library" OFF)
option_with_flags(ENABLE_XHOST "Enable processor-specific optimization" ON
                  "-xHost" "-march=native")
#option_with_print(ENABLE_PYTHON_MODULE "Enable pybind11 bindings of library" OFF)

######################### Process & Validate Options ###########################
include(autocmake_safeguards)
#include(autocmake_omp)
include(custom_static_library)

#find_package(TargetLAPACK REQUIRED)
find_package(libefp REQUIRED)

#  <<  Pybind11 & Python  >>
set(PYBIND11_CPP_STANDARD "-std=c++11") #${CMAKE_CXX_STANDARD}")
find_package(pybind11 2.0.0 CONFIG REQUIRED)
message(STATUS "${Cyan}Using pybind11${ColourReset}: ${pybind11_INCLUDE_DIR} (version ${pybind11_VERSION} for Py${PYTHON_VERSION_STRING} and ${PYBIND11_CPP_STANDARD})")
message(STATUS "${Cyan}Using Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE}")

################################# Main Project #################################
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(PN ${PROJECT_NAME})

# <<<  Build  >>>

#set(raw_sources_list ai_disp.c balance.c clapack.c disp.c efp.c elec.c
#                     elec_terms.c int.c log.c parse.c pol.c pol_direct.c
#                     stream.c swf.c util.c xr.c)
#set(src_prefix "src/")
#string(REGEX REPLACE "([^;]+)" "${src_prefix}\\1" sources_list "${raw_sources_list}")
#
## STATIC/SHARED on below governed by BUILD_SHARED_LIBS
#add_library(efp ${sources_list})
#set_target_properties(efp PROPERTIES POSITION_INDEPENDENT_CODE ${BUILD_FPIC}
#                                     COMPILE_FLAGS "-std=c99")
#set_source_files_properties(${src_prefix}int.c PROPERTIES COMPILE_FLAGS -O1)
#if(${BUILD_SHARED_LIBS})
#    target_link_libraries(efp PRIVATE ${LIBC_INTERJECT})
#    if(APPLE)
#        set_target_properties(efp PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
#    endif()
#endif()
#target_link_libraries(efp PRIVATE tgt::lapack)

# <<<  Python Module  >>>

pybind11_add_module(pylibefp pylibefp/core.cc) # ${source_list})
#    target_include_directories(pylibefp PRIVATE ${src_prefix})
target_include_directories(pylibefp PRIVATE /home/psilocaluser/gits/libefp/src)
target_link_libraries(pylibefp PRIVATE efp)
    #add_library(pylibefp SHARED pylibefp/pylibefp.cc ${sources_list})
    #target_link_libraries(pylibefp PUBLIC pybind11::module)



#set_target_properties(pylibefp PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"  # for python module
#                               OUTPUT_NAME pylibefp
#                               EXPORT_NAME pylibefp
#                               SUFFIX ".so")  # for python module
                               ##SUFFIX "${PYTHON_MODULE_EXTENSION}" # for python module
                               #INSTALL_RPATH "${psi4_RPATH}"
                               #BUILD_WITH_INSTALL_RPATH ON)

# <<<  Install  >>>

#install(FILES fraglib/makefp.inp
#        DESTINATION ${CMAKE_INSTALL_DATADIR}/${PN})
#install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/fraglib
#        DESTINATION ${CMAKE_INSTALL_DATADIR}/${PN})
#
## headers NOT namespace protected
#install(FILES ${src_prefix}/efp.h
#        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#install(TARGETS efp
#        EXPORT "${PN}Targets"
#        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
#
## <<<  Export Interface  >>>
#
#target_compile_definitions(efp INTERFACE USING_${PN})
#target_include_directories(efp INTERFACE
#                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
#
## <<<  Export Config  >>>
#
## GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
#set(CMAKECONFIG_INSTALL_DIR "share/cmake/${PN}")
#configure_package_config_file(cmake/${PN}Config.cmake.in
#                              "${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake"
#                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
#write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
#                                 VERSION ${${PN}_VERSION}
#                                 COMPATIBILITY AnyNewerVersion)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PN}Config.cmake
#              ${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake
#        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
#install(EXPORT "${PN}Targets"
#        NAMESPACE "${PN}::"
#        DESTINATION ${CMAKECONFIG_INSTALL_DIR})
