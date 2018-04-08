cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(pylibefp-core 
        LANGUAGES C CXX)

list(APPEND CMAKE_MODULE_PATH ${PYLIBEFP_ROOT}/cmake)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(custom_color_messages)
# GNUInstallDirs "DATADIR" wrong here; CMake search path wants "share".
set(CMAKECONFIG_INSTALL_DIR "share/cmake/pylibefp")

# <<<  Marshal Dependencies & Add-ons  >>>

# TODO put in libefpConfig instead?
#find_package(TargetLAPACK CONFIG REQUIRED)
#get_property(_ill TARGET tgt::lapk PROPERTY INTERFACE_LINK_LIBRARIES)
#list(GET _ill 0 _ill0)
#message(STATUS "${Cyan}Using LAPACK${_isMKL}${ColourReset}: ${_ill0};...")

#  <<  Pybind11 & Python  >>
set(PYBIND11_CPP_STANDARD "-std=c++${CMAKE_CXX_STANDARD}")
find_package(pybind11 2.2.2 EXACT CONFIG REQUIRED)
message(STATUS "${Cyan}Using pybind11${ColourReset}: ${pybind11_INCLUDE_DIR} (version ${pybind11_VERSION} for Py${PYTHON_VERSION_STRING} and ${PYBIND11_CPP_STANDARD})")
message(STATUS "${Cyan}Using Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}${ColourReset}: ${PYTHON_EXECUTABLE}")

set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)

find_package(libefp 1.5.0 CONFIG REQUIRED COMPONENTS shallow)
get_property(_loc TARGET libefp::efp PROPERTY LOCATION)
list(APPEND _addons ${_loc})
message(STATUS "${Cyan}Using libefp${ColourReset}: ${_loc} (version ${libefp_VERSION})")

# <<<  Build  >>>

pybind11_add_module(core NO_EXTRAS core.cc)
target_link_libraries(core PRIVATE libefp::efp)
target_link_libraries(core PRIVATE ${LIBC_INTERJECT})
target_link_libraries(core PUBLIC pybind11::module)
target_link_libraries(core PRIVATE Threads::Threads)

set(PN pylibefp)
install(TARGETS core
        EXPORT "${PN}Targets"
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pylibefp)
install(EXPORT "${PN}Targets"
        NAMESPACE "${PN}::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

set_target_properties(core PROPERTIES PREFIX "${PYTHON_MODULE_PREFIX}"  # for python module
                           OUTPUT_NAME core
                           EXPORT_NAME core
                           SUFFIX ".so"  # for python module
                           #SUFFIX "${PYTHON_MODULE_EXTENSION}" # for python module
                           INSTALL_RPATH "${psi4_RPATH}"
                           BUILD_WITH_INSTALL_RPATH ON)

# <<<  Version  >>>

add_custom_target(update_version ALL
                  COMMAND ${PYTHON_EXECUTABLE} versioner.py --metaout ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
                                                            --cmakeout ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
                  COMMAND cmake -DWTO="${CMAKE_CURRENT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}"
                                -DPN="pylibefp"
                                -P ${CMAKE_CURRENT_BINARY_DIR}/metadata.cmake
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  COMMENT "Generating version info")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/metadata.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pylibefp)

# <<<  Install  >>>

    # <<<  install lib/  >>>
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pylibefp)

install(FILES __init__.py extras.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pylibefp)

install(FILES ../tests/utils.py
              ../tests/addons.py
              ../tests/conftest.py
              ../tests/systems.py
              ../tests/test_efpefp.py
              ../tests/test_opts.py
              ../tests/test_coverage.py
              ../tests/test_dict.py
              ../tests/test_scf.py
              ../tests/test_psi.py
              ../tests/test_libefp.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pylibefp/tests/)

install(DIRECTORY wrapper
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pylibefp
        FILES_MATCHING PATTERN "*.py")

configure_file(wrapper/wrapper.py wrapper/wrapper.py @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/wrapper/wrapper.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/pylibefp/wrapper)

# <<<  Export Config  >>>

configure_package_config_file(
        ${PN}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/${PN}Config.cmake
        INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/${PN}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${CMAKECONFIG_INSTALL_DIR}/${PN}ConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR})

