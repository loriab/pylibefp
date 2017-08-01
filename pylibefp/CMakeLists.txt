
# <<<  Build  >>>

pybind11_add_module(core core.cc)
target_include_directories(core PRIVATE /home/psilocaluser/gits/libefp/src)
target_link_libraries(core PRIVATE efp)

# <<<  Install  >>>

    # <<<  install bin/  >>>
#configure_file(run_pylibefp.py psi4 @ONLY)
#install(PROGRAMS ${CMAKE_BINARY_DIR}/psi4
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
#
#    # <<<  install lib/  >>>
#file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)
#
#install(DIRECTORY share/psi4/
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4
#        PATTERN "*.pyc" EXCLUDE)
#
#install(DIRECTORY include/psi4/
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}/psi4)
#
#install(DIRECTORY wrapper/
#        DESTINATION ${CMAKE_BINARY_DIR}
#        FILES_MATCHING PATTERN "*.py")
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4
                       #EXCLUDE wrapper.py)

#install(FILES header.py
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)

configure_file(__init__.py __init__.py COPYONLY)
configure_file(extras.py extras.py COPYONLY)
configure_file(wrapper/__init__.py wrapper/__init__.py COPYONLY)
configure_file(wrapper/wrapper.py wrapper/wrapper.py @ONLY)
configure_file(wrapper/psiutil.py wrapper/psiutil.py COPYONLY)
configure_file(wrapper/exceptions.py wrapper/exceptions.py COPYONLY)

install(FILES ${CMAKE_BINARY_DIR}/__init__.py
              ${CMAKE_BINARY_DIR}/extras.py
              ${CMAKE_BINARY_DIR}/wrapper/wrapper.py
        DESTINATION ${CMAKE_BINARY_DIR})
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}/psi4)

