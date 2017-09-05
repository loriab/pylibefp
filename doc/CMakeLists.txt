# Settle docs-specific dependencies
find_package(Sphinx QUIET)
if(NOT SPHINX_FOUND)
   message(STATUS "No Sphinx, no docs. Pre-built documentation at http://psicode.org/psi4manual/master/index.html")
endif()

execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
    "import sphinx-psi4-theme cloud_sptheme; import astropy_helpers; import python-graphviz"
    OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE _exec_python_err)
if(_exec_python_err)
    message(STATUS "No theme, no docs (try: conda install sphinx-psi4-theme cloud_sptheme astropy-helpers graphviz python-graphviz -c conda-forge -c astropy -c psi4)")
    set(_sphinx_stuffs "ON")
endif()

# Setup and build Sphinx documentation
if(SPHINX_FOUND AND _sphinx_stuffs)
    set(CCSD "${CMAKE_CURRENT_SOURCE_DIR}")
    set(CCBD "${CMAKE_CURRENT_BINARY_DIR}")

    # Static-doc reST files
    set(STATICDOC index.rst api.rst abbr_accents.rst 
    introduction.rst
    )

    # * compute relative path btwn top_srcdir and objdir/doc/sphinxman
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
        "import sys; import os; print('/' + os.path.relpath('${PROJECT_SOURCE_DIR}', '${CCBD}') + '/../')"
        OUTPUT_VARIABLE SFNX_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
    # * add ipynb doc if build env can handle it
#    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c
#        "import nbsphinx; import jupyter_client; import ipykernel"
#        OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE _exec_python_err)
#    if(_exec_python_err)
#        message(STATUS "Some ipynb docs require extra modules (try: conda install nbsphinx jupyter_client ipykernel -c conda-forge)")
#        list(APPEND STATICDOC "psiapi.rst")
#        set(_jupy "")
#        set(_jupyconf "")
#    else()
#        list(APPEND STATICDOC "psiapi.ipynb")
#        set(_jupy "w/Jupyter ")
#        set(_jupyconf "'nbsphinx',")
#    endif()
    # configure Sphinx build
    configure_file(
        "${CCSD}/source/conf.py.in"
        "${CCBD}/source/conf.py"
        @ONLY)
    # * copy all over to objdir, substituting literalinclude along the way
    foreach(rst ${STATICDOC})
       configure_file(
        ${CCSD}/source/${rst}
        ${CCBD}/source/${rst}
        @ONLY)
    endforeach()
    # * duplicate index.rst as contents.rst for toc link
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CCBD}/source/index.rst ${CCBD}/source/contents.rst)

    # Static resource png and pdf files
    set(STATICRES VFPt_octupole_wires_py.png)
    # * copy image and pdf files all over to objdir (configure_file corrupts file)
    foreach(fl ${STATICRES})
        execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CCSD}/source/resources/${fl} ${CCBD}/source)
    endforeach()

    # Autodoc efp fragments
#    if(${ENABLE_libefp})
#        add_custom_command(
#            COMMAND ${PYTHON_EXECUTABLE} ${CCSD}/document_efpfrag.py ${CCSD}
#            OUTPUT  ${CCBD}/source/autodoc_available_efpfrag.rst
#            DEPENDS ${CCSD}/document_efpfrag.py
#                    ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/psi4/efpfrag/*efp
#            COMMENT "Autodoc efp fragments")
#    else()
#        file(WRITE ${CCBD}/source/autodoc_available_efpfrag.rst "")
#    endif()


    # Establishing dependencies for docs targets

    # * announcement
    message(STATUS "Documentation targets ${_jupy}available: sphinxman (html)")

    # * full, proper build
    add_custom_target(sphinxman
        DEPENDS #${CCBD}/source/autodoc_psi4bases.rst
                #${CCBD}/source/autodoc_available_efpfrag.rst
        COMMENT "Preparing Sphinx HTML documentation build")
    add_dependencies(sphinxman pylibefp-core)

    install(
        DIRECTORY ${CCBD}/html
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DOCDIR}
        OPTIONAL)

    # Establishing commands to build docs targets

    # * full, proper build with html builder
    add_custom_command(
        TARGET sphinxman POST_BUILD
        COMMAND ${SPHINX_EXECUTABLE}
            # -v  # Repeat for greater verbosity
            # -q  # Enable for quiet mode
            # -Q  # Enable for even more quiet mode
            # -E  # don't use a saved environment, always read all file
            # -a  # write all files; default is to only write new and changed files
            -b html               # builder
            -d ${CCBD}/_doctrees  # cache w/ pickled reST
            -c ${CCBD}/source     # whereabouts of conf.py
            ${CCBD}/source        # build from
            ${CCBD}/html          # build to
        COMMENT "Building html documentation ...")

endif()

