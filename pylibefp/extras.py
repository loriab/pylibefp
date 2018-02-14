#
# @BEGIN LICENSE
#
#   pylibefp/extras.py:
#
#   Copyright (c) 2017 The Psi4 Developers
#
#   All rights reserved. Use of this source code is governed by a
#   BSD-style license that can be found in the LICENSE file.
#
# @END LICENSE
#

import os
import atexit
import subprocess

from . import core

###   # Numpy place holder for files and cleanup
###   numpy_files = []
###   def register_numpy_file(filename):
###       if filename not in numpy_files:
###           numpy_files.append(filename)
###   
###   def clean_numpy_files():
###       for nfile in numpy_files:
###           os.unlink(nfile)
###   
###   atexit.register(clean_numpy_files)
###   
###   # Exit printing
###   def exit_printing():
###       if _success_flag_:
###           core.print_out( "\n*** Psi4 exiting successfully. Buy a developer a beer!\n")
###       else:
###           core.print_out( "\n*** Psi4 encountered an error. Buy a developer more coffee!\n")
###           core.print_out( "*** Resources and help at github.com/psi4/psi4.\n")

_success_flag_ = False


###   # Working directory
###   _input_dir_ = os.getcwd()
###   
###   def get_input_directory():
###       return _input_dir_


# Testing
def test():
    """Runs a smoke test suite through pytest."""

    try:
        import pytest
    except ImportError:
        raise RuntimeError('Testing module `pytest` is not installed. Run `conda install pytest`')
    abs_test_dir = os.path.sep.join([os.path.abspath(os.path.dirname(__file__)), "tests"])
    pytest.main(['-rws', '-v', '--capture=sys', abs_test_dir])
