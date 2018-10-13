#
# @BEGIN LICENSE
#
#   pylibefp/extras.py:
#
#   Copyright (c) 2017-2018 The Psi4 Developers
#
#   All rights reserved. Use of this source code is governed by a
#   BSD-style license that can be found in the LICENSE file.
#
# @END LICENSE
#

import os

_success_flag_ = False


###   # Working directory
###   _input_dir_ = os.getcwd()
###   
###   def get_input_directory():
###       return _input_dir_

# Testing
def test(extent='full', extras=None):
    """Runs a test suite through pytest.

    Parameters
    ----------
    extent : {'smoke', 'quick', 'full', 'long'}
        All choices are defined, but choices may be redundant in some projects.
        _smoke_ will be minimal "is-working?" test(s).
        _quick_ will be as much coverage as can be got quickly, approx. 1/3 tests.
        _full_ will be the whole test suite, less some exceedingly long outliers.
        _long_ will be the whole test suite.
    extras : list
        Additional arguments to pass to `pytest`.

    Returns
    -------
    int
        Return code from `pytest.main()`. 0 for pass, 1 for fail.

    """
    try:
        import pytest
    except ImportError:
        raise RuntimeError('Testing module `pytest` is not installed. Run `conda install pytest`')
    abs_test_dir = os.path.sep.join([os.path.abspath(os.path.dirname(__file__)), "tests"])

    command = ['-rws', '-v']
    if extent.lower() in ['smoke', 'quick']:
        command.extend(['-m', 'quick'])
    elif extent.lower() == 'full':
        command.extend(['-m', 'not long'])
    elif extent.lower() == 'long':
        pass
    if extras is not None:
        command.extend(extras)
    command.extend(['--capture=sys', abs_test_dir])

    retcode = pytest.main(command)
    return retcode
