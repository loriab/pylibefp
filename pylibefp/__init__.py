#
# @BEGIN LICENSE
#
#   pylibefp/__init__.py:
#
#   Copyright (c) 2017-2020 The Psi4 Developers
#
#   All rights reserved. Use of this source code is governed by a
#   BSD-style license that can be found in the LICENSE file.
#
# @END LICENSE
#

import os

# Init core
from . import core
from .exceptions import (
    EFPException,
    EFPSyntaxError,
    Fatal,
    FileNotFound,
    NoMemory,
    PolNotConverged,
    PyEFPSyntaxError,
    UnknownFragment,
)

# A few extraneous functions
from .extras import test
from .metadata import __version__, version_formatter

# Load driver and version paraphernalia
from .wrapper import extract_subsets, from_dict, to_dict

pylibefp_module_loc = os.path.dirname(os.path.abspath(__file__))
