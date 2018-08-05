#
# @BEGIN LICENSE
#
#   pylibefp/__init__.py:
#
#   Copyright (c) 2017-2018 The Psi4 Developers
#
#   All rights reserved. Use of this source code is governed by a
#   BSD-style license that can be found in the LICENSE file.
#
# @END LICENSE
#

import os
pylibefp_module_loc = os.path.dirname(os.path.abspath(__file__))

# Init core
from . import core

# Load driver and version paraphernalia
from .wrapper import *
from .metadata import __version__, version_formatter

# A few extraneous functions
from .extras import test
