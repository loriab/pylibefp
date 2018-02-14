#
# @BEGIN LICENSE
#
#   pylibefp/__init__.py:
#
#   Copyright (c) 2017 The Psi4 Developers
#
#   All rights reserved. Use of this source code is governed by a
#   BSD-style license that can be found in the LICENSE file.
#
# @END LICENSE
#


# Figure out psidatadir: envvar trumps staged/installed
import os
pylibefp_module_loc = os.path.dirname(os.path.abspath(__file__))
#pymod = os.path.normpath(os.path.sep.join(['@PYMOD_INSTALL_LIBDIR@', '@CMAKE_INSTALL_LIBDIR@', 'psi4']))
#if pymod.startswith(os.path.sep + os.path.sep):
#    pymod = pymod[1:]
#pymod_dir_step = os.path.sep.join(['..'] * pymod.count(os.path.sep))
#data_dir = os.path.sep.join([psi4_module_loc, pymod_dir_step, '@CMAKE_INSTALL_DATADIR@', 'psi4'])
#data_dir = '@libefp_FRAGLIB_DIRS@'

#if "EFPDATADIR" in os.environ.keys():
#    data_dir = os.path.expanduser(os.environ["EFPDATADIR"])
#elif "CMAKE_INSTALL_DATADIR" in data_dir:
#    data_dir = os.path.sep.join([os.path.abspath(os.path.dirname(__file__)), "share", "psi4"])
#print("efp_frags {} {}".format('@libefp_FRAGLIB_DIRS@', data_dir))

# Init core
from . import core
#print("__init__ core    {}".format(core.__file__))

#from psi4.core import set_output_file, get_variable, set_variable, get_num_threads, set_num_threads
#core.initialize()
#core.efp_init()

# Load driver and outfile paraphernalia
from .wrapper import *
#print("__init__ core    {}".format(core.__file__))
#print("__init__ wrapper {}".format(wrapper.__file__))
from .metadata import __version__, version_formatter

# A few extraneous functions
from .extras import test
