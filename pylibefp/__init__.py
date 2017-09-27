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

#data_dir = os.path.abspath(data_dir)
#if not os.path.isdir(data_dir):
#    raise KeyError("Unable to read the Psi4 Python folder - check the PSIDATADIR environmental variable"
#                    "      Current value of PSIDATADIR is %s" % data_dir)
#os.environ["PSIDATADIR"] = data_dir

# Init core
#try:
from . import core
print("__init__ core    {}".format(core.__file__))
#except ImportError as err:
#    if 'CXXABI' in str(err):
#        raise ImportError("{0}\nLikely cause: GCC >= 4.9 not in [DY]LD_LIBRARY_PATH".format(err))
#    else:
#        raise ImportError("{0}".format(err))

#from psi4.core import set_output_file, get_variable, set_variable, get_num_threads, set_num_threads
#core.initialize()
#core.efp_init()
#
#if "PSI_SCRATCH" in os.environ.keys():
#    envvar_scratch = os.environ["PSI_SCRATCH"]
#    if not os.path.isdir(envvar_scratch):
#        raise Exception("Passed in scratch is not a directory (%s)." % envvar_scratch)
#    core.IOManager.shared_object().set_default_path(envvar_scratch)
#
## Cleanup core at exit
#import atexit
#atexit.register(core.set_legacy_molecule, None)
#atexit.register(core.clean_options)
#atexit.register(core.clean)
#atexit.register(core.finalize)
#
## Make official plugins accessible in input
#from .driver import endorsed_plugins
#
## Manage threads. Must be after endorsed plugins, honestly.
#core.set_num_threads(1, quiet=True)
#
## Load driver and outfile paraphernalia
from .wrapper import *
print("__init__ core    {}".format(core.__file__))
print("__init__ wrapper {}".format(wrapper.__file__))
#from .driver import *
#from .header import print_header
from .metadata import __version__, version_formatter

# A few extraneous functions
#from .extras import get_input_directory, addons, test
from .extras import test

## Python portions of compiled-in Add-Ons
#import sys
#if "@ENABLE_PCMSolver@".upper() in ["1", "ON", "YES", "TRUE", "Y"]:
#    sys.path.insert(1, "@PCMSolver_PYMOD@")
#
