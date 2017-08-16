#
# @BEGIN LICENSE
#
#   pylibefp/wrapper/psiutil.py:
#
#   Copyright (c) 2017 The Psi4 Developers
#
#   All rights reserved. Use of this source code is governed by a
#   BSD-style license that can be found in the LICENSE file.
#
# @END LICENSE
#

from __future__ import absolute_import
from __future__ import print_function
import os


## {{{ http://code.activestate.com/recipes/52224/#c4
def search_file(filename, search_path, implicitExt=''):
    """Given a search_pathsep-delimited search_path string, find filename.
    Returns search_path to filename if found, otherwise None.
    Also allows for files with implicit extensions (eg, .exe), but
    always returning filename as was provided.

    >>> search_file('ls', '/usr/bin:/bin', implicitExt='.exe')
    '/bin/ls'

    """

    if (os.path.isfile(filename)
          or implicitExt and os.path.isfile(filename + implicitExt)):
        # Already absolute path.
        return filename
    for p in search_path.split(os.pathsep):
        fp = os.path.abspath(os.path.expandvars(os.path.expanduser(p)))
        candidate = os.path.join(fp, filename)
        if (os.path.isfile(candidate)
              or implicitExt and os.path.isfile(candidate + implicitExt)):
            return candidate
    return None
## end of http://code.activestate.com/recipes/52224/#c4 }}}
