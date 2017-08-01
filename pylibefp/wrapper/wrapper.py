#
# @BEGIN LICENSE
#
#   pylibefp/wrapper/wrapper.py:
#
#   Copyright (c) 2017 The Psi4 Developers
#
#   All rights reserved. Use of this source code is governed by a
#   BSD-style license that can be found in the LICENSE file.
#
# @END LICENSE
#
from __future__ import print_function

import os
from pylibefp import core
from . import psiutil
from .exceptions import *

try:
    basestring
except NameError:
    basestring = str


def _pywrapped_result_to_error(res, msg=''):

    if res == core.efp_result.EFP_RESULT_SUCCESS:
        return
    elif res == core.efp_result.EFP_RESULT_FATAL:
        raise Fatal(msg)
    elif res == core.efp_result.EFP_RESULT_NO_MEMORY:
        raise NoMemory(msg)
    elif res == core.efp_result.EFP_RESULT_FILE_NOT_FOUND:
        raise FileNotFound(msg)
    elif res == core.efp_result.EFP_RESULT_SYNTAX_ERROR:
        raise SyntaxError(msg)
    elif res == core.efp_result.EFP_RESULT_UNKNOWN_FRAGMENT:
        raise UnknownFragment(msg)
    elif res == core.efp_result.EFP_RESULT_POL_NOT_CONVERGED:
        raise PolNotConverged(msg)


def _pywrapped_efp_prepare(efpobj):
    res = efpobj.raw_prepare()
    _pywrapped_result_to_error(res)


def _pywrapped_efp_compute(efpobj):
    res = efpobj.raw_compute()
    _pywrapped_result_to_error(res)


def _pywrapped_add_potential(efpobj, potential, fragpath='LIBRARY', duplicates_ok=False):
    """Searches for EFP fragments and adds to *efpobj*.

    *potential* is a single fragment name or a list of fragments,
    with or without `.efp` extension.

    *fragpath* should be a string with :-separated paths that may include
    absolute paths, relative paths, $-marked environment variables,
    and the word LIBRARY, all of which will be expanded with LIBRARY
    expanded to the native libefp fragment library.

    *duplicates_ok* is a boolean indicating whether to return
    pylibefp.wrapper.exceptions.Fatal if asked to load a duplicate
    potential according to libefp default behavior (filtered where
    possible) or to continue.

    """
    # form unified path list for efpfrags
    paths = []
    for pth in fragpath.split(os.pathsep):
        if pth == 'LIBRARY':
            paths.append('@libefp_FRAGLIB_DIRS@')
        else:
            paths.append(os.path.expandvars(os.path.expanduser(pth)))
    paths = os.pathsep.join(paths)

    #        paths.append(libraryPath)
    #    elif pth in os.environ:
    #        envvar = os.environ.get('PSIDATADIR', None)
    #        paths.extend(envvar.split(os.pathsep))
    #    else:
    #        paths.extend(pth.split(os.pathsep))

    # locate efpfrags full path name
    abspath_pots = []
    if isinstance(potential, basestring):
        potential = [potential]
    uniq_pots = list(set(potential))
    for pot in uniq_pots:
        if not pot.endswith('.efp'):
            pot += '.efp'
        abspath_pots.append(psiutil.search_file(pot, paths))

    # load the potentials
    for ipot, pot in enumerate(abspath_pots):
        res = efpobj.add_potential(pot)
        try:
            _pywrapped_result_to_error(res, uniq_pots[ipot])
        except Fatal as e:
            if duplicates_ok:
                pass
            else:
                raise

        print(r"""  EFP fragment {} read from {}""".format(uniq_pots[ipot], pot))


def _pywrapped_add_fragment(efpobj, fragments):
    """Registers EFP fragments to *efpobj* in order.

    """
    if isinstance(fragments, basestring):
        fragments = [fragments]
    for frag in fragments:
        res = efpobj.add_fragment(frag)
        _pywrapped_result_to_error(res, frag)


def opts_summary(efpobj, labels='libefp'):

    opts = efpobj.get_opts()
    print(opts)
#    py::dict opts = opts_to_dict(efp);

    elec = 'Electrostatics'
    xr = 'Exchange'
    disp = 'Dispersion'
    if labels == 'libefp':
        pol = 'Polarization'
    elif labels == 'psi4':
        pol = 'Induction'

    text = ''
    text += '\n  ==> EFP/EFP Setup <==\n\n'
    text +=   '  {:<30} {:12}\n'.format(elec + ' enabled?:', 'true' if opts["elec"] else 'false')
    text +=   '  {:<30} {:12}\n'.format(pol + ' enabled?:', 'true' if opts["pol"] else 'false')
    text +=   '  {:<30} {:12}\n'.format(disp + ' enabled?:', 'true' if opts["disp"] else 'false')
    text +=   '  {:<30} {:12}\n'.format(xr + ' enabled?:', 'true' if opts["xr"] else 'false')
    text +=   '  {:<30} {:12}\n'.format('Charge-Transfer enabled?:', 'undefined')

    text += '\n  {:<30} {:12}\n'.format(elec + ' damping:', opts["elec_damp"])
    text +=   '  {:<30} {:12}\n'.format(pol + 'damping:', opts["pol_damp"])
    text +=   '  {:<30} {:12}\n'.format(disp + ' damping:', opts["disp_damp"])
    text +=   '  {:<30} {:12}\n'.format(pol + ' driver:', opts["pol_driver"])

    text += '\n  ==> QM/EFP Setup <==\n\n'
#//    sprintf(buffer, "  Number of QM fragments:  %12d\n", -1); //, nfrag_);
    text += '  Electrostatics enabled?:   {:12}\n'.format('true' if opts["ai_elec"] else 'false')
    text += '  Polarization enabled?:     {:12}\n'.format('true' if opts["ai_pol"] else 'false')
    text += '  Dispersion enabled?:       {:12}\n'.format('undefined')
    text += '  Exchange enabled?:         {:12}\n'.format('undefined')
    text += '  Charge-Transfer enabled?:  {:12}\n'.format('undefined')

    return text


def energy_summary(efpobj):

    opt = efpobj.get_opts()
    ene = efpobj.get_energy()

    print('ENERGY SUMM', opt, ene)

#{'pol_damp': 'tt', 'elec_damp': 'screen', 'enable_cutoff': 0, 'enable_pbc': 0, 'pol_driver': 'iterative',  'pol': False, 'swf_cutoff': 0.0, 'ai_pol': False, 'disp_damp': 'overlap', 'disp': False, 'xr': False

    def _enabled(tf, t='*', f=''):
        if tf:
            return t
        else:
            return f

    text = ''
    text += '\n'
    text += '\n    EFP Results\n'
    text +=   '  ------------------------------------------------------------\n'
    text +=   '    Electrostatics                {:20.12f} [Eh] {}\n'.format(
                                      ene['electrostatic'] +
                                      ene['charge_penetration'] +
                                      ene['electrostatic_point_charges'],
                                      _enabled(opt['elec'] or opt['ai_elec']))
    text +=   '      EFP/EFP                     {:20.12f} [Eh] {}\n'.format(
                                      ene['electrostatic'] + 
                                      ene['charge_penetration'],
                                      _enabled(opt['elec']))
    text +=   '      QM-Nuc/EFP                  {:20.12f} [Eh] {}\n'.format(
                                      ene['electrostatic_point_charges'],
                                      _enabled(opt['ai_elec']))
    text += '\n    Exchange                      {:20.12f} [Eh] {}\n'.format(
                                      ene['exchange_repulsion'],
                                      _enabled(opt['xr']))
    text +=   '      EFP/EFP                     {:20.12f} [Eh] {}\n'.format(
                                      ene['exchange_repulsion'],
                                      _enabled(opt['xr']))
    text +=   '      QM/EFP                      {:20.12f} [Eh] {}\n'.format(
                                      0.0,
                                      '')
    text += '\n    Induction                     {:20.12f} [Eh] {}\n'.format(
                                      ene['polarization'],
                                      _enabled(opt['pol'] or opt['ai_pol']))
    text +=   '      {:7}                     {:20.12f} [Eh] {}\n'.format(
                                      _enabled(opt['ai_pol'], t='QM/EFP', f='EFP/EFP'),
                                      ene['polarization'],
                                      _enabled(opt['pol'] or opt['ai_pol']))
    text += '\n    Dispersion                    {:20.12f} [Eh] {}\n'.format(
                                      ene['dispersion'],
                                      _enabled(opt['disp']))
    text +=   '      EFP/EFP                     {:20.12f} [Eh] {}\n'.format(
                                      ene['dispersion'],
                                      _enabled(opt['disp']))
    text +=   '      QM/EFP                      {:20.12f} [Eh] {}\n'.format(
                                      0.0,
                                      '')
    text += '\n    Total EFP                     {:20.12f} [Eh]\n'.format(
                                      ene['total'])

    return text


# only wrapped to throw Py exceptions
core.efp.prepare = _pywrapped_efp_prepare
core.efp.compute = _pywrapped_efp_compute

core.efp.add_potentials = _pywrapped_add_potential
core.efp.add_fragments = _pywrapped_add_fragment
core.efp.energy_summary = energy_summary
