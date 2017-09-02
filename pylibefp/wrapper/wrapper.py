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
import re
import math
from pylibefp import core
from . import psiutil
from .exceptions import *

try:
    basestring
except NameError:
    basestring = str


_lbtl = {
    'libefp': {},
    'psi': {
        'elec': 'elst',
        'pol': 'ind',
        'xr': 'exch',
        'elec_damp': 'elst_damp',
        'pol_damp': 'ind_damp',
        'pol_driver': 'ind_driver',
        'ai_elec': 'ai_elst',
        'ai_pol': 'ai_ind',
        'ai_xr': 'ai_exch',
    },
}


def _rekey(rawdict, label):
    newdict = rawdict.copy()
    for key in rawdict.keys():
        topic = _lbtl[label].get(key, key)
        newdict[topic] = newdict.pop(key)
    return newdict


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
        raise EFPSyntaxError(msg)
    elif res == core.efp_result.EFP_RESULT_UNKNOWN_FRAGMENT:
        raise UnknownFragment(msg)
    elif res == core.efp_result.EFP_RESULT_POL_NOT_CONVERGED:
        raise PolNotConverged(msg)


def _pywrapped_efp_prepare(efpobj):
    res = efpobj.raw_prepare()
    _pywrapped_result_to_error(res)


def _pywrapped_efp_compute(efpobj, do_gradient=False):
    res = efpobj.raw_compute(do_gradient=int(do_gradient))
    _pywrapped_result_to_error(res)


def _pywrapped_add_potential(efpobj, potential, fragpath='LIBRARY', duplicates_ok=False):
    """Searches for EFP fragments and adds to *efpobj*.

    Parameters
    ----------
    potential : str or list
        Single fragment name or a list of fragments, with or without
        ".efp" extension.
    fragpath : string, optional
        String with ":"-separated paths that may include absolute
        paths, relative paths, "$"-marked environment variables, and
        the word LIBRARY. Relative paths and environment variables will
        be expanded. "LIBRARY" will be expanded to the native libefp
        fragment library.
    duplicates_ok : bool, optional
        Whether to continue or to return pylibefp.Fatal if asked
        to load a duplicate potential according to libefp default
        behavior. The `potential` list is always filtered to avoid
        duplicates. Activating `duplicates_ok` additionally allows
        repeated calls to this function to add duplicate potentials.

    Returns
    -------
    None

    """
    # TODO handle casing of potential
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
        res = efpobj.raw_add_potential(pot)
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
        res = efpobj.raw_add_fragment(frag)
        _pywrapped_result_to_error(res, frag)


def _pywrapped_get_opts(efpobj, label='libefp'):
    """Returns the options state of *efpobj* as a dictionary.

    Parameters
    ----------
    label : str, optional
        Returned dictionary keys are identical to libefp efp_opts struct
        names unless custom renaming requested via `label`.

    Returns
    -------
    dict
        Current options state of `efpobj` translated into bools, strings,
        and floats, rather than libefp custom datatypes.

    """
    opts = core.efp_opts()
    res = efpobj.raw_get_opts(opts)
    _pywrapped_result_to_error(res)

    dopts = {}

    dopts['elec'] = bool(opts.terms & core.efp_term.EFP_TERM_ELEC)
    dopts['pol'] = bool(opts.terms & core.efp_term.EFP_TERM_POL)
    dopts['disp'] = bool(opts.terms & core.efp_term.EFP_TERM_DISP)
    dopts['xr'] = bool(opts.terms & core.efp_term.EFP_TERM_XR)
    dopts['chtr'] = bool(opts.terms & core.efp_term.EFP_TERM_CHTR)

    dopts['elec_damp'] = {
            core.EFP_ELEC_DAMP_SCREEN: 'screen',
            core.EFP_ELEC_DAMP_OVERLAP: 'overlap',
            core.EFP_ELEC_DAMP_OFF: 'off',
        }[opts.elec_damp]

    dopts['pol_damp'] = {
               core.EFP_POL_DAMP_TT: 'tt',
               core.EFP_POL_DAMP_OFF: 'off',
        }[opts.pol_damp]

    dopts['disp_damp'] = {
               core.EFP_DISP_DAMP_TT: 'tt',
               core.EFP_DISP_DAMP_OVERLAP: 'overlap',
               core.EFP_DISP_DAMP_OFF: 'off',
        }[opts.disp_damp]

    dopts['enable_pbc'] = bool(opts.enable_pbc)
    dopts['enable_cutoff'] = bool(opts.enable_cutoff)
    dopts['swf_cutoff'] = opts.swf_cutoff

    dopts['pol_driver'] = {
               core.EFP_POL_DRIVER_ITERATIVE: 'iterative',
               core.EFP_POL_DRIVER_DIRECT: 'direct',
        }[opts.pol_driver]

    dopts['ai_elec'] = bool(opts.terms & core.efp_term.EFP_TERM_AI_ELEC)
    dopts['ai_pol'] = bool(opts.terms & core.efp_term.EFP_TERM_AI_POL)
    dopts['ai_disp'] = bool(opts.terms & core.efp_term.EFP_TERM_AI_DISP)
    dopts['ai_xr'] = bool(opts.terms & core.efp_term.EFP_TERM_AI_XR)
    dopts['ai_chtr'] = bool(opts.terms & core.efp_term.EFP_TERM_AI_CHTR)

    return _rekey(dopts, label=label)


def _pywrapped_set_opts(efpobj, dopts, label='libefp', append='libefp'):
    """Sets the options state of *efpobj* from dictionary `dopts`.

    Parameters
    ----------
    dopts : dict
        Input dict with keys from libefp efp_opts (see `label`) and
        values bools, strings, floats, and ints, as appropriate, rather
        than libefp custom datatypes.
    label : str, optional
        Input `dopts` keys are read as libefp efp_opts struct names or
        by the custom translation set defined for `label`.
    append : str, optional
        When 'libefp', input `dopts` keys are applied to the default
        (generally OFF) efp_opts state. When 'psi', input `dopts`
        keys are applied to the default (generally ON) Psi efp_opts
        state. When 'append', input `dopts` keys are applied to the
        current *efpobj* opt_opts state.

    Returns
    -------
    dict
        After setting the options state, `efpobj` is queried as to the
        current options state, which is then returned.

    """
    # warn on stray dopts keys
    allowed = ['elec', 'pol', 'disp', 'xr', 'elec_damp', 'pol_damp', 'disp_damp',
               'enable_pbc', 'enable_cutoff', 'swf_cutoff', 'pol_driver',
               'ai_elec', 'ai_pol']
    label_allowed = [_lbtl[label].get(itm, itm) for itm in allowed]
    for key in dopts.keys():
        if key not in label_allowed:
            print('Warning: unrecognized key {}'.format(key))

    # prepare base options state for dopts
    opts = core.efp_opts()
    if append == 'libefp':
        pass
    elif append == 'psi':
        opts.terms |= core.efp_term.EFP_TERM_ELEC
        opts.terms |= core.efp_term.EFP_TERM_POL
        opts.terms |= core.efp_term.EFP_TERM_DISP
        opts.terms |= core.efp_term.EFP_TERM_XR
        opts.elec_damp = core.EFP_ELEC_DAMP_SCREEN
        opts.pol_damp = core.EFP_POL_DAMP_TT
        opts.disp_damp = core.EFP_DISP_DAMP_OVERLAP
        opts.terms |= core.efp_term.EFP_TERM_AI_ELEC
        opts.terms |= core.efp_term.EFP_TERM_AI_POL
    #elif append == 'qchem':
    #    # q-chem and psi4 have different defaults for at least this option
    #    opts.disp_damp = core.EFP_DISP_DAMP_TT
    elif append == 'append':
        res = efpobj.raw_get_opts(opts)
        _pywrapped_result_to_error(res)
    else:
        raise PyEFPSyntaxError('Unrecognized opts default set: {}'.format(append))

    # apply dopts to options state
    topic = _lbtl[label].get('elec', 'elec')
    if topic in dopts:
        if dopts[topic] is True:
            opts.terms |= core.efp_term.EFP_TERM_ELEC
        elif dopts[topic] is False:
            opts.terms &= ~core.efp_term.EFP_TERM_ELEC
        else:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [T/F] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('pol', 'pol')
    if topic in dopts:
        if dopts[topic] is True:
            opts.terms |= core.efp_term.EFP_TERM_POL
        elif dopts[topic] is False:
            opts.terms &= ~core.efp_term.EFP_TERM_POL
        else:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [T/F] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('disp', 'disp')
    if topic in dopts:
        if dopts[topic] is True:
            opts.terms |= core.efp_term.EFP_TERM_DISP
        elif dopts[topic] is False:
            opts.terms &= ~core.efp_term.EFP_TERM_DISP
        else:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [T/F] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('xr', 'xr')
    if topic in dopts:
        if dopts[topic] is True:
            opts.terms |= core.efp_term.EFP_TERM_XR
        elif dopts[topic] is False:
            opts.terms &= ~core.efp_term.EFP_TERM_XR
        else:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [T/F] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('chtr', 'chtr')  # may be enabled in a future libefp release

    topic = _lbtl[label].get('elec_damp', 'elec_damp')
    if topic in dopts:
        try:
            opts.elec_damp = {
                    'screen': core.EFP_ELEC_DAMP_SCREEN,
                    'overlap': core.EFP_ELEC_DAMP_OVERLAP,
                    'off': core.EFP_ELEC_DAMP_OFF,
                }[dopts[topic].lower()]
        except KeyError:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [screen/overlap/off] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('pol_damp', 'pol_damp')
    if topic in dopts:
        try:
            opts.pol_damp = {
                    'tt': core.EFP_POL_DAMP_TT,
                    'off': core.EFP_POL_DAMP_OFF,
                }[dopts[topic].lower()]
        except KeyError:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [tt/off] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('disp_damp', 'disp_damp')
    if topic in dopts:
        try:
            opts.disp_damp = {
                    'overlap': core.EFP_DISP_DAMP_OVERLAP,
                    'tt': core.EFP_DISP_DAMP_TT,
                    'off': core.EFP_DISP_DAMP_OFF,
                }[dopts[topic].lower()]
        except KeyError:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [overlap/tt/off] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('enable_pbc', 'enable_pbc')
    if topic in dopts:
        if dopts[topic] is True:
            opts.enable_pbc = 1
        elif dopts[topic] is False:
            opts.enable_pbc = 0
        else:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [T/F] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('enable_cutoff', 'enable_cutoff')
    if topic in dopts:
        if dopts[topic] is True:
            opts.enable_cutoff = 1
        elif dopts[topic] is False:
            opts.enable_cutoff = 0
        else:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [T/F] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('swf_cutoff', 'swf_cutoff')
    if topic in dopts:
        opts.swf_cutoff = float(dopts[topic])

    topic = _lbtl[label].get('pol_driver', 'pol_driver')
    if topic in dopts:
        try:
            opts.pol_driver = {
                    'iterative': core.EFP_POL_DRIVER_ITERATIVE,
                    'direct': core.EFP_POL_DRIVER_DIRECT,
                }[dopts[topic].lower()]
        except KeyError:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [iterative/direct] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('ai_elec', 'ai_elec')
    if topic in dopts:
        if dopts[topic] is True:
            opts.terms |= core.efp_term.EFP_TERM_AI_ELEC
        elif dopts[topic] is False:
            opts.terms &= ~core.efp_term.EFP_TERM_AI_ELEC
        else:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [T/F] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('ai_pol', 'ai_pol')
    if topic in dopts:
        if dopts[topic] is True:
            opts.terms |= core.efp_term.EFP_TERM_AI_POL
        elif dopts[topic] is False:
            opts.terms &= ~core.efp_term.EFP_TERM_AI_POL
        else:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [T/F] {}: {}'.format(topic, dopts[topic]))

    topic = _lbtl[label].get('ai_damp', 'ai_damp')  # may be enabled in a future libefp release
    topic = _lbtl[label].get('ai_xr', 'ai_xr')      # may be enabled in a future libefp release
    topic = _lbtl[label].get('ai_chtr', 'ai_chtr')  # may be enabled in a future libefp release

    # set computed options state
    res = efpobj.raw_set_opts(opts)
    _pywrapped_result_to_error(res)

    return efpobj.get_opts(label=label)



#def opts_summary(efpobj, labels='libefp'):
#
#    opts = efpobj.get_opts()
##    py::dict opts = opts_to_dict(efp);
#
#    elec = 'Electrostatics'
#    xr = 'Exchange'
#    disp = 'Dispersion'
#    if labels == 'libefp':
#        pol = 'Polarization'
#    elif labels == 'psi4':
#        pol = 'Induction'
#
#    text = ''
#    text += '\n  ==> EFP/EFP Setup <==\n\n'
#    text +=   '  {:<30} {:12}\n'.format(elec + ' enabled?:', 'true' if opts["elec"] else 'false')
#    text +=   '  {:<30} {:12}\n'.format(pol + ' enabled?:', 'true' if opts["pol"] else 'false')
#    text +=   '  {:<30} {:12}\n'.format(disp + ' enabled?:', 'true' if opts["disp"] else 'false')
#    text +=   '  {:<30} {:12}\n'.format(xr + ' enabled?:', 'true' if opts["xr"] else 'false')
#    text +=   '  {:<30} {:12}\n'.format('Charge-Transfer enabled?:', 'undefined')
#
#    text += '\n  {:<30} {:12}\n'.format(elec + ' damping:', opts["elec_damp"])
#    text +=   '  {:<30} {:12}\n'.format(pol + 'damping:', opts["pol_damp"])
#    text +=   '  {:<30} {:12}\n'.format(disp + ' damping:', opts["disp_damp"])
#    text +=   '  {:<30} {:12}\n'.format(pol + ' driver:', opts["pol_driver"])
#
#    text += '\n  ==> QM/EFP Setup <==\n\n'
##//    sprintf(buffer, "  Number of QM fragments:  %12d\n", -1); //, nfrag_);
#    text += '  Electrostatics enabled?:   {:12}\n'.format('true' if opts["ai_elec"] else 'false')
#    text += '  Polarization enabled?:     {:12}\n'.format('true' if opts["ai_pol"] else 'false')
#    text += '  Dispersion enabled?:       {:12}\n'.format('undefined')
#    text += '  Exchange enabled?:         {:12}\n'.format('undefined')
#    text += '  Charge-Transfer enabled?:  {:12}\n'.format('undefined')
#
#    return text


def _pywrapped_get_energy(efpobj, label='libefp'):
    """Gets the energy components from *efpobj* computation as a dictionary.

    Parameters
    ----------
    label : str, optional
        Returned dictionary keys are identical to libefp efp_energy struct
        names plus elec, pol, disp, xr, & total components unless custom renaming
        requested via `label`.

    Returns
    -------
    dict
        Individual terms, summed components, and total energies.

    """
    ene = core.efp_energy()
    res = efpobj.raw_get_energy(ene)
    _pywrapped_result_to_error(res)

    energies = {
        'electrostatic':               ene.electrostatic,
        'charge_penetration':          ene.charge_penetration,
        'electrostatic_point_charges': ene.electrostatic_point_charges,
        'polarization':                ene.polarization,
        'dispersion':                  ene.dispersion,
        'exchange_repulsion':          ene.exchange_repulsion,
        'total':                       ene.total,
        'elec':                        ene.electrostatic +
                                       ene.charge_penetration +
                                       ene.electrostatic_point_charges,
        'xr':                          ene.exchange_repulsion,
        'pol':                         ene.polarization,
        'disp':                        ene.dispersion,
    }

    return _rekey(energies, label=label)


def _pywrapped_get_frag_count(efpobj):
    """Gets the number of fragments in *efpobj* computation.

    Returns
    -------
    int
        Number of fragments in calculation.

    """
    (res, nfrag) = efpobj.cwrapped_get_frag_count()
    _pywrapped_result_to_error(res)

    return nfrag


def _pywrapped_get_multipole_count(efpobj):
    """Gets the number of multipoles in *efpobj* computation.

    Returns
    -------
    int
        Total number of multipoles from electrostatics calculation.

    """
    (res, nmult) = efpobj.cwrapped_get_multipole_count()
    _pywrapped_result_to_error(res)

    return nmult


def _pywrapped_get_multipole_coordinates(efpobj, quiet=False):
    """Gets the coordinates of *efpobj* electrostatics multipoles.

    Parameters
    ----------
    quiet : bool,optional
        Print out the multipole coordinates.

    Returns
    -------
    list
        3 x `n_mult` (flat) array of multipole locations.

    Examples
    --------
    >>> # Use with NumPy
    >>> n_mp = efpobj.get_multipole_count()
    >>> xyz_mp = np.asarray(efpobj.get_multipole_coordinates()).reshape(n_mp, 3)

    """
    nmult = efpobj.get_multipole_count()
    (res, xyz) = efpobj.cwrapped_get_multipole_coordinates(nmult)
    _pywrapped_result_to_error(res)

    if not quiet:
        xyz3 = list(map(list, zip(*[iter(xyz)] * 3)))

        text = '\n  ==>  EFP Multipole Coordinates  <==\n\n'
        for mu in range(nmult):
            text += '{:6d}   {:14.8f} {:14.8f} {:14.8f}\n'.format(
                mu, *xyz3[mu])
        print(text)

    return xyz


def _pywrapped_get_multipole_values(efpobj, quiet=False):
    """Gets the computed per-point multipoles of *efpobj*.

    Parameters
    ----------
    quiet : bool, optional
        Print out the multipole array.

    Returns
    -------
    list
        20 x `n_mult` (flat) array of per-point multipole values including
        charges + dipoles + quadrupoles + octupoles.
        dipoles stored as     x,y,z
        quadrupoles stored as xx,yy,zz,xy,xz,yz
        octupoles stored as   xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz

    Examples
    --------
    >>> # Use with NumPy
    >>> n_mp = efpobj.get_multipole_count()
    >>> val_mp = np.asarray(efpobj.get_multipole_values()).reshape(n_mp, 20)

    """
    nmult = efpobj.get_multipole_count()
    (res, mult) = efpobj.cwrapped_get_multipole_values(nmult)
    _pywrapped_result_to_error(res)

    if not quiet:
        mult20 = list(map(list, zip(*[iter(mult)] * 20)))

        text = '\n  ==>  EFP Multipoles: Charge & Dipole  <==\n\n'
        for mu in range(nmult):
            text += '{:6d}   {:14.8f}   {:14.8f} {:14.8f} {:14.8f}\n'.format(
                mu, *mult20[mu][:4])
        text += '\n  ==>  EFP Multipoles: Quadrupole  <==\n\n'
        for mu in range(nmult):
            text += '{:6d}   {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n'.format(
                mu, *mult20[mu][4:10])
        text += '\n  ==>  EFP Multipoles: Octupole  <==\n\n'
        for mu in range(nmult):
            text += '{:6d}   {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n'.format(
                mu, *mult20[mu][10:])
        print(text)

    return mult


def _pywrapped_get_induced_dipole_count(efpobj):
    """Gets the number of polarization induced dipoles in *efpobj* computation.

    Returns
    -------
    int
        Total number of polarization induced dipoles.

    """
    (res, ndip) = efpobj.cwrapped_get_induced_dipole_count()
    _pywrapped_result_to_error(res)

    return ndip


def _pywrapped_get_induced_dipole_coordinates(efpobj, quiet=False):
    """Gets the coordinates of *efpobj* induced dipoles.

    Parameters
    ----------
    quiet : bool, optional
        Print out the induced dipole coordinates.

    Returns
    -------
    list
        3 x n_dip (flat) array of induced dipole locations.

    """
    ndip = efpobj.get_induced_dipole_count()
    (res, xyz) = efpobj.cwrapped_get_induced_dipole_coordinates(ndip)
    _pywrapped_result_to_error(res)

    if not quiet:
        xyz3 = list(map(list, zip(*[iter(xyz)] * 3)))

        text = '\n  ==>  EFP Induced Dipole Coordinates  <==\n\n'
        for mu in range(ndip):
            text += '{:6d}   {:14.8f} {:14.8f} {:14.8f}\n'.format(
                mu, *xyz3[mu])
        print(text)

    return xyz


def _pywrapped_get_induced_dipole_values(efpobj, quiet=False):
    """Gets the values of polarization induced dipoles of *efpobj*.

    Parameters
    ----------
    quiet : bool, optional
        Print out the induced dipole array.

    Returns
    -------
    list
        3 x n_dip (flat) array of polarization induced dipole values.

    Examples
    --------
    Use with NumPy
    >>> n_dip = efpobj.get_induced_dipole_count()
    >>> val_dip = np.asarray(efpobj.get_induced_dipole_values()).reshape(n_dip, 3)

    """
    ndip = efpobj.get_induced_dipole_count()
    (res, vals) = efpobj.cwrapped_get_induced_dipole_values(ndip)
    _pywrapped_result_to_error(res)

    if not quiet:
        vals3 = list(map(list, zip(*[iter(vals)] * 3)))

        text = '\n  ==>  EFP Induced Dipoles  <==\n\n'
        for mu in range(ndip):
            text += '{:6d}   {:14.8f} {:14.8f} {:14.8f}\n'.format(
                mu, *vals3[mu])
        print(text)

    return vals


def _pywrapped_get_induced_dipole_conj_values(efpobj, quiet=False):
    """Gets the values of polarization conjugated induced dipoles of *efpobj*.

    Parameters
    ----------
    quiet : bool, optional
        Print out the induced dipole array.

    Returns
    -------
    list
        3 x n_dip (flat) array of conjugate induced dipole values.

    """
    ndip = efpobj.get_induced_dipole_count()
    (res, vals) = efpobj.cwrapped_get_induced_dipole_conj_values(ndip)
    _pywrapped_result_to_error(res)

    if not quiet:
        vals3 = list(map(list, zip(*[iter(vals)] * 3)))

        text = '\n  ==>  EFP Conj. Induced Dipoles  <==\n\n'
        for mu in range(ndip):
            text += '{:6d}   {:14.8f} {:14.8f} {:14.8f}\n'.format(
                mu, *vals3[mu])
        print(text)

    return vals


def _pywrapped_get_wavefunction_dependent_energy(efpobj):
    """Updates wavefunction-dependent energy terms for SCF.

    Returns
    -------
    float
        Wavefunction-dependent EFP energy.

    """
    (res, wde) = efpobj.cwrapped_get_wavefunction_dependent_energy()
    _pywrapped_result_to_error(res)

    return wde


def _pywrapped_get_gradient(efpobj, quiet=False):
    """Gets the computed per-fragment EFP energy gradient of *efpobj*.

    Parameters
    ----------
    quiet : bool, optional
        Print out the gradient array.

    Returns
    -------
    list
        6 x n_frag array of per-fragment negative force and torque.

    """
    nfrag = efpobj.get_frag_count()
    (res, grad) = efpobj.cwrapped_get_gradient(nfrag)
    _pywrapped_result_to_error(res)

    if not quiet:
        grad6 = list(map(list, zip(*[iter(grad)] * 6)))

        text = '\n  ==> EFP Gradient <==\n\n'
        for fr in range(nfrag):
            text += '{:14.8f} {:14.8f} {:14.8f}   {:14.8f} {:14.8f} {:14.8f}\n'.format(
                *grad6[fr])
        print(text)

    return grad


def energy_summary(efpobj, label='libefp', scfefp=None):
    """Forms summary of EFP and SCFEFP energy components from `efpobj`.

    Parameters
    ----------
    label : str, optional
        Text labels use libefp terms names unless custom renaming
        requested via `label`.
    scfefp : float, optional
        Total SCF energy (Hartrees) including EFP wavefunction dependent
        and wavefunction independent terms. Used for add'l printing.

    Returns
    -------
    str
        Summary suitable for printing indicating energy components
        (electrostatics, exchange, induction, dispersion, total), whether
        each are enabled in options, and breakdown into pure-EFP and
        QM-EFP, where available. If scfefp, includes a section on SCF
        iterated and SCF total.

    """
    opt = efpobj.get_opts()
    ene = efpobj.get_energy()

    def _enabled(tf, t='*', f=''):
        if tf:
            return t
        else:
            return f

    elec = 'Electrostatics'
    disp = 'Dispersion'
    if label == 'libefp':
        xr = 'Exchange-Repulsion'
        indc = 'Polarization'
    elif label == 'psi':
        xr = 'Exchange'
        indc = 'Induction'

    text = '\n'
    text += '\n    EFP Results\n'
    text +=   '  ------------------------------------------------------------\n'
    text +=   '    {:<30}{:20.12f} [Eh] {}\n'.format(elec, ene['electrostatic'] +
                                                           ene['charge_penetration'] +
                                                           ene['electrostatic_point_charges'],
                                                     _enabled(opt['elec'] or opt['ai_elec']))
    text += '      {:<28}{:20.12f} [Eh] {}\n'.format('EFP/EFP', ene['electrostatic'] +
                                                                ene['charge_penetration'],
                                                     _enabled(opt['elec']))
    text += '      {:<28}{:20.12f} [Eh] {}\n'.format('QM-Nuc/EFP',
                                                     ene['electrostatic_point_charges'],
                                                     _enabled(opt['ai_elec']))
    text += '\n    {:<30}{:20.12f} [Eh] {}\n'.format(xr, ene['exchange_repulsion'],
                                                     _enabled(opt['xr']))
    text += '      {:<28}{:20.12f} [Eh] {}\n'.format('EFP/EFP', ene['exchange_repulsion'],
                                                     _enabled(opt['xr']))
    text += '      {:<28}{:20.12f} [Eh] {}\n'.format('QM/EFP', 0.0, '')
    text += '\n    {:<30}{:20.12f} [Eh] {}\n'.format(indc, ene['polarization'],
                                                     _enabled(opt['pol'] or opt['ai_pol']))
    text += '      {:<28}{:20.12f} [Eh] {}\n'.format(_enabled(opt['ai_pol'], t='QM/EFP', f='EFP/EFP'),
                                                     ene['polarization'],
                                                     _enabled(opt['pol'] or opt['ai_pol']))
    text += '\n    {:<30}{:20.12f} [Eh] {}\n'.format(disp, ene['dispersion'], _enabled(opt['disp']))
    text += '      {:<28}{:20.12f} [Eh] {}\n'.format('EFP/EFP', ene['dispersion'],
                                                     _enabled(opt['disp']))
    text += '      {:<28}{:20.12f} [Eh] {}\n'.format('QM/EFP', 0.0, '')
    text += '\n    {:<30}{:20.12f} [Eh]\n'.format('Total EFP', ene['total'])
    if scfefp is not None:
        wie = ene['total'] - ene['pol']
        text +=   '    EFP excluding EFP {:<12}{:20.12f} [Eh]\n'.format(indc, wie)
        text +=   '    SCF including EFP {:<12}{:20.12f} [Eh]\n'.format(indc, scfefp - wie)
        text +=   '    Total SCF including Total EFP {:20.12f} [Eh]\n'.format(scfefp)

    return text


#def nuclear_repulsion_energy(efpobj):
#    """Computes nuclear repulsion energy."""
#
#    pyat = efpobj.get_atoms()
#    nre = 0.0
#    for iat1, at1 in enumerate(pyat['full_atoms']):
#        for iat2, at2 in enumerate(pyat['full_atoms']):
#            if iat2 < iat1:
#                ZZ = at1['Z'] * at2['Z']
#                dx = at1['x'] - at2['x']
#                dy = at1['y'] - at2['y']
#                dz = at1['z'] - at2['z']
#                dist = math.sqrt(dx * dx + dy * dy + dz * dz)
#                nre += ZZ / dist
#
#    return nre


#def _frag_idx_validation(efpobj, ifr):
#    nfr = efpobj.get_frag_count()
#    if (ifr < 0) or (ifr >= nfr):
#        raise PyEFPSyntaxError('Invalid fragment index for 0-indexed {}-fragment EFP: {}'.format(nfr, ifr))

def _pywrapped_set_frag_coordinates(efpobj, ifr, ctype, coord):
    """Set fragment orientation on `efpobj` from hint.

    Parameters
    ----------
    ifr : int
        Index of fragment (0-indexed).
    ctype : core.efp_coord_type or str
        Type of coodinates hint among `xyzabc`, `points`, & `rotmat`.
    coord : list of floats
        6-, 9-, or 12-element hint of coordinates.

    Returns
    -------
    None

    """
    if isinstance(ctype, basestring):
        try:
            ctype = {'xyzabc': core.EFP_COORD_TYPE_XYZABC,
                     'points': core.EFP_COORD_TYPE_POINTS,
                     'rotmat': core.EFP_COORD_TYPE_ROTMAT,
                    }[ctype.lower()]
        except KeyError:
            _pywrapped_result_to_error(core.efp_result.EFP_RESULT_SYNTAX_ERROR,
                                       'invalid value for [xyzabc/points/rotmat] {}: {}'.format('ctype', ctype))

    res = efpobj.cwrapped_set_frag_coordinates(ifr, ctype, coord)
    _pywrapped_result_to_error(res)


def _pywrapped_set_point_charges(efpobj, ptc, coord):

    if (len(ptc) * 3) != len(coord):
        print('Protest~')

    res = efpobj.cwrapped_set_point_charges(len(ptc), ptc, coord)
    _pywrapped_result_to_error(res)


def _pywrapped_get_frag_name(efpobj, ifr=None):
    """Gets system name on fragment(s) of `efpobj`.

    Parameters
    ----------
    ifr : int, optional
        Index of fragment (0-indexed) if not all.

    Returns
    -------
    str, list of str
        If `ifr`, name of fragment `ifr`. Otherwise, names of all
        fragments in list.

    """
    nfr = efpobj.get_frag_count()

    if ifr is None:
        frags = []
        for fr in range(nfr):
            (res, fname) = efpobj.cwrapped_get_frag_name(fr)
            _pywrapped_result_to_error(res)
            frags.append(fname)

        return frags

    else:
        if ifr in range(nfr):
            (res, fname) = efpobj.cwrapped_get_frag_name(ifr)
            _pywrapped_result_to_error(res)

            return fname
        else:
            raise PyEFPSyntaxError('Invalid fragment index for 0-indexed {}-fragment EFP: {}'.format(nfr, ifr))


def _pywrapped_get_frag_charge(efpobj, ifr=None, zero=1e-8):
    """Gets total charge on fragment(s) of `efpobj`.

    Parameters
    ----------
    ifr : int, optional
        Index of fragment (0-indexed) if not all.
    zero : float, optional
        Absolute value under which to zero charge.

    Returns
    -------
    str, list of str
        If `ifr`, charge of fragment `ifr`. Otherwise, charges of all
        fragments in list.

    """
    nfr = efpobj.get_frag_count()

    if ifr is None:
        frags = []
        for fr in range(nfr):
            (res, chg) = efpobj.cwrapped_get_frag_charge(fr)
            _pywrapped_result_to_error(res)

            if math.fabs(chg) < zero:
                frags.append(0.0)
            else:
                frags.append(chg)
        return frags

    else:
        if ifr in range(nfr):
            (res, chg) = efpobj.cwrapped_get_frag_charge(ifr)
            _pywrapped_result_to_error(res)

            if math.fabs(chg) < zero:
                return 0.0
            else:
                return chg
        else:
            raise PyEFPSyntaxError('Invalid fragment index for 0-indexed {}-fragment EFP: {}'.format(nfr, ifr))


def _pywrapped_get_frag_multiplicity(efpobj, ifr=None):
    """Gets spin multiplicity on fragment(s) of `efpobj`.

    Parameters
    ----------
    ifr : int, optional
        Index of fragment (0-indexed) if not all.

    Returns
    -------
    str, list of str
        If `ifr`, multiplicity of fragment `ifr`. Otherwise, multiplicity
        of all fragments in list.

    """
    nfr = efpobj.get_frag_count()

    if ifr is None:
        frags = []
        for fr in range(nfr):
            (res, mult) = efpobj.cwrapped_get_frag_multiplicity(fr)
            _pywrapped_result_to_error(res)
            frags.append(mult)
        return frags

    else:
        if ifr in range(nfr):
            (res, mult) = efpobj.cwrapped_get_frag_multiplicity(ifr)
            _pywrapped_result_to_error(res)

            return mult
        else:
            raise PyEFPSyntaxError('Invalid fragment index for 0-indexed {}-fragment EFP: {}'.format(nfr, ifr))


def _pywrapped_get_frag_atom_count(efpobj, ifr=None):
    """Gets atom count on fragment(s) of `efpobj`.

    Parameters
    ----------
    ifr : int, optional
        Index of fragment (0-indexed) if not all.

    Returns
    -------
    str, list of str
        If `ifr`, atom count of fragment `ifr`. Otherwise, atom counts
        of all fragments in list.

    """
    nfr = efpobj.get_frag_count()

    if ifr is None:
        frags = []
        for fr in range(nfr):
            (res, natom) = efpobj.cwrapped_get_frag_atom_count(fr)
            _pywrapped_result_to_error(res)
            frags.append(natom)
        return frags

    else:
        if ifr in range(nfr):
            (res, natom) = efpobj.cwrapped_get_frag_atom_count(ifr)
            _pywrapped_result_to_error(res)

            return natom
        else:
            raise PyEFPSyntaxError('Invalid fragment index for 0-indexed {}-fragment EFP: {}'.format(nfr, ifr))


def _pywrapped_get_frag_atoms(efpobj, ifr):
    """Gets geometry information for atoms modeled by fragment in `efpobj`.

    Parameters
    ----------
    ifr : int
        Index of fragment (0-indexed).

    Returns
    -------
    list of dict
        Each atom in fragment `ifr` has position, charge, and element
        fields below in a dictionary at list index `ifr`
        Z : float               nuclear charge.
        label : str             atom label from EFP file, e.g., A02H2.
        x : float               X coordinate of atom position.
        y : float               Y coordinate of atom position.
        z : float               Z coordinate of atom position.
        mass : float            atom mass [amu]
        symbol : str            atomic symbol extracted from label.

    """
    nfr = efpobj.get_frag_count()
    nat = efpobj.get_frag_atom_count(ifr)

    if ifr in range(nfr):
        (res, atoms) = efpobj.cwrapped_get_frag_atoms(ifr, nat)
        _pywrapped_result_to_error(res)

        for at in atoms:
            mobj = re.match(r'\AA\d*(?P<symbol>[A-Z]{1,3})\d*\Z', at['label'])
            if mobj:
                at['symbol'] = mobj.group('symbol').capitalize()
#            at['xyz'] = [at['x'], at['y'], at['z']]

        return atoms
    else:
        raise PyEFPSyntaxError('Invalid fragment index for 0-indexed {}-fragment EFP: {}'.format(nfr, ifr))



def py_get_atoms(efpobj):
    #enum efp_result res;
    #size_t frag_natom, natom=0;
    #double frag_chg;
    #int frag_mult;
    #py::list fr, frt, frcg, frmp, full_atoms;


    natom = 0
    frag_count = efpobj.get_frag_count()
    frag_natom = efpobj.get_frag_atom_count()
    fr = []
    full_atoms = []
    for ifr in range(frag_count):
        frat = frag_natom[ifr]

        fr.append([natom, natom + frat])
        natom += frat

        pyat = efpobj.get_frag_atoms(ifr)
        full_atoms.append(pyat)

    mol_info = {}
    mol_info["units"] = "Bohr"
    mol_info["input_units_to_au"] = 1.0
    mol_info["fix_com"] = True
    mol_info["fix_orientation"] = True
    mol_info["fix_symmetry"] = "c1"

    mol_info['fragments'] = fr
    mol_info['fragment_types'] = ['Real'] * frag_count
    mol_info['fragment_charges'] = efpobj.get_frag_charge()
    mol_info['fragment_multiplicities'] = efpobj.get_frag_multiplicity()
    mol_info['full_atoms'] = full_atoms

    return mol_info


def to_viz_dict(efpobj):

    pyat = efpobj.get_atoms()
    for at in pyat['full_atoms']:
        at['ghosted'] = False
        at['at_type'] = 'efpxyz'
        mobj = re.match(r'\AA\d*(?P<symbol>[A-Z]{1,3})\d*\Z', at['label'])
        if mobj:
            at['symbol'] = mobj.group('symbol').capitalize()
        at['charge'] = at['Z']
        #pyat['molecule']['fragment_charges'].append(efpobj.get_frag_charges(fr)
        #pyat['molecule']['fragment_multiplicities'].append(efpobj.get_frag_multiplicity(fr)

    return pyat


def to_dict(efpobj):
    pysys = {}
    pysys['full_fragments'] = []

    for fr in range(efpobj.get_frag_count()):
        (res, xyzabc) = efpobj.cwrapped_get_frag_xyzabc(fr)
        _pywrapped_result_to_error(res)

        pysys['full_fragments'].append({
            'coordinates_hint': xyzabc,
            'efp_type': 'xyzabc',
            'fragment_file': efpobj.get_frag_name(fr).lower(),
                   })

    pysys['molecule'] = {
        'fix_com': True,
        'fix_orientation': True,
        'fix_symmetry': 'c1',
        'fragment_charges': [],
        'fragment_multiplicities': [],
        'fragment_types': [],
        'fragments': [],
        'full_atoms': [],
        'input_units_to_au': 1.8897261328856432,
        'name': 'default',
        'units': 'Bohr'}

    return pysys


yuio = {'libefp': {'full_fragments':
 [{'coordinates_hint': [-0.5753870821672306,
                       -4.23695594520049,
                       -0.5552607051670226,
                       -0.642499,
                       1.534222,
                       -0.568147],
  'efp_type': 'xyzabc',
  'fragment_file': 'C6H6'},
 {'coordinates_hint': [-1.1352612324342508,
                       2.578405376972965,
                       1.4862284641766452,
                       3.137879,
                       1.557344,
                       -2.56855],
  'efp_type': 'xyzabc',
  'fragment_file': 'C6H6'}]},
 'molecule': {
 'fix_com': True,
 'fix_orientation': True,
 'fix_symmetry': 'c1',
 'fragment_charges': [],
 'fragment_multiplicities': [],
 'fragment_types': [],
 'fragments': [],
 'full_atoms': [],
 'input_units_to_au': 1.8897261328856432,
 'name': 'default',
 'units': 'Angstrom'}}

# only wrapped to throw Py exceptions
core.efp.prepare = _pywrapped_efp_prepare
core.efp.compute = _pywrapped_efp_compute

core.efp.add_potential = _pywrapped_add_potential
core.efp.add_fragment = _pywrapped_add_fragment
core.efp.get_opts = _pywrapped_get_opts
core.efp.set_opts = _pywrapped_set_opts
core.efp.get_frag_count = _pywrapped_get_frag_count
core.efp.get_energy = _pywrapped_get_energy
core.efp.get_gradient = _pywrapped_get_gradient
core.efp.energy_summary = energy_summary
#core.efp.nuclear_repulsion_energy = nuclear_repulsion_energy
core.efp.to_viz_dict = to_viz_dict
core.efp.to_dict = to_dict
core.efp.get_frag_name = _pywrapped_get_frag_name
core.efp.get_frag_charge = _pywrapped_get_frag_charge
core.efp.get_frag_multiplicity = _pywrapped_get_frag_multiplicity
core.efp.set_frag_coordinates = _pywrapped_set_frag_coordinates
core.efp.set_point_charges = _pywrapped_set_point_charges
core.efp.get_multipole_count = _pywrapped_get_multipole_count
core.efp.get_multipole_coordinates = _pywrapped_get_multipole_coordinates
core.efp.get_multipole_values = _pywrapped_get_multipole_values
core.efp.get_induced_dipole_count = _pywrapped_get_induced_dipole_count
core.efp.get_induced_dipole_coordinates = _pywrapped_get_induced_dipole_coordinates
core.efp.get_induced_dipole_values = _pywrapped_get_induced_dipole_values
core.efp.get_induced_dipole_conj_values = _pywrapped_get_induced_dipole_conj_values
core.efp.get_frag_atom_count = _pywrapped_get_frag_atom_count
core.efp.get_wavefunction_dependent_energy = _pywrapped_get_wavefunction_dependent_energy

core.efp.get_frag_atoms = _pywrapped_get_frag_atoms
core.efp.py_get_atoms = py_get_atoms


def from_dict(efp_init):

    sys = core.efp()

    for ifr, fr in enumerate(efp_init['full_fragments']):
        sys.add_potential(fr['fragment_file'], duplicates_ok=True)
        sys.add_fragment(fr['fragment_file'])
        sys.set_frag_coordinates(ifr, fr['efp_type'], fr['coordinates_hint'])

    sys.prepare()
    return sys
