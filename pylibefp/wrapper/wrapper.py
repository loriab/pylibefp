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

    print('GET BEFORE: ', dopts)
    def rekey(rawdict):
        newdict = rawdict.copy()
        for key in rawdict.keys():
            topic = _lbtl[label].get(key, key)
            newdict[topic] = newdict.pop(key)
        return newdict

    dopts = rekey(dopts)
    print('GET AFTER: ', dopts)

    return dopts


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
    elif append == 'append':
        res = efpobj.raw_get_opts(opts)
        _pywrapped_result_to_error(res)
    else:
        raise EFPSyntaxError('Unrecognized opts default set: {}'.format(append))

    print('SET BEFORE:', dopts)
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

    ropts = efpobj.get_opts(label=label)
    print('SET AFTER:', ropts)
    return ropts



def opts_summary(efpobj, labels='libefp'):

    opts = efpobj.get_opts()
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


def _pywrapped_get_energy(efpobj):

    ene = core.efp_energy()
    res = efpobj.raw_get_energy(ene)
    _pywrapped_result_to_error(res)

#    if (do_grad_) {
#        SharedMatrix smgrad(new Matrix("EFP Gradient", nfrag_, 6));
#        double ** psmgrad = smgrad->pointer();
#        if ((res = efp_get_gradient(efp_, psmgrad[0])))
#            throw PsiException("EFP::compute():efp_get_gradient(): " +
#                std::string (efp_result_to_string(res)),__FILE__,__LINE__);
#        smgrad->print_out();
#
#        outfile->Printf("  ==> EFP Gradient <==\n\n");
#
#        for (int i=0; i<nfrag_; i++) {
#            for (int j=0; j<6; j++) {
#                outfile->Printf("%14.6lf", psmgrad[i][j]);
#            }
#            outfile->Printf("\n");
#        }
#        outfile->Printf("\n");
#
#        torque_ = smgrad;
#
#        std::shared_ptr<Wavefunction> wfn = Process::environment.legacy_wavefunction();
#        wfn->set_efp_torque(smgrad);
#    }

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

    return energies


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


def nuclear_repulsion_energy(efpobj):
    """Computes nuclear repulsion energy."""

    pyat = efpobj.get_atoms()
    nre = 0.0
    for iat1, at1 in enumerate(pyat['full_atoms']):
        for iat2, at2 in enumerate(pyat['full_atoms']):
            if iat2 < iat1:
                ZZ = at1['Z'] * at2['Z']
                dx = at1['x'] - at2['x']
                dy = at1['y'] - at2['y']
                dz = at1['z'] - at2['z']
                dist = math.sqrt(dx * dx + dy * dy + dz * dz)
                nre += ZZ / dist

    return nre


def to_dict(efpobj):

    pyat = efpobj.get_atoms()
    for at in pyat['full_atoms']:
        at['ghosted'] = False
        at['at_type'] = 'efpxyz'
        mobj = re.match(r'\AA\d*(?P<symbol>[A-Z]{1,3})\d*\Z', at['label'])
        if mobj:
            at['symbol'] = mobj.group('symbol').capitalize()
        at['charge'] = at['Z']

    return pyat 


# only wrapped to throw Py exceptions
core.efp.prepare = _pywrapped_efp_prepare
core.efp.compute = _pywrapped_efp_compute

core.efp.add_potentials = _pywrapped_add_potential
core.efp.add_fragments = _pywrapped_add_fragment
core.efp.get_opts = _pywrapped_get_opts
core.efp.set_opts = _pywrapped_set_opts
core.efp.get_energy = _pywrapped_get_energy
core.efp.energy_summary = energy_summary
core.efp.nuclear_repulsion_energy = nuclear_repulsion_energy
core.efp.to_dict = to_dict


def from_dict(efp_init):

    sys = core.efp()
    sys.create()

    for ifr, fr in enumerate(efp_init['full_fragments']):
        sys.add_potentials(fr['fragment_file'])
        sys.add_fragments(fr['fragment_file'])
        sys.set_frag_coordinates(ifr, fr['efp_type'], fr['coordinates_hint'])

    sys.prepare()
    return sys


