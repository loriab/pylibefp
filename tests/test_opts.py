import sys
import pytest
import pprint
import pylibefp
from utils import *
from systems import *


def test_opts_libefp():
    asdf = system_1()

    ref = {'ai_elec': False, 'elec_damp': 'screen', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 0.0, 'enable_cutoff': 0, 'disp': False, 'ai_pol': False, 'pol': False, 'xr': False, 'pol_driver': 'iterative', 'ai_xr': False, 'elec': False, 'pol_damp': 'tt', 'disp_damp': 'overlap', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({})
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': blank'))

    ref = {'ai_elec': True, 'elec_damp': 'off', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 1.0, 'enable_cutoff': 0, 'disp': False, 'ai_pol': False, 'pol': False, 'xr': False, 'pol_driver': 'iterative', 'ai_xr': False, 'elec': True, 'pol_damp': 'tt', 'disp_damp': 'overlap', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({'elec_damp': 'OFF', 'swf_cutoff': 1.0, 'elec': True, 'ai_elec': True, 'enable_cutoff': 0})
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': setting'))

    ref = {'ai_elec': True, 'elec_damp': 'off', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 2.0, 'enable_cutoff': 0, 'disp': False, 'ai_pol': False, 'pol': False, 'xr': False, 'pol_driver': 'iterative', 'ai_xr': False, 'elec': False, 'pol_damp': 'tt', 'disp_damp': 'tt', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({'swf_cutoff': 2, 'elec': False, 'ai_elec': True, 'disp_damp': 'TT'}, append='append')
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': append setting'))

    ref = {'ai_elec': False, 'elec_damp': 'off', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 0.0, 'enable_cutoff': 0, 'disp': False, 'ai_pol': False, 'pol': False, 'xr': False, 'pol_driver': 'iterative', 'ai_xr': False, 'elec': True, 'pol_damp': 'tt', 'disp_damp': 'overlap', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({'elec_damp': 'OFF', 'swf_cutoff': 0.0, 'elec': True, 'enable_cutoff': 0}, append='libefp')
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': reset setting'))


def test_opts_fail_1():
    asdf = system_1()

    with pytest.raises(pylibefp.EFPSyntaxError) as e_info:
        ans = asdf.set_opts({'elec_damp': 'nonsense'})


def test_opts_fail_2():
    asdf = system_1()

    with pytest.raises(pylibefp.EFPSyntaxError) as e_info:
        ans = asdf.set_opts({'elec': 1})


def test_opts_psi():
    asdf = system_1()

    ref = {'ai_elst': False, 'elst_damp': 'screen', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 0.0, 'enable_cutoff': 0, 'disp': False, 'ai_ind': False, 'ind': False, 'exch': False, 'ind_driver': 'iterative', 'ai_exch': False, 'elst': False, 'ind_damp': 'tt', 'disp_damp': 'overlap', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({}, label='psi')
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': psi blank'))

    ref = {'ai_elst': True, 'elst_damp': 'off', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 1.0, 'enable_cutoff': 0, 'disp': False, 'ai_ind': False, 'ind': False, 'exch': False, 'ind_driver': 'iterative', 'ai_exch': False, 'elst': True, 'ind_damp': 'tt', 'disp_damp': 'overlap', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({'elst_damp': 'OFF', 'swf_cutoff': 1.0, 'elst': True, 'ai_elst': True, 'enable_cutoff': 0}, label='psi')
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': psi setting'))

    ref = {'ai_elst': True, 'elst_damp': 'off', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 2.0, 'enable_cutoff': 0, 'disp': False, 'ai_ind': False, 'ind': False, 'exch': False, 'ind_driver': 'iterative', 'ai_exch': False, 'elst': False, 'ind_damp': 'tt', 'disp_damp': 'tt', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({'swf_cutoff': 2, 'elst': False, 'ai_elst': True, 'disp_damp': 'TT'}, append='append', label='psi')
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': psi append setting'))

    ref = {'ai_elst': False, 'elst_damp': 'off', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 0.0, 'enable_cutoff': 0, 'disp': False, 'ai_ind': False, 'ind': False, 'exch': False, 'ind_driver': 'iterative', 'ai_exch': False, 'elst': True, 'ind_damp': 'tt', 'disp_damp': 'overlap', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({'elst_damp': 'OFF', 'swf_cutoff': 0.0, 'elst': True, 'enable_cutoff': 0}, append='libefp', label='psi')
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': psi reset setting'))


def test_opts_psi_dflt():
    asdf = system_1()

    ref = {'ai_elst': True, 'elst_damp': 'screen', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 0.0, 'enable_cutoff': 0, 'disp': True, 'ai_ind': True, 'ind': True, 'exch': True, 'ind_driver': 'iterative', 'ai_exch': False, 'elst': True, 'ind_damp': 'tt', 'disp_damp': 'overlap', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({}, label='psi', append='psi')
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': psi default blank'))

    ref = {'ai_elst': False, 'elst_damp': 'off', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 1.0, 'enable_cutoff': 0, 'disp': True, 'ai_ind': True, 'ind': True, 'exch': True, 'ind_driver': 'iterative', 'ai_exch': False, 'elst': True, 'ind_damp': 'tt', 'disp_damp': 'overlap', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({'elst_damp': 'OFF', 'swf_cutoff': 1.0, 'elst': True, 'ai_elst': False, 'enable_cutoff': 0}, label='psi', append='append')
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': psi default append setting'))

    ref = {'ai_elst': True, 'elst_damp': 'overlap', 'ai_disp': False, 'chtr': False, 'swf_cutoff': 2.0, 'enable_cutoff': 0, 'disp': True, 'ai_ind': True, 'ind': True, 'exch': True, 'ind_driver': 'iterative', 'ai_exch': False, 'elst': False, 'ind_damp': 'tt', 'disp_damp': 'overlap', 'enable_pbc': 0, 'ai_chtr': False}
    ans = asdf.set_opts({'elst_damp': 'OVERlap', 'swf_cutoff': 2.0, 'elst': False, 'enable_cutoff': 0}, append='psi', label='psi')
    assert(compare_dicts(ref, ans, 6, sys._getframe().f_code.co_name + ': psi default reset setting'))

