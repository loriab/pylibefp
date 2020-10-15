import sys
import pytest
import pprint
import pylibefp
from systems import system_1

from qcelemental.testing import compare_recursive


def test_opts_libefp():
    asdf = system_1()

    ref1 = {
        'ai_elec': False,
        'elec_damp': 'screen',
        'ai_disp': False,
        'chtr': False,
        'swf_cutoff': 0.0,
        'enable_cutoff': False,
        'disp': False,
        'ai_pol': False,
        'pol': False,
        'xr': False,
        'pol_driver': 'iterative',
        'ai_xr': False,
        'elec': False,
        'pol_damp': 'tt',
        'disp_damp': 'overlap',
        'enable_pbc': False,
        'ai_chtr': False
    }
    ans1 = asdf.set_opts({})
    assert compare_recursive(ref1, ans1, sys._getframe().f_code.co_name + ': blank', atol=1.e-6)

    ref2 = {
        'ai_elec': True,
        'elec_damp': 'off',
        'ai_disp': False,
        'chtr': False,
        'swf_cutoff': 1.0,
        'enable_cutoff': False,
        'disp': False,
        'ai_pol': False,
        'pol': False,
        'xr': False,
        'pol_driver': 'iterative',
        'ai_xr': False,
        'elec': True,
        'pol_damp': 'tt',
        'disp_damp': 'overlap',
        'enable_pbc': False,
        'ai_chtr': False
    }
    ans2 = asdf.set_opts({
        'elec_damp': 'OFF',
        'swf_cutoff': 1.0,
        'elec': True,
        'ai_elec': True,
        'enable_cutoff': False
    })
    assert compare_recursive(ref2, ans2, sys._getframe().f_code.co_name + ': setting', atol=1.e-6)

    ref3 = {
        'ai_elec': True,
        'elec_damp': 'off',
        'ai_disp': False,
        'chtr': False,
        'swf_cutoff': 2.0,
        'enable_cutoff': False,
        'disp': False,
        'ai_pol': False,
        'pol': False,
        'xr': False,
        'pol_driver': 'iterative',
        'ai_xr': False,
        'elec': False,
        'pol_damp': 'tt',
        'disp_damp': 'tt',
        'enable_pbc': False,
        'ai_chtr': False
    }
    ans3 = asdf.set_opts({'swf_cutoff': 2, 'elec': False, 'ai_elec': True, 'disp_damp': 'TT'}, append='append')
    assert compare_recursive(ref3, ans3, sys._getframe().f_code.co_name + ': append setting', atol=1.e-6)

    ref4 = {
        'ai_elec': False,
        'elec_damp': 'off',
        'ai_disp': False,
        'chtr': False,
        'swf_cutoff': 0.0,
        'enable_cutoff': False,
        'disp': False,
        'ai_pol': False,
        'pol': False,
        'xr': False,
        'pol_driver': 'iterative',
        'ai_xr': False,
        'elec': True,
        'pol_damp': 'tt',
        'disp_damp': 'overlap',
        'enable_pbc': False,
        'ai_chtr': False
    }
    ans4 = asdf.set_opts({
        'elec_damp': 'OFF',
        'swf_cutoff': 0.0,
        'elec': True,
        'enable_cutoff': False
    },
                         append='libefp')
    assert compare_recursive(ref4, ans4, sys._getframe().f_code.co_name + ': reset setting', atol=1.e-6)


def test_opts_fail_1():
    asdf = system_1()

    asdf.set_opts({'nonsense_key': 'harmless'})
    with pytest.raises(pylibefp.EFPSyntaxError):
        asdf.set_opts({'elec_damp': 'nonsense'})


def test_opts_fail_2():
    asdf = system_1()

    with pytest.raises(pylibefp.EFPSyntaxError):
        asdf.set_opts({'elec': 'yEs'})


def test_opts_psi():
    asdf = system_1()

    ref = {
        'qm_elst': False,
        'elst_damping': 'screen',
        'qm_disp': False,
        'chtr': False,
        'swf_cutoff': 0.0,
        'enable_cutoff': False,
        'disp': False,
        'qm_ind': False,
        'ind': False,
        'exch': False,
        'ind_driver': 'iterative',
        'qm_exch': False,
        'elst': False,
        'ind_damping': 'tt',
        'disp_damping': 'overlap',
        'enable_pbc': False,
        'qm_chtr': False
    }
    ans = asdf.set_opts({}, label='psi')
    assert compare_recursive(ref, ans, sys._getframe().f_code.co_name + ': psi blank', atol=1.e-6)

    ref = {
        'qm_elst': True,
        'elst_damping': 'off',
        'qm_disp': False,
        'chtr': False,
        'swf_cutoff': 1.0,
        'enable_cutoff': False,
        'disp': False,
        'qm_ind': False,
        'ind': False,
        'exch': False,
        'ind_driver': 'iterative',
        'qm_exch': False,
        'elst': True,
        'ind_damping': 'tt',
        'disp_damping': 'overlap',
        'enable_pbc': False,
        'qm_chtr': False
    }
    ans = asdf.set_opts(
        {
            'elst_damping': 'OFF',
            'swf_cutoff': 1.0,
            'elst': True,
            'qm_elst': True,
            'enable_cutoff': False
        }, label='psi')
    assert compare_recursive(ref, ans, sys._getframe().f_code.co_name + ': psi setting', atol=1.e-6)

    ref = {
        'qm_elst': True,
        'elst_damping': 'off',
        'qm_disp': False,
        'chtr': False,
        'swf_cutoff': 2.0,
        'enable_cutoff': False,
        'disp': False,
        'qm_ind': False,
        'ind': False,
        'exch': False,
        'ind_driver': 'iterative',
        'qm_exch': False,
        'elst': False,
        'ind_damping': 'tt',
        'disp_damping': 'tt',
        'enable_pbc': False,
        'qm_chtr': False
    }
    ans = asdf.set_opts({
        'swf_cutoff': 2,
        'elst': False,
        'qm_elst': True,
        'disp_damping': 'TT'
    },
                        append='append',
                        label='psi')
    assert compare_recursive(ref, ans, sys._getframe().f_code.co_name + ': psi append setting', atol=1.e-6)

    ref = {
        'qm_elst': False,
        'elst_damping': 'off',
        'qm_disp': False,
        'chtr': False,
        'swf_cutoff': 0.0,
        'enable_cutoff': False,
        'disp': False,
        'qm_ind': False,
        'ind': False,
        'exch': False,
        'ind_driver': 'iterative',
        'qm_exch': False,
        'elst': True,
        'ind_damping': 'tt',
        'disp_damping': 'overlap',
        'enable_pbc': False,
        'qm_chtr': False
    }
    ans = asdf.set_opts({
        'elst_damping': 'OFF',
        'swf_cutoff': 0.0,
        'elst': True,
        'enable_cutoff': False
    },
                        append='libefp',
                        label='psi')
    assert compare_recursive(ref, ans, sys._getframe().f_code.co_name + ': psi reset setting', atol=1.e-6)


def test_opts_psi_dflt():
    asdf = system_1()

    ref = {
        'qm_elst': True,
        'elst_damping': 'screen',
        'qm_disp': False,
        'chtr': False,
        'swf_cutoff': 0.0,
        'enable_cutoff': False,
        'disp': True,
        'qm_ind': True,
        'ind': True,
        'exch': True,
        'ind_driver': 'iterative',
        'qm_exch': False,
        'elst': True,
        'ind_damping': 'tt',
        'disp_damping': 'overlap',
        'enable_pbc': False,
        'qm_chtr': False
    }
    ans = asdf.set_opts({}, label='psi', append='psi')
    assert compare_recursive(ref, ans, sys._getframe().f_code.co_name + ': psi default blank', atol=1.e-6)

    ref = {
        'qm_elst': False,
        'elst_damping': 'off',
        'qm_disp': False,
        'chtr': False,
        'swf_cutoff': 1.0,
        'enable_cutoff': False,
        'disp': True,
        'qm_ind': True,
        'ind': True,
        'exch': True,
        'ind_driver': 'iterative',
        'qm_exch': False,
        'elst': True,
        'ind_damping': 'tt',
        'disp_damping': 'overlap',
        'enable_pbc': False,
        'qm_chtr': False
    }
    ans = asdf.set_opts(
        {
            'elst_damping': 'OFF',
            'swf_cutoff': 1.0,
            'elst': True,
            'qm_elst': False,
            'enable_cutoff': False
        },
        label='psi',
        append='append')
    assert compare_recursive(ref, ans, sys._getframe().f_code.co_name + ': psi default append setting', atol=1.e-6)

    ref = {
        'qm_elst': True,
        'elst_damping': 'overlap',
        'qm_disp': False,
        'chtr': False,
        'swf_cutoff': 2.0,
        'enable_cutoff': True,
        'disp': True,
        'qm_ind': True,
        'ind': True,
        'exch': True,
        'ind_driver': 'iterative',
        'qm_exch': False,
        'elst': False,
        'ind_damping': 'tt',
        'disp_damping': 'overlap',
        'enable_pbc': False,
        'qm_chtr': False
    }
    ans = asdf.set_opts({
        'elst_damping': 'OVERlap',
        'swf_cutoff': 2.0,
        'elst': False,
        'enable_cutoff': True
    },
                        append='psi',
                        label='psi')
    assert compare_recursive(ref, ans, sys._getframe().f_code.co_name + ': psi default reset setting', atol=1.e-6)
