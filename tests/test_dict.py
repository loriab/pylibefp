import sys
import pytest
import pprint
import pylibefp
from systems import *

from qcelemental.testing import compare, compare_values


def test_dict_1():
    sys1 = system_1()

    dict1 = sys1.to_dict()
    print('DICT1')
    pprint.pprint(dict1)
    sys1p = pylibefp.from_dict(dict1)

    sys1p.set_opts({'elec': True, 'elec_damp': 'screen', 'xr': True, 'pol': True, 'disp': True, 'disp_damp': 'tt'})
    sys1p.compute()
    ene = sys1p.get_energy()

    pprint.pprint(ene)
    print('<<< get_opts():  ', sys1p.get_opts(), '>>>')
    print('<<< get_energy():', ene, '>>>')
    print('<<< get_atoms(): ', sys1p.get_atoms(), '>>>')
    print(sys1p.energy_summary())
    print(sys1p.geometry_summary(units_to_bohr=b2a))
    print(sys1p.geometry_summary(units_to_bohr=1.0))

    assert compare(2, sys1p.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag')
    assert compare_values(0.0001922903, ene['total'], sys._getframe().f_code.co_name + ': ene', atol=1.e-6)


def test_dict_2a():
    sys1 = system_2()
    sys1p = pylibefp.from_dict(sys1.to_dict())

    sys1p.set_opts({'elec': True, 'pol': True, 'disp': True, 'xr': True, 'elec_damp': 'screen', 'disp_damp': 'tt'})
    sys1p.compute()
    ene = sys1p.get_energy()

    assert compare(5, sys1p.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag')
    assert compare_values(0.0007440865, ene['total'], sys._getframe().f_code.co_name, atol=1.e-6)


def test_dict_3a():
    sys1 = system_3()
    sys1p = pylibefp.from_dict(sys1.to_dict())

    sys1p.set_opts({
        'elec': True,
        'pol': True,
        'disp': True,
        'xr': True,
        'elec_damp': 'screen',
        'disp_damp': 'tt',
        'pol_damp': 'tt'
    })
    sys1p.compute()
    ene = sys1p.get_energy()

    assert compare(9, sys1p.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag')
    assert compare_values(0.0061408841, ene['total'], sys._getframe().f_code.co_name, atol=1.e-5)


def test_dict_4a():
    sys1 = system_4()
    sys1p = pylibefp.from_dict(sys1.to_dict())

    sys1p.set_opts({
        'elec': True,
        'pol': True,
        'disp': True,
        'xr': True,
        'elec_damp': 'screen',
        'disp_damp': 'tt',
        'pol_damp': 'tt'
    })
    sys1p.compute()
    ene = sys1p.get_energy()

    assert compare(12, sys1p.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag')
    assert compare_values(-0.0095597483, ene['total'], sys._getframe().f_code.co_name, atol=1.e-5)


def test_dict_5():
    dsys = {'units': 'Angstrom', 'fragment_files': [], 'hint_types': [], 'geom_hints': []}

    sys = pylibefp.from_dict(dsys)

    with pytest.raises(pylibefp.PolNotConverged) as e_info:
        sys.compute()
    assert sys.get_frag_count() == 0


def test_subset_1():
    sys2 = system_2()
    ss15_files = [sys2.get_frag_name(0).lower(), sys2.get_frag_name(4).lower()]
    ss15_xyzabc = [sys2.get_frag_xyzabc(0), sys2.get_frag_xyzabc(4)]

    ss15 = pylibefp.extract_subsets(sys2, [1,5])
    dss15 = ss15.to_dict()

    assert compare(ss15_files, dss15["fragment_files"], "ss_fragfiles")
    assert compare_values(ss15_xyzabc, dss15["geom_hints"], "ss_xyzabc", atol=1.e-6)


@pytest.mark.parametrize("args", [
    {"reals": 0},
    {"reals": [2, 6]},
    {"reals": 2, "ghosts": 3},
])
def test_subset_error(args):
    sys2 = system_2()
    with pytest.raises(pylibefp.PyEFPSyntaxError):
        pylibefp.extract_subsets(sys2, **args)