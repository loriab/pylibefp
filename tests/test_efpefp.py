import sys
import pytest
import pprint
import pylibefp
from utils import *
from systems import *


def blank_ene():
    fields = ['charge_penetration', 'disp', 'dispersion', 'elec',
              'electrostatic', 'electrostatic_point_charges',
              'exchange_repulsion', 'pol', 'polarization', 'xr']
    ene = {f: 0.0 for f in fields}
    return ene


def test_elec_1a():
    asdf = system_1()
    asdf.set_opts({'elec': True, 'elec_damp': 'screen'})
    asdf.compute()
    ene = asdf.get_energy()

    expected_ene = blank_ene()
    expected_ene['elec'] = expected_ene['electrostatic'] = expected_ene['total'] = 0.0002900482
    assert(compare_dicts(expected_ene, ene, 6,  sys._getframe().f_code.co_name + ': ene'))


def test_elec_1b():
    asdf = system_1()
    asdf.set_opts({'elec': True, 'elec_damp': 'overlap'})
    asdf.compute()
    ene = asdf.get_energy()

    elst = 0.0002910961
    cp = -8.066354689359154e-07
    expected_ene = blank_ene()
    expected_ene['elec'] = expected_ene['total'] = elst
    expected_ene['charge_penetration'] = cp
    expected_ene['electrostatic'] = elst - cp
    assert(compare_dicts(expected_ene, ene, 6,  sys._getframe().f_code.co_name + ': ene'))


def test_pol_1a():
    asdf = system_1()
    opts = {'elec': True, 'pol': True, 'elec_damp': 'screen'}
    asdf.set_opts(opts)
    asdf.compute()
    ene = asdf.get_energy()

    elec = 0.0002900482
    pol = 0.0002777238 - elec
    expected_ene = blank_ene()
    expected_ene['elec'] = expected_ene['electrostatic'] = elec
    expected_ene['pol'] = expected_ene['polarization'] = pol
    expected_ene['total'] = elec + pol
    pprint.pprint(opts)
    assert(compare_dicts(expected_ene, ene, 6,  sys._getframe().f_code.co_name + ': ene'))


def test_pol_1b():
    asdf = system_1()
    asdf.set_opts({'pol': True, 'elec_damp': 'screen', 'elec': True, 'pol_driver': 'direct'})
    asdf.compute()
    ene = asdf.get_energy()

    elec = 0.0002900478
    pol = 0.0002777238 - elec
    expected_ene = blank_ene()
    expected_ene['elec'] = expected_ene['electrostatic'] = elec
    expected_ene['pol'] = expected_ene['polarization'] = pol
    expected_ene['total'] = elec + pol
    assert(compare_dicts(expected_ene, ene, 6,  sys._getframe().f_code.co_name + ': ene'))


def test_disp_1a():
    asdf = system_1()
    asdf.set_opts({'disp': True, 'disp_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()

    expected_ene = blank_ene()
    expected_ene['disp'] = expected_ene['dispersion'] = expected_ene['total'] = -0.0000989033
    assert(compare_dicts(expected_ene, ene, 6,  sys._getframe().f_code.co_name + ': ene'))


def test_disp_1b():

    asdf = system_1()
    asdf.set_opts({'disp': True, 'disp_damp': 'overlap'})
    asdf.compute()
    ene = asdf.get_energy()

    expected_ene = blank_ene()
    expected_ene['disp'] = expected_ene['dispersion'] = expected_ene['total'] = -0.0001007275
    assert(compare_dicts(expected_ene, ene, 6,  sys._getframe().f_code.co_name + ': ene'))


def test_xr_1():
    asdf = system_1()
    asdf.set_opts({'xr': True})
    asdf.compute()
    ene = asdf.get_energy()

    expected_ene = blank_ene()
    expected_ene['xr'] = expected_ene['exchange_repulsion'] = expected_ene['total'] = 0.0000134716
    assert(compare_dicts(expected_ene, ene, 6,  sys._getframe().f_code.co_name + ': ene'))


def test_total_1a():
    asdf = system_1()
    asdf.set_opts({'elec': True, 'elec_damp': 'screen',
                     'xr': True,
                    'pol': True, # 'pol_damp': 'tt',
                   'disp': True, 'disp_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
    pprint.pprint(ene)
    print('<<< get_opts():  ', asdf.get_opts(), '>>>')
    #print('<<< summary():   ', asdf.summary(), '>>>')
    print('<<< get_energy():', ene, '>>>')
    print('<<< get_atoms(): ', asdf.get_atoms(), '>>>')
    print(asdf.energy_summary())
    print(asdf.print_geometry(units_to_bohr=b2a))
    print(asdf.print_geometry(units_to_bohr=1.0))

    expected_ene = blank_ene()
    expected_ene['elec'] = expected_ene['electrostatic'] = 0.0002900482
    expected_ene['xr'] = expected_ene['exchange_repulsion'] = 0.0000134716
    expected_ene['pol'] = expected_ene['polarization'] = 0.0002777238 - expected_ene['electrostatic']
    expected_ene['disp'] = expected_ene['dispersion'] = -0.0000989033
    expected_ene['total'] = 0.0001922903
    assert(compare_integers(2, asdf.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag'))
    assert(compare_dicts(expected_ene, ene, 6,  sys._getframe().f_code.co_name + ': ene'))


def test_elec_2a():
    asdf = system_2()
    asdf.set_opts({'elec': True, 'elec_damp': 'screen'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0015865516, ene['elec'], 6, sys._getframe().f_code.co_name))


def test_elec_2b():
    asdf = system_2()
    asdf.set_opts({'elec': True, 'elec_damp': 'overlap'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0017049246, ene['elec'], 6, sys._getframe().f_code.co_name))


def test_pol_2a():
    asdf = system_2()
    asdf.set_opts({'elec': True, 'pol': True, 'elec_damp': 'screen'})
    asdf.compute()
    ene = asdf.get_energy()
    pprint.pprint(ene)
    assert(compare_values(0.0013685212, ene['total'], 6, sys._getframe().f_code.co_name))


def test_pol_2b():
    asdf = system_2()
    asdf.set_opts({'elec': True, 'pol': True, 'elec_damp': 'screen', 'pol_driver': 'direct'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0013685212, ene['total'], 6, sys._getframe().f_code.co_name))


def test_disp_2a():
    asdf = system_2()
    asdf.set_opts({'disp': True, 'disp_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0014688094, ene['disp'], 6, sys._getframe().f_code.co_name))


def test_disp_2b():
    asdf = system_2()
    asdf.set_opts({'disp': True, 'disp_damp': 'overlap'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0015801770, ene['disp'], 6, sys._getframe().f_code.co_name))


def test_xr_2():
    asdf = system_2()
    asdf.set_opts({'xr': True})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0008443933, ene['xr'], 6, sys._getframe().f_code.co_name))


def test_total_2a():
    asdf = system_2()
    asdf.set_opts({'elec': True, 'pol': True, 'disp': True, 'xr': True, 'elec_damp': 'screen', 'disp_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_integers(5, asdf.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag'))
    assert(compare_values(0.0007440865, ene['total'], 6, sys._getframe().f_code.co_name))


def test_elec_3a():
    asdf = system_3()
    asdf.set_opts({'elec': True, 'elec_damp': 'screen'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0039531505, ene['elec'], 6, sys._getframe().f_code.co_name))


def test_elec_3b():
    asdf = system_3()
    asdf.set_opts({'elec': True, 'elec_damp': 'overlap'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0023592829, ene['elec'], 6, sys._getframe().f_code.co_name))


def test_pol_3a():
    asdf = system_3()
    asdf.set_opts({'elec': True, 'pol': True, 'elec_damp': 'screen', 'pol_damp': 'off'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0066095992, ene['total'], 6, sys._getframe().f_code.co_name))


def test_pol_3b():
    asdf = system_3()
    asdf.set_opts({'elec': True, 'pol': True, 'elec_damp': 'screen', 'pol_damp': 'off', 'pol_driver': 'direct'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0066095992, ene['total'], 6, sys._getframe().f_code.co_name))


def test_disp_3a():
    asdf = system_3()
    asdf.set_opts({'disp': True, 'disp_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0173897265, ene['disp'], 6, sys._getframe().f_code.co_name))


def test_disp_3b():
    asdf = system_3()
    asdf.set_opts({'disp': True, 'disp_damp': 'overlap'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0220107872, ene['disp'], 6, sys._getframe().f_code.co_name))


def test_xr_3():
    asdf = system_3()
    asdf.set_opts({'xr': True})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0301402098, ene['xr'], 5, sys._getframe().f_code.co_name))


def test_total_3a():
    asdf = system_3()
    asdf.set_opts({'elec': True, 'pol': True, 'disp': True, 'xr': True, 'elec_damp': 'screen', 'disp_damp': 'tt', 'pol_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_integers(9, asdf.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag'))
    assert(compare_values(0.0061408841, ene['total'], 5, sys._getframe().f_code.co_name))


def test_total_4a():
    asdf = system_4()
    asdf.set_opts({'elec': True, 'pol': True, 'disp': True, 'xr': True, 'elec_damp': 'screen', 'disp_damp': 'tt', 'pol_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_integers(12, asdf.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag'))
    assert(compare_values(-0.0095597483, ene['total'], 5, sys._getframe().f_code.co_name))


def test_total_4b():
    asdf = system_4()
    asdf.set_opts({'elec': True, 'pol': True, 'disp': True, 'xr': True, 'elec_damp': 'overlap', 'disp_damp': 'overlap', 'pol_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0092400662, ene['total'], 5, sys._getframe().f_code.co_name))


def test_total_4c():
    asdf = system_4()
    asdf.set_opts({'elec': True, 'pol': True, 'disp': True, 'xr': True, 'elec_damp': 'off', 'disp_damp': 'off', 'pol_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0091278725, ene['total'], 5, sys._getframe().f_code.co_name))


def test_total_4d():
    asdf = system_4()
    asdf.set_opts({'elec': True, 'pol': True, 'disp': True, 'xr': True, 'elec_damp': 'screen', 'disp_damp': 'tt', 'pol_damp': 'tt', 'pol_driver': 'direct'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0095597483, ene['total'], 5, sys._getframe().f_code.co_name))


def test_efpefptorque():
    asdf = system_3()
    asdf.set_opts({'disp': False, 'exch': False, 'elst_damp': 'screen', 'ind_damp': 'off'}, label='psi', append='psi')  # V equiv
    #asdf.set_opts({'elec': True, 'pol': True, 'elec_damp': 'screen', 'pol_damp': 'off', 'dertype': 'first'})            # ^ equiv
    asdf.compute(do_gradient=True)
    ene = asdf.get_energy()
    torq = asdf.get_gradient()

    ref = {'torque': [
        -0.0014557485,    -0.0024650113,    -0.0007420245,     0.0018487317,    -0.0065430367,    -0.0003612802,
        -0.0024798509,    -0.0002766252,     0.0029343456,    -0.0033124877,    -0.0048014449,    -0.0046442270,
         0.0021341431,     0.0023700691,     0.0015655930,    -0.0005188401,    -0.0004406075,    -0.0016388193,
        -0.0020017801,     0.0045394287,     0.0001140076,    -0.0011159049,     0.0021766586,    -0.0035556589,
        -0.0004997047,     0.0037416773,    -0.0017226579,     0.0108138324,     0.0056465424,    -0.0031926302,
        -0.0004161161,    -0.0046891120,    -0.0017098053,    -0.0023800599,     0.0042322597,     0.0105675357,
         0.0007828963,     0.0001744122,    -0.0006861146,     0.0003752826,    -0.0032331154,    -0.0011471607,
         0.0038830634,    -0.0039883720,    -0.0001194227,     0.0012427711,    -0.0026362462,    -0.0005023332,
         0.0000530976,     0.0005935332,     0.0003660789,    -0.0015382262,    -0.0048146666,     0.0026841256 ]}

    assert(compare_values(-0.0066095987170644, ene['total'], 5, sys._getframe().f_code.co_name + ': ene'))
    assert(compare_dicts(ref, {'torque': torq}, 6, sys._getframe().f_code.co_name + ': torq'))


def test_efpefp_bz2():
    """psi4/test/libefp/qchem-efp-sp"""
    b2a = 0.529177
    a2b = 1.0 / b2a

    asdf = pylibefp.core.efp()
    asdf.create()

    frags = ['c6h6', 'c6h6']
    asdf.add_potential(frags)
    asdf.add_fragment(frags)
    asdf.set_frag_coordinates(0, 'xyzabc', [-0.30448173 * a2b, -2.24210052 * a2b, -0.29383131 * a2b, -0.642499, 1.534222, -0.568147])
    asdf.set_frag_coordinates(1, 'xyzabc', [-0.60075437 * a2b,  1.36443336 * a2b,  0.78647823 * a2b,  3.137879, 1.557344, -2.568550])
    asdf.prepare()

    asdf.set_opts({'disp_damp': 'tt'}, append='psi')
    asdf.compute()
    ene = asdf.get_energy(label='psi')

    # values copied from q-chem output file
    assert(compare_values(-0.006945881265, ene['elst'], 6, sys._getframe().f_code.co_name + ': ene elst'))
    assert(compare_values( 0.046915489574, ene['exch'], 6, sys._getframe().f_code.co_name + ': ene exch'))
    assert(compare_values(-0.000675030191, ene['ind'], 6, sys._getframe().f_code.co_name + ': ene ind'))
    assert(compare_values(-0.021092526180, ene['disp'], 6, sys._getframe().f_code.co_name + ': ene disp'))
    assert(compare_values( 0.018202051938, ene['total'], 6, sys._getframe().f_code.co_name + ': ene'))

