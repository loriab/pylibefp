import sys
import pylibefp
from testing_utils import *

b2a = 0.529177
a2b = 1.0 / b2a


def system_1():
    sys = pylibefp.core.efp()
    sys.create()
    
    frags = ['h2o', 'nh3']
    sys.add_potentials(frags)
    sys.add_fragments(frags)
    #sys.add_potential('../../fraglib/h2o.efp')
    #sys.add_potential('../../fraglib/nh3.efp')
    #sys.add_fragment('h2o_l')
    sys.set_frag_coordinates(0, 'xyzabc', [0.0 * a2b, 0.0 * a2b, 0.0 * a2b, 1.0, 2.0, 3.0])
    #sys.add_fragment('nh3_l')
    sys.set_frag_coordinates(1, 'xyzabc', [5.0 * a2b, 0.0 * a2b, 0.0 * a2b, 5.0, 2.0, 8.0])

    sys.prepare()
    return sys


def system_2():
    sys = pylibefp.core.efp()
    sys.create()
    
    frags = ['h2o', 'nh3', 'h2o', 'h2o', 'nh3']
    sys.add_potentials(frags)
    sys.add_fragments(frags)
    #sys.add_potential('../../fraglib/h2o.efp')
    #sys.add_potential('../../fraglib/nh3.efp')
    #sys.add_fragment('h2o_l')
    #sys.add_fragment('nh3_l')
    #sys.add_fragment('h2o_l')
    #sys.add_fragment('h2o_l')
    #sys.add_fragment('nh3_l')
    sys.set_frag_coordinates(0, 'xyzabc', [-1.0 * a2b,   3.7 * a2b,   0.4 * a2b,  -1.3,   0.0,   7.0])
    sys.set_frag_coordinates(1, 'xyzabc', [ 0.4 * a2b,  -0.9 * a2b,  -0.7 * a2b,   4.0,   1.6,  -2.3])
    sys.set_frag_coordinates(2, 'xyzabc', [ 1.7 * a2b,   2.0 * a2b,   3.3 * a2b,  -1.2,  -2.0,   6.2])
    sys.set_frag_coordinates(3, 'xyzabc', [ 0.0 * a2b,   3.9 * a2b,  -3.4 * a2b,   1.3,   5.2,  -3.0])
    sys.set_frag_coordinates(4, 'xyzabc', [-3.5 * a2b,   0.0 * a2b,  -0.7 * a2b,   0.0,  -2.7,   2.7])

    sys.prepare()
    return sys


def system_3():
    sys = pylibefp.core.efp()
    sys.create()

    frags = ['h2o', 'nh3', 'nh3', 'nh3', 'ch3oh', 'h2o', 'h2o', 'ch3oh', 'h2o']
    sys.add_potentials(frags)
    sys.add_fragments(frags)
    #sys.add_potential('../../fraglib/nh3.efp')
    #sys.add_potential('../../fraglib/ch3oh.efp')
    #sys.add_potential('../../fraglib/h2o.efp')
    #sys.add_fragment('h2o_l')
    #sys.add_fragment('nh3_l')
    #sys.add_fragment('nh3_l')
    #sys.add_fragment('nh3_l')
    #sys.add_fragment('ch3oh_l')
    #sys.add_fragment('h2o_l')
    #sys.add_fragment('h2o_l')
    #sys.add_fragment('ch3oh_l')
    #sys.add_fragment('h2o_l')
    sys.set_frag_coordinates(0, 'points', [  -3.394 * a2b,  -1.900 * a2b,  -3.700 * a2b, -3.524 * a2b,  -1.089 * a2b,  -3.147 * a2b, -2.544 * a2b,  -2.340 * a2b,  -3.445 * a2b])
    sys.set_frag_coordinates(1, 'points', [  -5.515 * a2b,   1.083 * a2b,   0.968 * a2b, -5.161 * a2b,   0.130 * a2b,   0.813 * a2b, -4.833 * a2b,   1.766 * a2b,   0.609 * a2b])
    sys.set_frag_coordinates(2, 'points', [   1.848 * a2b,   0.114 * a2b,   0.130 * a2b,  1.966 * a2b,   0.674 * a2b,  -0.726 * a2b,  0.909 * a2b,   0.273 * a2b,   0.517 * a2b])
    sys.set_frag_coordinates(3, 'points', [  -1.111 * a2b,  -0.084 * a2b,  -4.017 * a2b, -1.941 * a2b,   0.488 * a2b,  -3.813 * a2b, -0.292 * a2b,   0.525 * a2b,  -4.138 * a2b])
    sys.set_frag_coordinates(4, 'points', [  -2.056 * a2b,   0.767 * a2b,  -0.301 * a2b, -2.999 * a2b,  -0.274 * a2b,  -0.551 * a2b, -1.201 * a2b,   0.360 * a2b,   0.258 * a2b])
    sys.set_frag_coordinates(5, 'points', [  -0.126 * a2b,  -2.228 * a2b,  -0.815 * a2b,  0.310 * a2b,  -2.476 * a2b,   0.037 * a2b,  0.053 * a2b,  -1.277 * a2b,  -1.011 * a2b])
    sys.set_frag_coordinates(6, 'points', [  -1.850 * a2b,   1.697 * a2b,   3.172 * a2b, -1.050 * a2b,   1.592 * a2b,   2.599 * a2b, -2.666 * a2b,   1.643 * a2b,   2.614 * a2b])
    sys.set_frag_coordinates(7, 'points', [   1.275 * a2b,  -2.447 * a2b,  -4.673 * a2b,  0.709 * a2b,  -3.191 * a2b,  -3.592 * a2b,  2.213 * a2b,  -1.978 * a2b,  -4.343 * a2b])
    sys.set_frag_coordinates(8, 'points', [  -5.773 * a2b,  -1.738 * a2b,  -0.926 * a2b, -5.017 * a2b,  -1.960 * a2b,  -1.522 * a2b, -5.469 * a2b,  -1.766 * a2b,   0.014 * a2b])

    sys.prepare()
    return sys


def system_4():
    sys = pylibefp.core.efp()
    sys.create()

    frags = ['acetone', 'c2h5oh', 'c6h6', 'ccl4', 'ch3oh', 'ch4', 'cl2', 'dcm', 'dmso', 'h2', 'h2o', 'nh3']
    sys.add_potentials(frags)
    sys.add_fragments(frags)
    #for fr in ['acetone', 'c2h5oh', 'c6h6', 'ccl4', 'ch3oh', 'ch4', 'cl2', 'dcm', 'dmso', 'h2', 'h2o', 'nh3']:
    #    sys.add_potential('../../fraglib/' + fr + '.efp')
    #sys.add_fragment('acetone_l')
    #sys.add_fragment('c2h5oh_l')
    #sys.add_fragment('c6h6_l')
    #sys.add_fragment('ccl4_l')
    #sys.add_fragment('ch3oh_l')
    #sys.add_fragment('ch4_l')
    #sys.add_fragment('cl2_l')
    #sys.add_fragment('dcm_l')
    #sys.add_fragment('dmso_l')
    #sys.add_fragment('h2_l')
    #sys.add_fragment('h2o_l')
    #sys.add_fragment('nh3_l')
    sys.set_frag_coordinates(0, 'xyzabc', [   0.0 * a2b,   0.0 * a2b,   0.0 * a2b,   0.0,   0.2,   0.3])
    sys.set_frag_coordinates(1, 'xyzabc', [   7.0 * a2b,   0.0 * a2b,   0.0 * a2b,   0.0,   2.0,   3.7])
    sys.set_frag_coordinates(2, 'xyzabc', [  14.0 * a2b,   0.0 * a2b,   0.0 * a2b,   3.1,   0.8,   2.0])
    sys.set_frag_coordinates(3, 'xyzabc', [  21.0 * a2b,   0.0 * a2b,   0.0 * a2b,   0.0,   8.0,   0.0])
    sys.set_frag_coordinates(4, 'xyzabc', [   0.0 * a2b,   6.0 * a2b,   0.0 * a2b,   0.7,   2.0,   1.0])
    sys.set_frag_coordinates(5, 'xyzabc', [   7.0 * a2b,   6.0 * a2b,   0.0 * a2b,   0.6,   0.0,   4.7])
    sys.set_frag_coordinates(6, 'xyzabc', [  14.0 * a2b,   6.0 * a2b,   0.0 * a2b,   0.0,   0.0,   0.3])
    sys.set_frag_coordinates(7, 'xyzabc', [  21.0 * a2b,   6.0 * a2b,   0.0 * a2b,   0.0,   0.4,   0.3])
    sys.set_frag_coordinates(8, 'xyzabc', [   0.0 * a2b,  12.0 * a2b,   0.0 * a2b,   0.8,   0.0,   0.0])
    sys.set_frag_coordinates(9, 'xyzabc', [   7.0 * a2b,  12.0 * a2b,   0.0 * a2b,   8.0,   0.7,   0.8])
    sys.set_frag_coordinates(10, 'xyzabc', [  14.0 * a2b,  12.0 * a2b,   0.0 * a2b,   0.0,   0.0,   0.0])
    sys.set_frag_coordinates(11, 'xyzabc', [  21.0 * a2b,  12.0 * a2b,   0.0 * a2b,   0.0,   2.0,   0.0])

    sys.prepare()
    return sys


def test_elec_1a():
    asdf = system_1()
    asdf.set_opts({'elec': True, 'elec_damp': 'screen'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0002900482, ene['elec'], 6, sys._getframe().f_code.co_name))


def test_elec_1b():

    asdf = system_1()
    asdf.set_opts({'elec': True, 'elec_damp': 'overlap'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0002910961, ene['elec'], 6, sys._getframe().f_code.co_name))


def test_pol_1a():
    asdf = system_1()
    asdf.set_opts({'pol': True, 'elec_damp': 'screen'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0002777238, ene['pol'], 6, sys._getframe().f_code.co_name))


def test_pol_1b():

    asdf = system_1()
    asdf.set_opts({'pol': True, 'elec_damp': 'screen', 'elec': True, 'pol_driver': 'direct'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0002777238, ene['pol'], 6, sys._getframe().f_code.co_name))


def test_disp_1a():
    asdf = system_1()
    asdf.set_opts({'disp': True, 'disp_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0000989033, ene['disp'], 6, sys._getframe().f_code.co_name))


def test_disp_1b():

    asdf = system_1()
    asdf.set_opts({'disp': True, 'disp_damp': 'overlap'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0001007275, ene['disp'], 6, sys._getframe().f_code.co_name))


def test_xr_1():
    asdf = system_1()
    asdf.set_opts({'xr': True})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0000134716, ene['xr'], 6, sys._getframe().f_code.co_name))


def test_total_1a():
    asdf = system_1()
    asdf.set_opts({'elec': True, 'elec_damp': 'screen',
                     'xr': True,
                    'pol': True, # 'pol_damp': 'tt',
                   'disp': True, 'disp_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
    print('<<< get_opts():  ', asdf.get_opts(), '>>>')
    print('<<< summary():   ', asdf.summary(), '>>>')
    print('<<< get_energy():', ene, '>>>')
    print('<<< get_atoms(): ', asdf.get_atoms(), '>>>')
    print(asdf.print_geometry(units_to_bohr=b2a))
    assert(compare_values(0.0001922903, ene['total'], 6, sys._getframe().f_code.co_name))


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
    asdf.set_opts({'pol': True, 'elec_damp': 'screen'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0013685212, ene['pol'], 6, sys._getframe().f_code.co_name))


def test_pol_2b():
    asdf = system_2()
    asdf.set_opts({'pol': True, 'elec_damp': 'screen', 'pol_driver': 'direct'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(0.0013685212, ene['pol'], 6, sys._getframe().f_code.co_name))


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
    assert(compare_values(-0.0066095992, ene['pol'], 6, sys._getframe().f_code.co_name))


def test_pol_3b():
    asdf = system_3()
    asdf.set_opts({'pol': True, 'elec_damp': 'screen', 'pol_damp': 'off', 'pol_driver': 'direct'})
    asdf.compute()
    ene = asdf.get_energy()
    assert(compare_values(-0.0066095992, ene['pol'], 6, sys._getframe().f_code.co_name))


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
    assert(compare_values(0.0061408841, ene['total'], 5, sys._getframe().f_code.co_name))


def test_total_4a():
    asdf = system_4()
    asdf.set_opts({'elec': True, 'pol': True, 'disp': True, 'xr': True, 'elec_damp': 'screen', 'disp_damp': 'tt', 'pol_damp': 'tt'})
    asdf.compute()
    ene = asdf.get_energy()
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

