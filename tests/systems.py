import pylibefp
#from utils import *

b2a = 0.529177
a2b = 1.0 / b2a


def system_1():
    sys = pylibefp.core.efp()

    frags = ['h2o', 'nh3']
    sys.add_potential(frags)
    sys.add_fragment(frags)
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

    frags = ['h2o', 'nh3', 'h2o', 'h2o', 'nh3']
    sys.add_potential(frags)
    sys.add_fragment(frags)
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

    frags = ['h2o', 'nh3', 'nh3', 'nh3', 'ch3oh', 'h2o', 'h2o', 'ch3oh', 'h2o']
    sys.add_potential(frags)
    sys.add_fragment(frags)
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

    frags = ['acetone', 'c2h5oh', 'c6h6', 'ccl4', 'ch3oh', 'ch4', 'cl2', 'dcm', 'dmso', 'h2', 'h2o', 'nh3']
    sys.add_potential(frags)
    sys.add_fragment(frags)
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


def system_5():
    sys = pylibefp.core.efp()

    frags = ['h2o', 'nh3']
    sys.add_potential(frags)
    sys.add_fragment(frags)
    sys.set_frag_coordinates(0, 'xyzabc', [   0.0 * a2b,   0.0 * a2b,   0.0 * a2b,   3.0,   0.0,   7.0])
    sys.set_frag_coordinates(1, 'xyzabc', [  18.0 * a2b,  18.0 * a2b,  18.0 * a2b,   5.0,   4.0,   6.0])

    sys.prepare()
    return sys


def system_6():
    sys = pylibefp.core.efp()

    frags = ['h2o', 'ch3oh', 'h2o', 'ch3oh', 'nh3']
    sys.add_potential(frags)
    sys.add_fragment(frags)
    sys.set_frag_coordinates(0, 'xyzabc', [   0.0 * a2b,   0.0 * a2b,   0.0 * a2b,   0.0,   0.0,   0.0])
    sys.set_frag_coordinates(1, 'xyzabc', [  19.0 * a2b,   0.0 * a2b,   0.0 * a2b,   0.0,   0.0,   0.0])
    sys.set_frag_coordinates(2, 'xyzabc', [   0.0 * a2b,  19.0 * a2b,   0.0 * a2b,   0.0,   0.0,   0.0])
    sys.set_frag_coordinates(3, 'xyzabc', [   0.0 * a2b,   0.0 * a2b,  19.0 * a2b,   0.0,   0.0,   0.0])
    sys.set_frag_coordinates(4, 'xyzabc', [  18.0 * a2b,  18.0 * a2b,  18.0 * a2b,   0.0,   0.0,   0.0])

    sys.prepare()
    return sys


def system_qm1():
    sys = pylibefp.core.efp()

    frags = ['h2o', 'c6h6', 'nh3']
    sys.add_potential(frags)
    sys.add_fragment(frags)
    sys.set_frag_coordinates(0, 'xyzabc', [  -1.6 * a2b,   4.7 * a2b,   1.4 * a2b,  -1.3,   0.1,   7.0])
    sys.set_frag_coordinates(1, 'xyzabc', [   0.4 * a2b,  -0.9 * a2b,  -0.7 * a2b,   2.3,   1.6,  -2.3])
    sys.set_frag_coordinates(2, 'xyzabc', [  -3.5 * a2b,  -2.0 * a2b,  -0.7 * a2b,   0.0,   2.2,   2.7])

    sys.prepare()
    return sys


def system_qm2():
    sys = pylibefp.core.efp()

    frags = ['ch3oh', 'dmso', 'dmso', 'acetone', 'dcm', 'acetone', 'acetone']
    sys.add_potential(frags)
    sys.add_fragment(frags)
    sys.set_frag_coordinates(0, 'xyzabc', [   0.0 * a2b,  -1.0 * a2b,   0.0 * a2b,   0.0,   1.1,   2.0])
    sys.set_frag_coordinates(1, 'xyzabc', [  -5.0 * a2b,  12.0 * a2b,   0.0 * a2b,   3.0,   0.2,   5.0])
    sys.set_frag_coordinates(2, 'xyzabc', [   0.0 * a2b,  -3.0 * a2b,   5.0 * a2b,   6.0,   2.3,   8.0])
    sys.set_frag_coordinates(3, 'xyzabc', [  -5.0 * a2b,  -4.0 * a2b,  -5.0 * a2b,   9.0,   0.4,   1.0])
    sys.set_frag_coordinates(4, 'xyzabc', [  -9.0 * a2b,  -5.0 * a2b,   1.0 * a2b,   2.0,   1.5,   4.0])
    sys.set_frag_coordinates(5, 'xyzabc', [   7.0 * a2b,  -2.0 * a2b,  11.0 * a2b,   5.0,   0.6,   7.0])
    sys.set_frag_coordinates(6, 'xyzabc', [  -9.0 * a2b,  -7.0 * a2b,  -9.0 * a2b,   8.0,   2.7,   0.0])

    sys.prepare()
    return sys
