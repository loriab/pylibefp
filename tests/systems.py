import pylibefp
#from utils import *

b2a = 0.529177
a2b = 1.0 / b2a


def system_1():
    sys = pylibefp.core.efp()

    frags = ['h2o', 'nh3']
    sys.add_potential(frags)
    sys.add_fragment(frags)
    sys.set_frag_coordinates(0, 'xyzabc', [0.0 * a2b, 0.0 * a2b, 0.0 * a2b, 1.0, 2.0, 3.0])
    sys.set_frag_coordinates(1, 'xyzabc', [5.0 * a2b, 0.0 * a2b, 0.0 * a2b, 5.0, 2.0, 8.0])

    sys.prepare()
    return sys


def system_2():
    sys = pylibefp.core.efp()

    frags = ['h2o', 'nh3', 'h2o', 'h2o', 'nh3']
    sys.add_potential(frags)
    sys.add_fragment(frags)
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
    sys.set_frag_coordinates(0, 'points', [i * a2b for i in [ -3.394,  -1.900,  -3.700, -3.524,  -1.089,  -3.147, -2.544,  -2.340,  -3.445]])
    sys.set_frag_coordinates(1, 'points', [i * a2b for i in [ -5.515,   1.083,   0.968, -5.161,   0.130,   0.813, -4.833,   1.766,   0.609]])
    sys.set_frag_coordinates(2, 'points', [i * a2b for i in [  1.848,   0.114,   0.130,  1.966,   0.674,  -0.726,  0.909,   0.273,   0.517]])
    sys.set_frag_coordinates(3, 'points', [i * a2b for i in [ -1.111,  -0.084,  -4.017, -1.941,   0.488,  -3.813, -0.292,   0.525,  -4.138]])
    sys.set_frag_coordinates(4, 'points', [i * a2b for i in [ -2.056,   0.767,  -0.301, -2.999,  -0.274,  -0.551, -1.201,   0.360,   0.258]])
    sys.set_frag_coordinates(5, 'points', [i * a2b for i in [ -0.126,  -2.228,  -0.815,  0.310,  -2.476,   0.037,  0.053,  -1.277,  -1.011]])
    sys.set_frag_coordinates(6, 'points', [i * a2b for i in [ -1.850,   1.697,   3.172, -1.050,   1.592,   2.599, -2.666,   1.643,   2.614]])
    sys.set_frag_coordinates(7, 'points', [i * a2b for i in [  1.275,  -2.447,  -4.673,  0.709,  -3.191,  -3.592,  2.213,  -1.978,  -4.343]])
    sys.set_frag_coordinates(8, 'points', [i * a2b for i in [ -5.773,  -1.738,  -0.926, -5.017,  -1.960,  -1.522, -5.469,  -1.766,   0.014]])

    sys.prepare()
    return sys


def system_4():
    sys = pylibefp.core.efp()

    frags = ['acetone', 'c2h5oh', 'c6h6', 'ccl4', 'ch3oh', 'ch4', 'cl2', 'dcm', 'dmso', 'h2', 'h2o', 'nh3']
    sys.add_potential(frags)
    sys.add_fragment(frags)
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
    sys.set_frag_coordinates(10, 'xyzabc', [ 14.0 * a2b,  12.0 * a2b,   0.0 * a2b,   0.0,   0.0,   0.0])
    sys.set_frag_coordinates(11, 'xyzabc', [ 21.0 * a2b,  12.0 * a2b,   0.0 * a2b,   0.0,   2.0,   0.0])

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
