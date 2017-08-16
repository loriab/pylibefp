import sys
import pytest
import pprint
import pylibefp
from utils import *
from systems import *

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
    print(sys1p.print_geometry(units_to_bohr=b2a))
    print(sys1p.print_geometry(units_to_bohr=1.0))

    assert(compare_integers(2, sys1p.get_frag_count(), sys._getframe().f_code.co_name + ': nfrag'))
    assert(compare_values(0.0001922903, ene['total'], 6,  sys._getframe().f_code.co_name + ': ene'))


    
