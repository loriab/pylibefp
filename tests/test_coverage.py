import sys
import pytest
import pylibefp
from utils import *
from systems import *


def test_grad_fail():
    asdf = system_1()
    asdf.compute(do_gradient=False)

    with pytest.raises(pylibefp.Fatal) as e_info:
        grad = asdf.get_gradient()

