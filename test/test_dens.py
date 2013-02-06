import numpy as np
from daltools import dens

def setup():
    pass

def teardown():
    pass

def assert_(this, ref):
    print this, ref
    assert np.allclose(this, ref)

def test_h1diag():
    #D1 = dens.h1diag(1,1)
    assert True

if __name__ == "__main__":
    test_h1diag()
