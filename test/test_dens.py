import os
import numpy as np
from ..util import full
from .. import dens

def setup():
    global tmpdir
    n, e = os.path.splitext(__file__)
    tmpdir = n + ".d"
    print tmpdir

def teardown():
    pass

def assert_(this, ref):
    print this, ref
    assert np.allclose(this, ref)

def test_h1diag():
    diref = full.init(
            [[ 0.53709620, 0.33768834, 0.],
             [ 0.33768834, 0.21231469, 0.],
             [ 0.        ,  0.       , 0.]]
            )

    di, da = dens.h1diag(1, 1, filename=os.path.join(tmpdir, "AOONEINT"))
    assert_(di, diref)

def test_Dab():
    na = 2; nb = 1
    C = full.unit(3)
    da, db = dens.C2Dab(C, na, nb)
    daref = full.init(
        [[1., 0., 0], 
        [0., 1., 0.], 
        [0., 0., 0.]]
        )
    dbref = full.init(
        [[1., 0., 0], 
        [0., 0., 0.], 
        [0., 0., 0.]]
        )
    assert_(da, daref)
    assert_(db, dbref)
    

if __name__ == "__main__":
    setup()
    test_h1diag()
