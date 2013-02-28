import os
import numpy as np
from daltools import one
from util import blocked

def setup():
    global tmpdir
    n, e = os.path.splitext(__file__)
    tmpdir = n + ".d"

def teardown():
    pass

def assert_(this, ref):
    print this, ref
    assert np.allclose(this, ref)

def test_overlap():
    Sref = blocked.triangular([2, 1])
    Sref.subblock[0][0, 0] = 1.0
    Sref.subblock[0][1, 0] = 1.24636433    
    Sref.subblock[0][1, 1] = 2.92555541
    Sref.subblock[1][0, 0] = 1.0
    S = one.read('OVERLAP', tmpdir + '/AOONEINT')
    for s, sref in zip(S.subblock, Sref.subblock):
        assert_(s, sref)

if __name__ == "__main__":
    test_h1diag()
