import os
import numpy as np
from daltools import one
from util import blocked

def setup():
    global tmpdir
    n, e = os.path.splitext(__file__)
    tmpdir = n + ".d"
    global aooneint
    aooneint = os.path.join(tmpdir, 'AOONEINT')

def teardown():
    pass

def assert_(this, ref):
    print this, ref
    assert np.allclose(this, ref)

def test_header():
    head = one.readhead(aooneint)
    assert "B term (MCD) components of H3+" in head["ttitle"]
    assert head["naos"] == (2, 1, 0, 0)
    assert head["nsym"] == 4
    assert_(head["potnuc"], 1.82903817207)

def test_isordk():
    isordk = one.readisordk(aooneint)
    assert isordk["nucdep"] == 3
    #assert isordk["mxcent"] == 120
    assert_(isordk["chrn"][:3], (1., 1., 1.))
    assert_(isordk["cooo"][0::120], (0, -0.224906, 0))
    assert_(isordk["cooo"][1::120], (1, 0.899624, 0))
    assert_(isordk["cooo"][2::120], (-1, 0.899624, 0))

def test_scfinp():
    scfinp = one.readscfinp(aooneint)
    assert scfinp["nsym"] == 4
    assert_(scfinp["cooo"], (0.0, -0.224905893, 0.0, 1.0, 0.899623572, 0.0, -1.0, 0.899623572, 0.0))

def test_overlap():
    Sref = blocked.triangular([2, 1])
    Sref.subblock[0][0, 0] = 1.0
    Sref.subblock[0][1, 0] = 1.24636433    
    Sref.subblock[0][1, 1] = 2.92555541
    Sref.subblock[1][0, 0] = 1.0
    S = one.read("OVERLAP", aooneint)
    for s, sref in zip(S.subblock, Sref.subblock):
        assert_(s, sref)

if __name__ == "__main__":
    setup()
    test_overlap()
