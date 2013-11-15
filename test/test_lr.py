import os
import numpy as np
from daltools.lr import LR

def assert_(this, ref):
    print this
    print ref
    assert np.allclose(this, ref)

def setup():
    global tmpdir
    n, e = os.path.splitext(__file__)
    tmpdir = n + ".d"

def test_XX():
    print tmpdir
    XX = LR('XDIPLEN', 'XDIPLEN', 0, tmpdir)
    XXref = -2.461169664950
    assert_(XX, XXref)

def test_YY():
    print tmpdir
    YY = LR('YDIPLEN', 'YDIPLEN', 0, tmpdir)
    YYref = -6.184500121159
    assert_(YY, YYref)

def test_YZ():
    print tmpdir
    YZ = LR('YDIPLEN', 'ZDIPLEN', 0, tmpdir)
    YZref = 0.016158638023
    assert_(YZ, YZref)

def test_ZZ():
    print tmpdir
    ZZ = LR('ZDIPLEN', 'ZDIPLEN', 0, tmpdir)
    ZZref = -10.308181624834
    assert_(ZZ, ZZref)

if __name__ == "__main__":
    setup()
    test_XX()
