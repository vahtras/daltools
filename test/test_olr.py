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
    XX = LR('XDIPLEN', 'XDIPLEN', tmpdir)
    XXref = -5.455606903637
    assert_(XX, XXref)

def test_YY():
    print tmpdir
    YY = LR('YDIPLEN', 'YDIPLEN', tmpdir)
    YYref = -10.31180304740
    assert_(YY, YYref)

def test_YZ():
    print tmpdir
    YZ = LR('YDIPLEN', 'ZDIPLEN', tmpdir)
    YZref = 0.2324294799056
    assert_(YZ, YZref)

def test_ZZ():
    print tmpdir
    ZZ = LR('ZDIPLEN', 'ZDIPLEN', tmpdir)
    ZZref = -3.124432068117
    assert_(ZZ, ZZref)

if __name__ == "__main__":
    setup()
