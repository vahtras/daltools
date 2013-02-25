import os
import numpy as np
from daltools.lr import LR

def setup():
    global tmpdir
    n, e = os.path.splitext(__file__)
    tmpdir = n + ".d"

def test_XX():
    print tmpdir
    XX = LR('XDIPLEN', 'XDIPLEN', tmpdir)
    XXref = -2.461169664950
    assert np.allclose(XX, XXref)

def test_YY():
    print tmpdir
    YY = LR('YDIPLEN', 'YDIPLEN', tmpdir)
    YYref = -6.184500121159
    assert np.allclose(YY, YYref)

def test_YZ():
    print tmpdir
    YZ = LR('YDIPLEN', 'ZDIPLEN', tmpdir)
    YZref = 0.016158638023
    assert np.allclose(YZ, YZref)

def test_ZZ():
    print tmpdir
    ZZ = LR('ZDIPLEN', 'ZDIPLEN', tmpdir)
    ZZref = -10.308181624834
    assert np.allclose(ZZ, ZZref)

if __name__ == "__main__":
    setup()
    test_XY()
