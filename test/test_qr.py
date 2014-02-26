# 2014.02.26 13:49:28 CET
import os
import numpy as np
from daltools.qr import QR

def assert_(this, ref):
    print this
    print ref
    assert np.allclose(this, ref)



def setup():
    global tmpdir
    (n, e,) = os.path.splitext(__file__)
    tmpdir = n + '.d'



def test_XXX_A2():
    XXX = QR('XDIPLEN', 'XDIPLEN', 'XDIPLEN', tmpdir=tmpdir, a2test=True)
    XXXref = 9.652e-05
    assert_(XXX, XXXref)



def test_XXX():
    XXX = QR('XDIPLEN', 'XDIPLEN', 'XDIPLEN', tmpdir=tmpdir)
    XXXref = 0.00013323
    assert_(XXX, XXXref)

def test_ZXX():
    ZXX = QR('ZDIPLEN', 'XDIPLEN', 'XDIPLEN', tmpdir=tmpdir)
    ZXXref = -1.68075251
    assert_(ZXX, ZXXref)


if __name__ == '__main__':
    setup()
    test_XXX_A2()

#+++ okay decompyling test_qr.pyc 
# decompiled 1 files: 1 okay, 0 failed, 0 verify failed
# 2014.02.26 13:49:28 CET
