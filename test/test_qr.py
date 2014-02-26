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
    print tmpdir
    XXX = QR('XDIPLEN', 'XDIPLEN', 'XDIPLEN', tmpdir=tmpdir, a2test=True)
    XXXref = 9.652e-05
    assert_(XXX, XXXref)



def test_XXX():
    print tmpdir
    XX = QR('XDIPLEN', 'XDIPLEN', 'XDIPLEN', tmpdir=tmpdir)
    XXref = 0.00013323
    assert_(XX, XXref)


if __name__ == '__main__':
    setup()
    test_XXX_A2()

#+++ okay decompyling test_qr.pyc 
# decompiled 1 files: 1 okay, 0 failed, 0 verify failed
# 2014.02.26 13:49:28 CET
