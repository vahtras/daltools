import os
import numpy as np
from daltools import rspvec, sirifc

def setup():
    global tmpdir, RSPVEC, ifc
    n, e = os.path.splitext(__file__)
    tmpdir = n + ".d"
    RSPVEC = os.path.join(tmpdir, 'RSPVEC')
    SIRIFC = os.path.join(tmpdir, 'SIRIFC')
    ifc = sirifc.sirifc(name=SIRIFC)

def assert_(this, ref):
    print this, ref
    assert np.allclose(this, ref)

def test_read():
    Nx = rspvec.read("XDIPLEN", RSPVEC)
    this = Nx[12]
    ref = -0.75732690
    assert_(this, ref)

def test_tomat():
    Nx = rspvec.read("XDIPLEN", RSPVEC)
    kx = rspvec.tomat(Nx, ifc, tmpdir=tmpdir)
    this = kx[8,3]
    ref = 0.75732690
    assert_(this, ref)

def test_jwop():
    LUINDF = os.path.join(tmpdir, 'LUINDF')
    this = rspvec.jwop(LUINDF)
    ref = [ (i, j) for i in range(1, 9) for j in range(9,13)]
    assert_(this, ref)

if __name__ == "__main__":
    setup()
