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
    print this
    print ref
    try:
        assert np.allclose(this, ref)
    except AssertionError, e:
        from pdb import set_trace; set_trace()

def test_read():
    Nx = rspvec.read("XDIPLEN", propfile=RSPVEC)[0][0]
    this = Nx[5]
    ref = -2.34009730
    assert_(this, ref)

def test_tomat():
    Nx = rspvec.read("XDIPLEN", propfile=RSPVEC)[0][0]
    kx = rspvec.tomat(Nx, ifc, tmpdir=tmpdir)
    this = kx[5, 7]
    ref = -2.34009730
    assert_(this, ref)

def test_tovec():
    ref = rspvec.read("XDIPLEN", propfile=RSPVEC)[0][0]
    kx = rspvec.tomat(ref, ifc, tmpdir=tmpdir)
    this = rspvec.tovec(kx, ifc, tmpdir=tmpdir)
    assert_(this, ref)

def test_jwop():
    this = list(rspvec.jwop(ifc))
    ref = [ (i, j) for i in range(7) for j in range(7, 8)] + \
          [ (i, j) for i in range(7) for j in range(8, 11)] + \
          [ (i, j) for i in range(7,8) for j in range(8, 11)]
    assert_(this, ref)

if __name__ == "__main__":
    setup()
