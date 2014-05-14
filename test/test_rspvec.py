import os
import numpy as np
from nose.tools import raises
from .. import rspvec, sirifc

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
    assert np.allclose(this, ref)

def test_read():
    Nx = rspvec.read("XDIPLEN", propfile=RSPVEC)[0][0]
    this = Nx[12]
    ref = -0.75732690
    assert_(this, ref)

def test_read_w():
    Nx = rspvec.read("XDIPLEN", freqs=(0.5,), propfile=RSPVEC)[0][0]
    this = Nx[12]
    ref = -2.242435
    assert_(this, ref)

@raises(rspvec.RspVecError)
def test_read_missing():
    Nx = rspvec.read("WRONGLAB", propfile=RSPVEC)

def test_read_all():
    N1, N2 = rspvec.readall("XDIPLEN", RSPVEC)
    f1, f2 = N1[1], N2[1]
    assert_([f1, f2], [0, 0.5])
    this = (N1[0][12], N2[0][12])
    ref = (-0.75732690, -2.242435)
    assert_(this, ref)

def test_tomat():
    Nx = rspvec.read("XDIPLEN", propfile=RSPVEC)[0][0]
    kx = rspvec.tomat(Nx, ifc, tmpdir=tmpdir)
    this = kx[8, 3]
    ref = 0.75732690
    assert_(this, ref)

def test_tovec():
    Nx = rspvec.read("XDIPLEN", propfile=RSPVEC)[0][0]
    kx = rspvec.tomat(Nx, ifc, tmpdir=tmpdir)
    Nx = rspvec.tovec(kx, ifc, tmpdir=tmpdir)
    this = Nx[44]
    ref = 0.75732690
    assert_(this, ref)

def test_jwop():
    this = list(rspvec.jwop(ifc))
    ref = [ (i, j) for i in range(8) for j in range(8, 12)]
    assert_(this, ref)

if __name__ == "__main__":
    setup()
