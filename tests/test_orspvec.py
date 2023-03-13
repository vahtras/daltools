import numpy as np
from pytest import approx

from daltools import rspvec, sirifc

from . import tmpdir

class TestOpenRspVec:

    def setup_method(self):
        self.tmpdir = tmpdir(__file__)
        self.RSPVEC = self.tmpdir/"RSPVEC"
        self.SIRIFC = self.tmpdir/"SIRIFC"
        self.ifc = sirifc.sirifc(name=self.SIRIFC)

    def test_read(self):
        Nx = rspvec.read("XDIPLEN", propfile=self.RSPVEC)["XDIPLEN"]
        this = Nx[5]
        ref = -2.34009730
        assert this == approx(ref)

    def test_tomat(self):
        Nx = rspvec.read("XDIPLEN", propfile=self.RSPVEC)["XDIPLEN"]
        kx = rspvec.tomat(Nx, self.ifc, tmpdir=self.tmpdir)
        this = kx[5, 7]
        ref = -2.34009730
        assert this == approx(ref)

    def test_tovec(self):
        ref = rspvec.read("XDIPLEN", propfile=self.RSPVEC)["XDIPLEN"]
        kx = rspvec.tomat(ref, self.ifc, tmpdir=self.tmpdir)
        this = rspvec.tovec(kx, self.ifc, tmpdir=self.tmpdir)
        np.testing.assert_almost_equal(this, ref)

    def test_jwop(self):
        this = list(rspvec.jwop(self.ifc))
        ref = (
            [(i, j) for i in range(7) for j in range(7, 8)]
            + [(i, j) for i in range(7) for j in range(8, 11)]
            + [(i, j) for i in range(7, 8) for j in range(8, 11)]
        )
        np.testing.assert_almost_equal(this, ref)
