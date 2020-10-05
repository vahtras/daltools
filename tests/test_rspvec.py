import numpy as np
from pytest import approx, raises

from daltools import rspvec, sirifc
from . import tmpdir


class TestRspVec:

    def setup(self):
        self.tmpdir = tmpdir(__file__)
        self.RSPVEC = self.tmpdir/"RSPVEC"
        self.SIRIFC = self.tmpdir/"SIRIFC"
        self.E3VEC = self.tmpdir/"E3VEC"
        self.ifc = sirifc.sirifc(name=self.SIRIFC)

    def test_read(self):
        Nx = rspvec.read("XDIPLEN", propfile=self.RSPVEC)["XDIPLEN"]
        this = Nx[12]
        ref = -0.75732690
        assert this == approx(ref)

    def test_read_w(self):
        Nx = rspvec.read("XDIPLEN", freqs=(0.5,), propfile=self.RSPVEC)[
            ("XDIPLEN", 0.5)
        ]
        this = Nx[12]
        ref = -2.242435
        assert this == approx(ref)

    def test_read_e3(self):
        Nx = rspvec.read("XDIPLEN XDIPLEN", propfile=self.E3VEC)["XDIPLEN XDIPLEN"]
        this = Nx[1]
        ref = 0.03874785
        assert this == approx(ref)

    def test_read_e3_w(self):
        Nx = rspvec.read("XDIPLEN XDIPLEN", freqs=(0.5,), propfile=self.E3VEC)[
            ("XDIPLEN XDIPLEN", 0.5)
        ]
        this = Nx[1]
        ref = 0.05668560
        assert this == approx(ref)

    def test_read_missing(self):
        with raises(rspvec.RspVecError):
            rspvec.read("WRONGLAB", propfile=self.RSPVEC)

    def test_read_all(self):
        N1, N2 = rspvec.readall("XDIPLEN", self.RSPVEC)[:4:2]
        f1, f2 = N1[1], N2[1]
        assert [f1, f2] == [0, 0.5]
        this = (N1[0][12], N2[0][12])
        ref = (-0.75732690, -2.242435)
        np.testing.assert_almost_equal(this, ref)

    def test_tomat(self):
        Nx = rspvec.read("XDIPLEN", propfile=self.RSPVEC)["XDIPLEN"]
        kx = rspvec.tomat(Nx, self.ifc, tmpdir=self.tmpdir)
        this = kx[8, 3]
        ref = 0.75732690
        assert this == approx(ref)

    def test_tovec(self):
        Nx = rspvec.read("XDIPLEN", propfile=self.RSPVEC)["XDIPLEN"]
        kx = rspvec.tomat(Nx, self.ifc, tmpdir=self.tmpdir)
        Nx = rspvec.tovec(kx, self.ifc, tmpdir=self.tmpdir)
        this = Nx[44]
        ref = 0.75732690
        assert this == approx(ref)

    def test_jwop(self):
        this = list(rspvec.jwop(self.ifc))
        ref = [(i, j) for i in range(8) for j in range(8, 12)]
        assert  this == ref