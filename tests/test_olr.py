import numpy as np
from pytest import approx, raises

from daltools.lr import LR
from . import tmpdir


class TestOpenLR:

    def setup(self):
        self.tmpdir = tmpdir(__file__)

    def test_XX(self):
        XX = LR("XDIPLEN", "XDIPLEN", 0, self.tmpdir)
        XXref = -5.455606903637
        assert XX == approx(XXref)

    def test_YY(self):
        YY = LR("YDIPLEN", "YDIPLEN", 0, self.tmpdir)
        YYref = -10.31180304740
        assert  YY == approx(YYref)

    def test_YZ(self):
        YZ = LR("YDIPLEN", "ZDIPLEN", 0, self.tmpdir)
        YZref = 0.2324294799056
        assert YZ == approx(YZref)

    def test_ZZ(self):
        ZZ = LR("ZDIPLEN", "ZDIPLEN", 0, self.tmpdir)
        ZZref = -3.124432068117
        assert  ZZ == approx(ZZref)