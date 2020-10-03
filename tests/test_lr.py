import unittest.mock as mock

from daltools.lr import LR, main
from pytest import approx

from .common_tests import tmpdir


class TestLR:

    def setup(self):
        self.tmpdir = tmpdir(__file__)

    def test_XX(self):
        XX = LR("XDIPLEN", "XDIPLEN", 0, self.tmpdir)
        XXref = -2.461169664950
        assert XX == approx(XXref)

    def test_YY(self):
        YY = LR("YDIPLEN", "YDIPLEN", 0, self.tmpdir)
        YYref = -6.184500121159
        assert YY == approx(YYref)

    def test_YZ(self):
        YZ = LR("YDIPLEN", "ZDIPLEN", 0, self.tmpdir)
        YZref = 0.016158638023
        assert YZ == approx(YZref)

    def test_ZZ(self):
        ZZ = LR("ZDIPLEN", "ZDIPLEN", 0, self.tmpdir)
        ZZref = -10.308181624834
        assert ZZ == approx(ZZref)

    def test_main_ZZ(self):
        import sys

        sys.argv[1:] = ["XDIPLEN", "XDIPLEN", "-t", str(self.tmpdir)]
        with mock.patch("daltools.lr.print") as mock_print:
            main()
            mock_print.assert_called_once_with("-2.461170")
