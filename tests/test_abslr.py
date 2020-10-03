import numpy.testing as npt
from pytest import approx

from daltools.rspvec import read
from daltools.lr import LR

from .common_tests import tmpdir


class TestALR:
    def setup(self):
        self.tmpdir = tmpdir(__file__)

    def test_first_xvec(self):
        vecs = read("XDIPLEN", freqs=(0.4425,), propfile=self.tmpdir/"ABSVECS")
        npt.assert_allclose(
            vecs[("XDIPLEN", 0.4425, 0.0)][:10],
            [
                -1.7563977677089527e-09,
                0.0015755554234863397,
                -5.504742822192291e-10,
                0.002779217059896619,
                -0.0034886672773900523,
                7.051844231481651e-10,
                -0.005062366558336546,
                -0.006285259584578197,
                -0.008674385885382397,
                -7.986196390011771e-11,
            ],
        )

    def test_first_yvec(self):
        vecs = read(
            "XDIPLEN", "YDIPLEN", freqs=(0.4425,), propfile=self.tmpdir / "ABSVECS"
        )
        npt.assert_allclose(
            vecs[("YDIPLEN", 0.4425, 0.0)][:10],
            [
                -3.062472422126098e-08,
                0.0027792414636958196,
                -1.1898942484099601e-07,
                -0.0015755667065659468,
                -0.005062426714257009,
                1.0479887117851179e-07,
                0.00348869833460343,
                -0.008674366708527649,
                0.006285236842617945,
                2.145201913436897e-07,
            ],
        )

    def test_number_vecs(self):
        vecs = read(
            "XDIPLEN",
            "YDIPLEN",
            "ZDIPLEN",
            freqs=(0.4425, 0.49),
            propfile=self.tmpdir/"ABSVECS",
            lr_vecs=True,
        )
        assert len(vecs) == 6

    def test_real_alpha_xx(self):
        w = 0.4425
        re_alpha, _ = LR(
            "XDIPLEN", "XDIPLEN", w, tmpdir=self.tmpdir, absorption=True
        )
        assert -re_alpha == approx(30.854533)

    def test_real_alpha_yy(self):
        w = 0.4450
        re_alpha, _ = LR(
            "YDIPLEN", "YDIPLEN", w, tmpdir=self.tmpdir, absorption=True
        )
        assert -re_alpha == approx(31.544221)

    def test_real_alpha_zz(self):
        w = 0.4475
        re_alpha, _ = LR(
            "ZDIPLEN", "ZDIPLEN", w, tmpdir=self.tmpdir, absorption=True
        )
        assert -re_alpha == approx(32.275296)

    def test_real_alpha_xy(self):
        w = 0.4425
        re_alpha, _ = LR(
            "XDIPLEN", "YDIPLEN", w, tmpdir=self.tmpdir, absorption=True
        )
        assert -round(re_alpha, 6) == approx(-0.000004)

    def test_real_alpha_xz(self):
        w = 0.4425
        re_alpha, _ = LR(
            "XDIPLEN", "ZDIPLEN", w, tmpdir=self.tmpdir, absorption=True
        )
        assert -re_alpha == approx(0.000000)

    def test_real_alpha_yz(self):
        w = 0.4425
        re_alpha, _ = LR(
            "XDIPLEN", "ZDIPLEN", w, tmpdir=self.tmpdir, absorption=True
        )
        assert -re_alpha == approx(0.000000)

    def test_real_im(self):
        w = 0.4425
        _, im_alpha = LR(
            "XDIPLEN", "XDIPLEN", w, tmpdir=self.tmpdir, absorption=True
        )
        assert -im_alpha == approx(1.228334)
