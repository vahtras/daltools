import os
import unittest

import numpy as np
from pytest import approx

from daltools import one

from . import tmpdir

class TestOne:

    def setup(self):
        self.tmpdir = tmpdir(__file__)
        self.aooneint = self.tmpdir/"AOONEINT"
        self.header = one.readhead(self.aooneint)

    def test_header_title(self):
        assert "CH2O" in  self.header["ttitle"]

    def test_header_naos(self):
        assert self.header["naos"] == (7, 2, 3, 0)

    def test_header_nsym(self):
        assert self.header["nsym"] == 4

    def test_header_potnuc(self):
        assert self.header["potnuc"] == approx(33.941796143849)

    def test_isordk_nucdep(self):
        isordk = one.readisordk(self.aooneint)
        assert  isordk["nucdep"] == 4

    def test_isordk_chrn(self):
        isordk = one.readisordk(self.aooneint)
        assert isordk["chrn"][:3] == (6.0, 8.0, 1.0)

    # @unittest.skip('Wrong format here')
    def test_isordk_cooo(self):
        isordk = one.readisordk(self.aooneint)
        C = [0.0000000000, 0.0000000000, -0.8802560791]
        O = [0.0000000000, 0.0000000000, 1.0094700459]
        H1 = [0.0000000000, 1.8897261250, -2.7699822041]
        H2 = [0.0000000000, -1.8897261250, -2.7699822041]
        np.testing.assert_almost_equal(isordk["cooo"][0::500], C)
        np.testing.assert_almost_equal(isordk["cooo"][1::500], O)
        np.testing.assert_almost_equal(isordk["cooo"][2::500], H1)
        np.testing.assert_almost_equal(isordk["cooo"][3::500], H2)

    def test_scfinp(self):
        scfinp = one.readscfinp(self.aooneint)
        assert  scfinp["nsym"] == 4
        coor_bohr = (
            0.0000000000,
            0.0000000000,
            -0.8802560791,
            0.0000000000,
            0.0000000000,
            1.0094700459,
            0.0000000000,
            1.8897261250,
            -2.7699822041,
            0.0000000000,
            -1.8897261250,
            -2.7699822041,
        )

        np.testing.assert_almost_equal(scfinp["cooo"], coor_bohr)

    def test_overlap(self):
        Sref = [
            [
                1.00000000,
                0.24836239,
                1.00000000,
                0.00000000,
                0.00000000,
                1.00000000,
                0.00005599,
                0.05942391,
                0.10228506,
                1.00000000,
                0.07053810,
                0.49046481,
                0.53543774,
                0.23670394,
                1.00000000,
                -0.11599611,
                -0.37963762,
                -0.27726502,
                -0.00000000,
                0.00000000,
                1.00000000,
                0.05926125,
                0.64602901,
                -0.48763174,
                0.00548289,
                0.09221295,
                -0.09620875,
                2.23843306,
            ],
            [1.00000000, 0.31956952, 1.00000000],
            [1.00000000, 0.31956952, 1.00000000, 0.48763174, 0.04810438, 1.76156694],
        ]

        S = one.read("OVERLAP", self.aooneint)
        np.testing.assert_almost_equal(np.array(S.subblock[0]), Sref[0])
        np.testing.assert_almost_equal(np.array(S.subblock[1]), Sref[1])
        np.testing.assert_almost_equal(np.array(S.subblock[2]), Sref[2])


if __name__ == "__main__":  # pragma no cover
    unittest.main()
