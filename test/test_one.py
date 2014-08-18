import unittest
import os
import numpy as np
from .. import one
from ..util import blocked

class TestOne(unittest.TestCase):

    def setUp(self):
        n, _ = os.path.splitext(__file__)
        self.tmpdir = n + ".d"
        self.aooneint = os.path.join(self.tmpdir, 'AOONEINT')

    def test_header(self):
        head = one.readhead(self.aooneint)
        self.assertIn("B term (MCD) components of H3+", head["ttitle"])
        self.assertTupleEqual(head["naos"], (2, 1, 0, 0))
        self.assertEqual(head["nsym"], 4)
        self.assertAlmostEqual(head["potnuc"], 1.82903817207)

    def test_isordk(self):
        isordk = one.readisordk(self.aooneint)
        assert isordk["nucdep"] == 3
        #assert isordk["mxcent"] == 120
        self.assertTupleEqual(isordk["chrn"][:3], (1., 1., 1.))
        np.testing.assert_almost_equal(isordk["cooo"][0::120], (0, -0.224906, 0))
        np.testing.assert_almost_equal(isordk["cooo"][1::120], (1, 0.8996236, 0))
        np.testing.assert_almost_equal(isordk["cooo"][2::120], (-1, 0.8996236, 0))

    def test_scfinp(self):
        scfinp = one.readscfinp(self.aooneint)
        self.assertEqual(scfinp["nsym"], 4)
        self.assertTupleEqual(scfinp["cooo"], (0.0, -0.224905893, 0.0, 1.0, 0.899623572, 0.0, -1.0, 0.899623572, 0.0))

    def test_overlap(self):
        Sref = blocked.triangular([2, 1])
        Sref.subblock[0][0, 0] = 1.0
        Sref.subblock[0][1, 0] = 1.24636433    
        Sref.subblock[0][1, 1] = 2.92555541
        Sref.subblock[1][0, 0] = 1.0
        S = one.read("OVERLAP", self.aooneint)
        for s, sref in zip(S.subblock, Sref.subblock):
            np.testing.assert_almost_equal(s, sref)

if __name__ == "__main__":
    unittest.main()
