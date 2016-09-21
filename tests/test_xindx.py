import unittest
import os
import sys
import numpy
from . import daltools
from daltools.sirrst import SiriusRestart, main
from daltools.sirifc import sirifc
from util.blocked import BlockDiagonalMatrix

class TestSirRst(unittest.TestCase):

    def setUp(self):
        self.suppdir = os.path.splitext(__file__)[0] + ".d"
        self.daltgz = os.path.join(self.suppdir, 'mc_h2.tar.gz')
        self.sirrst = SiriusRestart(tgz=self.daltgz)
        self.sirifc = sirifc(os.path.join(self.suppdir, 'SIRIFC'))

    def test_cmo(self):
        numpy.testing.assert_almost_equal(
            self.sirrst.cmo[0],
            [[0.6212021, -.8425698], [0.6212021, 0.8425698]]
            )

    def test_xindx(self):
        numpy.testing.assert_almost_equal(
            self.sirrst.ci,
            [0.9491331152, 0, 0, -0.314874991]
            )

    def test_dets1(self):
        self.assertTupleEqual(
            list(self.sirifc.xindx())[0],
            ((0,), (0,))
            )

    def test_dets2(self):
        self.assertTupleEqual(
            list(self.sirifc.xindx())[3],
            ((1,), (1,))
            )

if __name__ == "__main__":
    unittest.main()
