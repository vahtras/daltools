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
        self.sirifc3 = sirifc(os.path.join(self.suppdir, 'SIRIFC3'))
        self.sirifc4 = sirifc(os.path.join(self.suppdir, 'alt4/SIRIFC'))
        self.xindx4 = list(self.sirifc4.xindx())

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

    def test_dets3(self):
        self.assertTupleEqual(
            list(self.sirifc3.xindx())[0],
            ((0, 1,),())
            )

    def test_dets4(self):
        self.assertListEqual(self.xindx4, [
            ((0, 1), (0, 1)),
            ((0, 2), (0, 1)),
            ((1, 2), (0, 1)),
            ((0, 3), (0, 1)),
            ((1, 3), (0, 1)),
            ((2, 3), (0, 1)),
            ((0, 1), (0, 2)),
            ((0, 2), (0, 2)),
            ((1, 2), (0, 2)),
            ((0, 3), (0, 2)),
            ((1, 3), (0, 2)),
            ((2, 3), (0, 2)),
            ((0, 1), (1, 2)),
            ((0, 2), (1, 2)),
            ((1, 2), (1, 2)),
            ((0, 3), (1, 2)),
            ((1, 3), (1, 2)),
            ((2, 3), (1, 2)),
            ((0, 1), (0, 3)),
            ((0, 2), (0, 3)),
            ((1, 2), (0, 3)),
            ((0, 3), (0, 3)),
            ((1, 3), (0, 3)),
            ((2, 3), (0, 3)),
            ((0, 1), (1, 3)),
            ((0, 2), (1, 3)),
            ((1, 2), (1, 3)),
            ((0, 3), (1, 3)),
            ((1, 3), (1, 3)),
            ((2, 3), (1, 3)) ,
            ((0, 1), (2, 3)),
            ((0, 2), (2, 3)),
            ((1, 2), (2, 3)),
            ((0, 3), (2, 3)),
            ((1, 3), (2, 3)),
            ((2, 3), (2, 3)) 
            ])


if __name__ == "__main__":
    unittest.main()
