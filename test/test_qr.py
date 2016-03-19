import unittest
import os
import pdb
import numpy as np
from ..qr import QR

class NewTest(unittest.TestCase):

    def setUp(self):
        (n, e,) = os.path.splitext(__file__)
        self.tmpdir = n + '.d'

    def tearDown(self):
        pass

    def assert_(self, this, ref):
        assert np.allclose(this, ref)

    def test_XXX_A2(self):
        XXX = QR('XDIPLEN', 'XDIPLEN', 'XDIPLEN', tmpdir=self.tmpdir, a2test=True)
        XXXref = 9.652e-05
        self.assert_(XXX, XXXref)



    def test_XXX(self):
        XXX = QR('XDIPLEN', 'XDIPLEN', 'XDIPLEN', tmpdir=self.tmpdir)
        XXXref = 0.00013323
        self.assert_(XXX, XXXref)

    def test_ZXX(self):
        ZXX = QR('ZDIPLEN', 'XDIPLEN', 'XDIPLEN', tmpdir=self.tmpdir)
        ZXXref = -1.68075251
        self.assert_(ZXX, ZXXref)

    def test_XZX(self):
        XZX = QR('XDIPLEN', 'ZDIPLEN', 'XDIPLEN', tmpdir=self.tmpdir)
        XZXref = -1.68075251
        self.assert_(XZX, XZXref)

    def test_XXZ(self):
        XXZ = QR('XDIPLEN', 'XDIPLEN', 'ZDIPLEN', tmpdir=self.tmpdir)
        XXZref = -1.68075251
        self.assert_(XXZ, XXZref)


if __name__ == "__main__":#pragma no cover
    unittest.main()
