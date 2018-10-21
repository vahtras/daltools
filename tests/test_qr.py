import unittest
import unittest.mock as mock
import os
import sys
import pdb
import numpy as np
from . import daltools
from daltools.qr import QR, main

class NewTest(unittest.TestCase):

    def setUp(self):
        (n, e,) = os.path.splitext(__file__)
        self.tmpdir = n + '.d'

    def tearDown(self):
        pass

    def assert_(self, this, ref, **kwargs):
        np.testing.assert_allclose(this, ref, **kwargs)

    def test_XXX_A2(self):
        XXX = QR('XDIPLEN', 'XDIPLEN', 'XDIPLEN', tmpdir=self.tmpdir, a2test=True)
        XXXref = 9.652e-05
        self.assert_(XXX, XXXref, rtol=1e-4)



    def test_XXX(self):
        XXX = QR('XDIPLEN', 'XDIPLEN', 'XDIPLEN', tmpdir=self.tmpdir)
        XXXref = 0.00013323
        self.assert_(XXX, XXXref, rtol=1e-5)

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

    def test_mainXXX(self):
        sys.argv[1:] = ['XDIPLEN', 'XDIPLEN', 'XDIPLEN', '-t', self.tmpdir]
        with mock.patch('daltools.qr.print') as mock_print:
            main()
        mock_print.assert_called_once_with("  0.00013323")

if __name__ == "__main__":#pragma no cover
    unittest.main()
