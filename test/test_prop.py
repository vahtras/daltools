import unittest
import os
import numpy as np
from .. import prop, sirifc

class TestProp(unittest.TestCase):

    def setUp(self):
        n, e = os.path.splitext(__file__)
        tmpdir = n + ".d"
        self.propfile = os.path.join(tmpdir, 'AOPROPER')
        self.ifcfile = os.path.join(tmpdir, 'SIRIFC')

    def test_xdiplen(self):
        x, = prop.read('XDIPLEN', filename=self.propfile, unpack=False)
        xref = [0, 0, 0, 0.62318216, 2.00000000, 0.00000000]
        np.testing.assert_almost_equal(x, xref)

    def test_ydiplen(self):
        y, = prop.read('YDIPLEN', filename=self.propfile, unpack=False)
        yref = [-0.22490589, 0.42047202, 2.63189861, 0.00000000, 0.00000000, 0.96659568]
        np.testing.assert_almost_equal(y, yref)

    def test_zdiplen(self):
        z, = prop.read('ZDIPLEN', filename=self.propfile, unpack=False)
        zref = [0, 0, 0, 0, 0, 0]
        np.testing.assert_almost_equal(z, zref)

if __name__ == "__main__":
    unittest.main()

