import unittest
import os
import numpy as np
from .. import sirifc
from ..util import full

class TestSirIfc(unittest.TestCase):

    def setUp(self):
        n, _ = os.path.splitext(__file__)
        tmpdir = n + ".d"
        self.ifc = sirifc.sirifc(os.path.join(tmpdir, 'SIRIFC'))

    def test_potnuc(self):
        self.assertAlmostEqual(self.ifc.potnuc, 1.8290382)

    def test_emcscf(self):
        self.assertAlmostEqual(self.ifc.emcscf, -1.2492926)

    def test_dimensions(self):
        this = (self.ifc.nisht, self.ifc.nasht, self.ifc.norbt, self.ifc.nbast, self.ifc.nsym)
        self.assertTupleEqual(this, (0, 3, 3, 3, 4))

    def test_cmo(self):
        cmo1 = np.array([[0.45236534, 1.38834090], [0.36316622, -0.77259529]])
        np.testing.assert_almost_equal(self.ifc.cmo[0], cmo1)

    def test_dvself(self):
        dv = np.array([1.97086440, 0, 0.00940855, 0, 0, 0.01972706])
        np.testing.assert_almost_equal(self.ifc.dv.view(full.matrix), dv)

if __name__ == "__main__":
    unittest.main()
