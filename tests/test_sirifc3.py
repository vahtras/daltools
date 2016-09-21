import unittest
import os
import numpy as np
from daltools import sirifc
from util import full

class TestSirIfc(unittest.TestCase):

    def tmpdir(self, name=""):
        n, _ = os.path.splitext(__file__)
        dir_ = n + ".d"
        return os.path.join(dir_, name)
        

    def setUp(self):
        self.ifc = sirifc.sirifc(self.tmpdir('SIRIFC'))

    def test_wrong_file_header(self):
        with self.assertRaises(RuntimeError):
            wrong = sirifc.sirifc(name=self.tmpdir('AOONEINT'))

    def test_potnuc(self):
        self.assertAlmostEqual(self.ifc.potnuc, 31.249215315972)

    def test_emy(self):
        self.assertAlmostEqual(self.ifc.emy, -142.1249229)

    def test_eactive(self):
        self.assertAlmostEqual(self.ifc.eactive, -1.5308779)

    def test_emcscf(self):
        self.assertAlmostEqual(self.ifc.emcscf, -112.4065855)

    def test_istate(self):
        self.assertEqual(self.ifc.istate, 1)

    def test_ispin(self):
        self.assertEqual(self.ifc.ispin, 1)

    def test_nactel(self):
        self.assertEqual(self.ifc.nactel, 2)

    def test_lsym(self):
        self.assertEqual(self.ifc.lsym, 1)

    def test_nisht(self):
        self.assertEqual(self.ifc.nisht, 7)

    def test_nasht(self):
        self.assertEqual(self.ifc.nasht, 2)

    def test_nocct(self):
        self.assertEqual(self.ifc.nocct, 9)

    def test_norbt(self):
        self.assertEqual(self.ifc.norbt, 12)

    def test_nbast(self):
        self.assertEqual(self.ifc.nbast, 12)

    def test_nsym(self):
        self.assertEqual(self.ifc.nsym, 1)

    def test_pv(self):
        pv = self.ifc.pv
        np.testing.assert_allclose(
            pv.diagonal(), [1.87279201, -0.24404616, 0.12720799],
            rtol=1e-7
            )
