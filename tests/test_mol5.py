import unittest
import os
import sys
from . import daltools
from daltools import mol


class MolTest(unittest.TestCase):

    def setUp(self):
        pwd = os.path.dirname(__file__)
        dalton_bas = os.path.join(pwd, 'test_mol5.d', 'DALTON.BAS')
        self.bas = mol.readin(dalton_bas)
        self.maxDiff = None

    def tearDown(self):
        pass

    def test_pass(self):
        with self.assertRaises(NotImplementedError):
            mol.occupied_per_atom(self.bas)
