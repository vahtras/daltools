"""Test BASINFO data on Sirius restart files (SIRIUS.RST)"""
import unittest
import sys
import os

try:
    import mock
except ImportError:
    import unittest.mock as mock

import numpy

from daltools.basinfo import BasInfo, main
from daltools import basinfo


class TestBasInfo(unittest.TestCase):
    """Test class for BASINFO data"""

    def setUp(self):
        suppdir = os.path.splitext(__file__)[0] + ".d"
        self.bas_info = BasInfo(os.path.join(suppdir, "SIRIUS.RST"))

    def test_nsym(self):
        """Check number of symmetries NSYM"""
        self.assertEqual(self.bas_info.nsym, 1)

    def test_nbas(self):
        """Check number of basis functions/symmetry NBAS"""
        numpy.testing.assert_equal(self.bas_info.nbas, [5])

    def test_nbast(self):
        """Check total number of basis functions"""
        self.assertEqual(self.bas_info.nbast, 5)

    def test_norb(self):
        """Check number of molecular orbitals/symmetry, NORB"""
        numpy.testing.assert_equal(self.bas_info.norb, [5])

    def test_norbt(self):
        """Check total number of molecular orbitals NORBT"""
        self.assertEqual(self.bas_info.norbt, 5)

    def test_ncmot(self):
        """Check MO-coefficient dimensions"""
        self.assertEqual(self.bas_info.ncmot, 25)

    def test_str(self):
        """Check output text in main"""
        ref = """\
NSYM   :   1
NBAS   :   5
NORB   :   5
NRHF   :   1
IOPRHF :   0
"""
        self.assertEqual(str(self.bas_info), ref)

    @mock.patch.object(basinfo, "BasInfo")
    def test_main_with_arg(self, mock_BasInfo):
        sys.argv[1:] = ["one"]
        main()
        mock_BasInfo.assert_called_once_with("one")

    def test_main_without_arg(self):
        sys.argv[1:] = []
        with self.assertRaises(SystemExit):
            main()


if __name__ == "__main__":
    unittest.main()
