import unittest
import unittest.mock as mock
import sys

import numpy as np
from pytest import approx, raises

from daltools import one
from . import tmpdir


class TestOne:
    def setup_method(self):
        self.tmpdir = tmpdir(__file__)
        self.aooneint = self.tmpdir/"AOONEINT"
        self.aoproper = self.tmpdir/"AOPROPER"
        self.header = one.readhead(self.aooneint)
        self.maxDiff = None

    def test_header_title(self):
        assert "CH2O" in self.header["ttitle"]

    def test_header_naos(self):
        assert self.header["naos"] == (12,)

    def test_header_nsym(self):
        assert self.header["nsym"] == 1

    def test_header_potnuc(self):
        assert self.header["potnuc"] == approx(31.249215316217)

    def test_isordk_nucdep(self):
        isordk = one.readisordk(self.aooneint)
        assert isordk["nucdep"] == 4

    def test_isordk_chrn(self):
        isordk = one.readisordk(self.aooneint)
        assert isordk["chrn"][:3] == (6.0, 8.0, 1.0)

    def test_isordk_cooo(self):
        isordk = one.readisordk(self.aooneint)
        C = [-3.0015786160, -1.4563174382, 0.0550080378]
        O = [-3.1314330364, 0.8240509816, -0.0184248297]
        H1 = [-1.1728925345, -2.4468589606, 0.1025195320]
        H2 = [-4.7395143797, -2.6116033945, 0.0761219478]
        np.testing.assert_almost_equal(isordk["cooo"][0::120], C)
        np.testing.assert_almost_equal(isordk["cooo"][1::120], O)
        np.testing.assert_almost_equal(isordk["cooo"][2::120], H1)
        np.testing.assert_almost_equal(isordk["cooo"][3::120], H2)

    def test_scfinp(self):
        scfinp = one.readscfinp(self.aooneint)
        assert scfinp["nsym"] == 1
        coor_angstrom = (
            -1.588367,
            -0.770650,
            0.029109,
            -1.657083,
            0.436069,
            -0.009750,
            -0.620668,
            -1.294822,
            0.054251,
            -2.508043,
            -1.382001,
            0.040282,
        )
        coor_bohr = [i / 0.52917721 for i in coor_angstrom]
        np.testing.assert_almost_equal(scfinp["cooo"], coor_bohr)

    def test_overlap(self):
        Sref = [
            1.00000000,
            0.24836239,
            1.00000000,
            0.00000000,
            0.00000000,
            1.00000000,
            0.00000000,
            0.00000000,
            0.00000000,
            1.00000000,
            0.00000000,
            0.00000000,
            0.00000000,
            0.00000000,
            1.00000000,
            0.00000126,
            0.03708896,
            -0.00354693,
            0.06228751,
            -0.00200579,
            1.00000000,
            0.03664911,
            0.36526353,
            -0.02523128,
            0.44308559,
            -0.01426833,
            0.23670394,
            1.00000000,
            0.00349314,
            0.01832934,
            0.21019458,
            0.02979663,
            -0.00095952,
            0.00000000,
            0.00000000,
            1.00000000,
            -0.06134287,
            -0.32188081,
            0.02979663,
            -0.31136609,
            0.01685004,
            0.00000000,
            0.00000000,
            0.00000000,
            1.00000000,
            0.00197538,
            0.01036527,
            -0.00095952,
            0.01685004,
            0.21134872,
            0.00000000,
            0.00000000,
            0.00000000,
            0.00000000,
            1.00000000,
            0.06072046,
            0.48453953,
            0.40747211,
            -0.22071478,
            0.01058662,
            0.00476429,
            0.07308063,
            0.04174833,
            -0.06972286,
            0.00257806,
            1.00000000,
            0.06021809,
            0.48250496,
            -0.38496913,
            -0.25590672,
            0.00467693,
            0.00488694,
            0.07467580,
            -0.03512957,
            -0.07505408,
            0.00206544,
            0.14255017,
            1.00000000,
        ]

        S = one.read("OVERLAP", self.aooneint)
        np.testing.assert_almost_equal(np.array(S.subblock[0]), Sref)

    def test_main(self):
        sys.argv[1:] = [str(self.aooneint)]
        with mock.patch("daltools.one.print") as mock_print:
            one.main()
        mock_print.assert_not_called()

    def test_main_head(self):
        sys.argv[1:] = [str(self.aooneint), "--head"]
        ref_output = """\
Header on AOONEINT
ttitle CH2O                                                                    ------------------------                                                
nsym 1
naos (12,)
potnuc   31.24922
int_fmt q
float_fmt d"""
        with mock.patch("daltools.one.print") as mock_print:
            one.main()
        # mock_print.assert_called_once_with(ref_output)
        calls = [mock.call(s) for s in ref_output.split("\n")]
        mock_print.assert_has_calls(calls)

    def test_main_isordk(self):
        sys.argv[1:] = [str(self.aooneint), "--isordk"]
        ref_output = """\
nucdep=4 mxcent=120

 (4,) 
              Column   1
       1      6.00000000
       2      8.00000000
       3      1.00000000
       4      1.00000000


 (3, 4) 
              Column   1    Column   2    Column   3    Column   4
       1     -3.00157862   -3.13143304   -1.17289253   -4.73951438
       2     -1.45631744    0.82405098   -2.44685896   -2.61160339
       3      0.05500804   -0.01842483    0.10251953    0.07612195"""

        with mock.patch("daltools.one.print") as mock_print:
            one.main()
        calls = [mock.call(s) for s in ref_output.split("\n")]
        mock_print.assert_has_calls([])

    def test_main_scfinp(self):
        sys.argv[1:] = [str(self.aooneint), "--scfinp"]
        ref_output = """\
ttitle CH2O                                                                    ------------------------                                                
nsym 1
naos (12,)
potnuc 31.249215
kmax 8
ncent (1, 1, 1, 2, 2, 2, 3, 4)
nbasis 12
jtran (1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
itran (1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 0, 0, 11, 0, 0, 0, 0, 0, 0, 0, 12, 0, 0, 0, 0, 0, 0, 0)
ctran (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
inamn (2314885530818453571, 2314885530818453571, 2314885530818453571, 2314885530818453571, 2314885530818453571, 2314885530818453583, 2314885530818453583, 2314885530818453583, 2314885530818453583, 2314885530818453583, 2314885530818453576, 2314885530818453576)
iptyp (1, 1, 2, 3, 4, 1, 1, 2, 3, 4, 1, 1)
dpnuc (-48.973342901183, -7.203959131642938, 0.3612910686591494)
nucdep 4
cooo (-3.001578615977693, -1.4563174382263098, 0.0550080377724384, -3.1314330363827527, 0.8240509815998296, -0.0184248297186875, -1.1728925345475214, -2.4468589606164497, 0.10251953200702724, -4.739514379707297, -2.611603394467266, 0.07612194776699177)
ifxyz (0, 0, 0)
dummy 1e+20
qpol (1e+20, 1e+20, 1e+20, 1e+20, 1e+20, 1e+20)
qq (1e+20, 1e+20, 1e+20)
jfxyz (-9999999, -9999999, -9999999)"""

        with mock.patch("daltools.one.print") as mock_print:
            one.main()
        calls = [mock.call(s) for s in ref_output.split("\n")]
        mock_print.assert_has_calls(calls)

    def test_main_label(self):
        sys.argv[1:] = [str(self.aooneint), "--label", "OVERLAP", "-v"]
        ref_output = """\
OVERLAP 
Block 1

    1.00000000
    0.24836239    1.00000000
    0.00000000    0.00000000    1.00000000
    0.00000000    0.00000000    0.00000000    1.00000000
    0.00000000    0.00000000    0.00000000    0.00000000    1.00000000
    0.00000126    0.03708896   -0.00354693    0.06228751   -0.00200579    1.00000000
    0.03664911    0.36526353   -0.02523128    0.44308559   -0.01426833    0.23670394    1.00000000
    0.00349314    0.01832934    0.21019458    0.02979663   -0.00095952    0.00000000    0.00000000    1.00000000
   -0.06134287   -0.32188081    0.02979663   -0.31136609    0.01685004    0.00000000    0.00000000    0.00000000    1.00000000
    0.00197538    0.01036527   -0.00095952    0.01685004    0.21134872    0.00000000    0.00000000    0.00000000    0.00000000    1.00000000
    0.06072046    0.48453953    0.40747211   -0.22071478    0.01058662    0.00476429    0.07308063    0.04174833   -0.06972286    0.00257806    1.00000000
    0.06021809    0.48250496   -0.38496913   -0.25590672    0.00467693    0.00488694    0.07467580   -0.03512957   -0.07505408    0.00206544    0.14255017    1.00000000
"""

        with mock.patch("daltools.one.print") as mock_print:
            one.main()
        mock_print.assert_called_once_with(ref_output)

    def test_main_label_unpack(self):
        sys.argv[1:] = [str(self.aooneint), "--label", "OVERLAP", "-v", "-u"]
        ref_output = """\
OVERLAP 
 (12, 12)
              Column   1    Column   2    Column   3    Column   4    Column   5
       1      1.00000000    0.24836239    0.00000000    0.00000000    0.00000000
       2      0.24836239    1.00000000    0.00000000    0.00000000    0.00000000
       3      0.00000000    0.00000000    1.00000000    0.00000000    0.00000000
       4      0.00000000    0.00000000    0.00000000    1.00000000    0.00000000
       5      0.00000000    0.00000000    0.00000000    0.00000000    1.00000000
       6      0.00000126    0.03708896   -0.00354693    0.06228751   -0.00200579
       7      0.03664911    0.36526353   -0.02523128    0.44308559   -0.01426833
       8      0.00349314    0.01832934    0.21019458    0.02979663   -0.00095952
       9     -0.06134287   -0.32188081    0.02979663   -0.31136609    0.01685004
      10      0.00197538    0.01036527   -0.00095952    0.01685004    0.21134872
      11      0.06072046    0.48453953    0.40747211   -0.22071478    0.01058662
      12      0.06021809    0.48250496   -0.38496913   -0.25590672    0.00467693

              Column   6    Column   7    Column   8    Column   9    Column  10
       1      0.00000126    0.03664911    0.00349314   -0.06134287    0.00197538
       2      0.03708896    0.36526353    0.01832934   -0.32188081    0.01036527
       3     -0.00354693   -0.02523128    0.21019458    0.02979663   -0.00095952
       4      0.06228751    0.44308559    0.02979663   -0.31136609    0.01685004
       5     -0.00200579   -0.01426833   -0.00095952    0.01685004    0.21134872
       6      1.00000000    0.23670394    0.00000000    0.00000000    0.00000000
       7      0.23670394    1.00000000    0.00000000    0.00000000    0.00000000
       8      0.00000000    0.00000000    1.00000000    0.00000000    0.00000000
       9      0.00000000    0.00000000    0.00000000    1.00000000    0.00000000
      10      0.00000000    0.00000000    0.00000000    0.00000000    1.00000000
      11      0.00476429    0.07308063    0.04174833   -0.06972286    0.00257806
      12      0.00488694    0.07467580   -0.03512957   -0.07505408    0.00206544

              Column  11    Column  12
       1      0.06072046    0.06021809
       2      0.48453953    0.48250496
       3      0.40747211   -0.38496913
       4     -0.22071478   -0.25590672
       5      0.01058662    0.00467693
       6      0.00476429    0.00488694
       7      0.07308063    0.07467580
       8      0.04174833   -0.03512957
       9     -0.06972286   -0.07505408
      10      0.00257806    0.00206544
      11      1.00000000    0.14255017
      12      0.14255017    1.00000000
"""

        with mock.patch("daltools.one.print") as mock_print:
            one.main()
        mock_print.assert_called_once_with(ref_output)

    def test_read_wrong_file(self):
        with raises(RuntimeError):
            one.readhead(self.aoproper)

    def test_wrong_integer_format(self):
        class Dummy(object):
            reclen = 7

            def __len__(self):
                return self.reclen

        with raises(RuntimeError):
            i = one._get_integer_format(Dummy())
