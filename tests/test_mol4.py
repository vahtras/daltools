import unittest.mock as mock

from daltools import mol
from pytest import approx

from . import tmpdir


class TestMol:

    def setup(self):
        self.tmpdir = tmpdir(__file__)
        dalton_bas = self.tmpdir/"DALTON.BAS"
        self.bas = mol.readin(dalton_bas)
        self.maxDiff = None

    def tearDown(self):
        pass

    def test_pass(self):
        pass

    def test_dist(self):
        assert \
            mol.dist(self.bas[0]["center"][0], self.bas[1]["center"][0]) \
            == approx(2.2852428069)

    def test_num(self):
        assert mol.nuc(self.bas) == approx(56.16649212051651)

    def test_opa(self):
        assert mol.occupied_per_atom(self.bas) == \
            [[0, 1, 2, 3, 4], [0, 1, 2, 3, 4, 5, 6, 7, 8], [0]]

    def test_cpa(self):
        assert mol.contracted_per_atom(self.bas) == [5, 9, 1]

    def test_cpal(self):
        assert mol.contracted_per_atom_l(self.bas) == [[2, 3], [3, 6], [1]]

    def test_print_atoms(self):
        output = """\
Atom type 1 charge 6.000000
center 1 [0.0, -0.5155927197089, -1.4709704593027]
s-functions
    71.616837 [0.15432897, 0.0]
    13.045096 [0.53532814, 0.0]
    3.5305122 [0.44463454, 0.0]
    2.9412494 [0.0, -0.09996723]
    0.6834831 [0.0, 0.39951283]
    0.2222899 [0.0, 0.70011547]
p-functions
    2.9412494 [0.15591627]
    0.6834831 [0.60768372]
    0.2222899 [0.39195739]
Atom type 2 charge 16.000000
center 1 [0.0, 0.0627329390583, 0.7398832489011]
s-functions
    674.4463 [0.154329, 0.0, 0.0]
    122.8512 [0.535328, 0.0, 0.0]
    33.2485 [0.444635, 0.0, 0.0]
    45.16425 [0.0, -0.099967, 0.0]
    10.49518 [0.0, 0.399513, 0.0]
    3.413366 [0.0, 0.700115, 0.0]
    2.621366 [0.0, 0.0, -0.21962]
    0.731354 [0.0, 0.0, 0.225595]
    0.286247 [0.0, 0.0, 0.900398]
p-functions
    45.16425 [0.155916, 0.0]
    10.49518 [0.607684, 0.0]
    3.413366 [0.391957, 0.0]
    2.621366 [0.0, 0.010588]
    0.731354 [0.0, 0.595167]
    0.286247 [0.0, 0.462001]
Atom type 3 charge 2.000000
center 1 [0.0, 0.9194403989461, -2.9770257645433]
s-functions
    6.3624214 [0.15432897]
    1.158923 [0.53532814]
    0.3136498 [0.44463454]
"""

        with mock.patch("daltools.mol.print") as mock_print:
            mol.printbasis(self.bas)
        mock_print.assert_called_once_with(output)
