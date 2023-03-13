import unittest.mock as mock

import numpy as np
from pytest import approx, mark

from daltools import mol

from . import tmpdir


def _test_data(label, charge, prim, cont):
    return {
        'charge': charge,
        'label': [label],
        'center': [[0, 0, 0]],
        'basis': [
            [
                {
                    'prim': np.random.random(p).tolist(),
                    'cont': np.random.random((p, c)).tolist(),
                }
            ]
            for p, c in zip(prim, cont)
        ],
    }


class TestMol:

    def setup_method(self):
        self.tmpdir = tmpdir(__file__)
        dalton_bas = self.tmpdir/"DALTON.BAS"
        self.bas = mol.readin(dalton_bas)
        self.maxDiff = None

    def tearDown(self):
        pass

    def test_pass(self):
        pass

    def test_dist(self):
        assert mol.dist(
            self.bas[0]["center"][0], self.bas[1]["center"][0]
            ) == approx(2.2852428069)

    def test_num(self):
        assert mol.nuc(self.bas) == approx(31.2492153162)

    def test_opa(self):
        assert mol.occupied_per_atom(self.bas) == \
            [[0, 1, 2, 3, 4], [0, 1, 2, 3, 4], [0], [0]]

    @mark.parametrize(
        'molecule, expected',
        [
            ([_test_data('H', 1, [1], [1])], [[0]]),
            ([_test_data('He', 2, [2], [2])], [[0]]),
            ([_test_data('Li', 3, [2], [2])], [[0, 1]]),
            ([_test_data('B', 4, [2, 1], [2, 1])], [[0, 1]]),
            ([_test_data('C', 6, [3, 2], [3, 2])], [[0, 1, 3, 4, 5]]),
            (
                [_test_data('Na', 11, [8, 4, 2, 1], [4, 3, 2, 1])],
                [[0, 1, 2, 4, 5, 6]]
            ),
            (
                [_test_data('Cl', 17, [8, 4, 2, 1], [4, 3, 2, 1])],
                [[0, 1, 2, 4, 5, 6, 7, 8, 9]]
            ),
            (
                [_test_data('Ca', 20, [8, 4, 2, 1], [4, 3, 2, 1])],
                [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]
            ),
            (
                [_test_data('Mn', 25, [12, 10, 8, 4], [5, 4, 3, 2])],
                [[
                    0, 1, 2, 3,
                    5, 6, 7, 8, 9, 10,
                    17, 18, 19, 20, 21
                ]]
            ),
            (
                [_test_data('Br', 35, [12, 10, 8, 4], [5, 4, 3, 2])],
                [[
                    0, 1, 2, 3,
                    5, 6, 7, 8, 9, 10, 11, 12, 13,
                    17, 18, 19, 20, 21
                ]]
            ),
        ],
        ids=['H', 'He', 'Li', 'B', 'C', 'Na', 'Cl', 'Ca', 'Mn', 'Br']
    )
    def test_atomic_opas(self, molecule, expected):
        assert mol.occupied_per_atom(molecule) == expected

    def test_cpa(self):
        assert mol.contracted_per_atom(self.bas) == [5, 5, 1, 1]

    def test_cpal(self):
        assert mol.contracted_per_atom_l(self.bas) == \
            [[2, 3], [2, 3], [1], [1]]

    def test_print_atoms(self):
        output = """\
Atom type 1 charge 6.000000
center 1 [1.74063211e-05, 0.0010502765856, -1.1458244562083]
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
Atom type 2 charge 8.000000
center 1 [0.0, 0.0010582718101, 1.1394183506149]
s-functions
    130.70932 [0.15432897, 0.0]
    23.808861 [0.53532814, 0.0]
    6.4436083 [0.44463454, 0.0]
    5.0331513 [0.0, -0.09996723]
    1.1695961 [0.0, 0.39951283]
    0.380389 [0.0, 0.70011547]
p-functions
    5.0331513 [0.15591627]
    1.1695961 [0.60768372]
    0.380389 [0.39195739]
Atom type 3 charge 1.000000
center 1 [-7.11056897e-05, 1.7705033753955, -2.2396975555016]
center 2 [-6.98661118e-05, -1.7998043816988, -2.2005635940297]
s-functions
    3.4252509 [0.15432897]
    0.6239137 [0.53532814]
    0.1688554 [0.44463454]
"""
        with mock.patch("daltools.mol.print") as mock_print:
            mol.printbasis(self.bas)
        mock_print.assert_called_once_with(output)
