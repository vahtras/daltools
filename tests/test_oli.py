import os
from util.full import unit, init, matrix
from daltools import oli
from .common_tests import assert_


def setup():
    global suppdir
    n, e = os.path.splitext(__file__)
    suppdir = n + ".d"


def test_e2n_S():
    refe2 = init([[0.78643356, -0.50624296], [-0.50624296, 0.78643356]])
    e2 = [oli.e2n(n, tmpdir=suppdir) for n in unit(2)]
    assert_(e2, refe2)


def test_s2n_S():
    refs2 = matrix.diag([2.0, -2.0])
    s2 = [oli.s2n(n, tmpdir=suppdir) for n in unit(2)]
    assert_(s2, refs2)


if __name__ == "__main__":
    setup()
    test_e2n_S()
