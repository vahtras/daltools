import numpy as np
from pytest import approx, raises

from daltools import sirifc
from util import full

from . import tmpdir


class TestSirIfc:

    def setup_method(self):
        self.tmpdir = tmpdir(__file__)
        self.ifc = sirifc.sirifc(self.tmpdir/"SIRIFC")

    def test_wrong_file_header(self):
        with raises(RuntimeError):
            sirifc.sirifc(name=self.tmpdir/"AOONEINT")

    def test_potnuc(self):
        assert self.ifc.potnuc == approx(31.249215315972)

    def test_emy(self):
        assert self.ifc.emy == approx(-142.1249229)

    def test_eactive(self):
        assert self.ifc.eactive == approx(-1.5308779)

    def test_emcscf(self):
        assert self.ifc.emcscf == approx(-112.4065855)

    def test_istate(self):
        assert self.ifc.istate == 1

    def test_ispin(self):
        assert self.ifc.ispin == 1

    def test_nactel(self):
        assert self.ifc.nactel == 2

    def test_lsym(self):
        assert self.ifc.lsym == 1

    def test_nisht(self):
        assert self.ifc.nisht == 7

    def test_nasht(self):
        assert self.ifc.nasht == 2

    def test_nocct(self):
        assert self.ifc.nocct == 9

    def test_norbt(self):
        assert self.ifc.norbt == 12

    def test_nbast(self):
        assert self.ifc.nbast == 12

    def test_nsym(self):
        assert self.ifc.nsym == 1

    def test_pv(self):
        pv = self.ifc.pv
        np.testing.assert_allclose(
            pv.diagonal(), [1.87279201, -0.24404616, 0.12720799], rtol=1e-7
        )
