import numpy as np
import numpy.testing as npt

import daltools.huckel
import daltools.prop

from . import tmpdir


class TestHuckel:

    def setup_method(self):
        self.tmpdir = tmpdir(__file__)

    def test_huckelmat(self):
        huckel, = np.array(daltools.prop.read("HUCKEL", tmpdir=self.tmpdir))
        huckel_ref = np.loadtxt(self.tmpdir/"huckel.ref")
        npt.assert_allclose(huckel,  huckel_ref)

    def test_huckelovlp(self):
        huckovlp, = np.array(daltools.prop.read("HUCKOVLP", tmpdir=self.tmpdir))
        huckovlp_ref = np.loadtxt(self.tmpdir/"huckovlp.ref")
        npt.assert_allclose(huckovlp,  huckovlp_ref)

    def test_huckel_mo(self):
        cmo = np.array(daltools.huckel.get_mo(self.tmpdir))
        cmo_ref = np.loadtxt(self.tmpdir/"huckel_mo.ref")
        npt.assert_allclose(cmo, cmo_ref)
