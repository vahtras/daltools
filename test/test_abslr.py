import unittest
import os
import numpy as np
from ..rspvec import read
from ..lr import LR

class TestALR(unittest.TestCase):

    def setUp(self):
        n, _ = os.path.splitext(__file__)
        self.tmpdir = n + ".d"

    def tmp(self, file_):
        return os.path.join(self.tmpdir, file_)

    def test_first_vec(self):
        vecs = read('XDIPLEN', freqs=(0.4425,), propfile=self.tmp('ABSVECS'))
        np.testing.assert_allclose(
            vecs[('XDIPLEN', 0.4425, 0.0)][:10],
            [
             0.002520494065626863,    0.005717264921524636,
           3.8637208682448187e-11,  -1.003278319247039e-08,
             -0.08647573654795868,     -0.1961522566144851,
           -1.120999714728522e-09,  1.6920107171040828e-07,
               0.9497764738826885,      1.2962042579404978,
            ])

    def test_real_alpha(self):
        w = 0.4425
        vecs = read('XDIPLEN', freqs=(w,), propfile=self.tmp('ABSVECS'))
        re_alpha, _ = LR('XDIPLEN', 'XDIPLEN', w, tmpdir=self.tmpdir, absorption=True)
        self.assertAlmostEqual(-re_alpha, 7.096440, delta=1e-6)

    def test_real_im(self):
        w = 0.4425
        vecs = read('XDIPLEN', freqs=(w,), propfile=self.tmp('ABSVECS'))
        _, im_alpha = LR('XDIPLEN', 'XDIPLEN', w, tmpdir=self.tmpdir, absorption=True)
        self.assertAlmostEqual(-im_alpha, 0.045945, delta=1e-6)

            
if __name__ == "__main__":
    unittest.main()
