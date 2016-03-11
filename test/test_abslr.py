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

    def test_first_xvec(self):
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

    def test_first_yvec(self):
        vecs = read('XDIPLEN', 'YDIPLEN', freqs=(0.4425,), propfile=self.tmp('ABSVECS'))
        np.testing.assert_allclose(
            vecs[('YDIPLEN', 0.4425, 0.0)][:10],
            [
             0.005717270062672075,  -0.0025204990937031293,
          -2.0506203046993165e-07,  1.1264726491134807e-07,
              -0.1961521205392148,     0.08647607337762692,
            4.930847097955585e-06, -4.7055564716285557e-07,
              -1.0263058257613134,      0.8307850928096588,
            ])

    def test_number_vecs(self):
        vecs = read(
            'XDIPLEN', 'YDIPLEN', 'ZDIPLEN', freqs=(0.4425, 0.49), propfile=self.tmp('ABSVECS'),
            lr_vecs=True
            )
        self.assertEqual(len(vecs), 6)

    def test_real_alpha(self):
        w = 0.4425
        re_alpha, _ = LR('XDIPLEN', 'XDIPLEN', w, tmpdir=self.tmpdir, absorption=True)
        self.assertAlmostEqual(-re_alpha, 7.096440, delta=1e-6)

    def test_real_im(self):
        w = 0.4425
        _, im_alpha = LR('XDIPLEN', 'XDIPLEN', w, tmpdir=self.tmpdir, absorption=True)
        self.assertAlmostEqual(-im_alpha, 0.045945, delta=1e-6)

            
if __name__ == "__main__":
    unittest.main()
