import os
import numpy as np
from .. import sirifc
from ..util import full

def setup():
    global tmpdir
    n, e = os.path.splitext(__file__)
    tmpdir = n + ".d"
    global ifc
    ifc = sirifc.sirifc(os.path.join(tmpdir, 'SIRIFC'))

def teardown():
    global ifc
    del(ifc)

def assert_(this, ref):
    print this, ref
    assert np.allclose(this, ref)

def test_potnuc():
    assert_(ifc.potnuc, 1.829038)

def test_emcscf():
    assert_(ifc.emcscf, -1.249293)

def test_dimensions():
    assert ifc.nisht == 0
    assert ifc.nasht == 3
    assert ifc.norbt == 3
    assert ifc.nbast == 3
    assert ifc.nsym == 4

def test_cmo():
    cmo1 = np.array([[0.45236534, 1.38834090], [0.36316622, -0.77259529]])
    assert_(ifc.cmo[0], cmo1)

def test_dv():
    dv = np.array([1.97086440, 0, 0.00940855, 0, 0, 0.01972706])
    assert_(ifc.dv.view(full.matrix), dv)
