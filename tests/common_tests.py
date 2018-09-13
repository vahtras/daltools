import numpy

def assert_(this,ref):
    print(this)
    print(ref)
    assert numpy.allclose(this, ref)
