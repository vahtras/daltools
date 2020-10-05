import pathlib
import numpy


def tmpdir(f):
    f = pathlib.Path(f)
    return f.parent/(f.stem + ".d")


def assert_(this, ref):
    print(this)
    print(ref)
    assert numpy.allclose(this, ref)
