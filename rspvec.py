#!/usr/bin/env python
"""Module for reading response vector data from RSPVEC"""

import numpy as np

def read(property_label, propfile="RSPVEC", timing=False):
    """Read response vector given property"""
    import time, numpy
    from util import full, unformatted
    if timing:
        t0 = time.clock()
    rspvec = unformatted.FortranBinary(propfile)
    rspvec.find(property_label)
    rspvec.next()
    kzyvar = rspvec.reclen / 8
    buffer_ = rspvec.readbuf(kzyvar,'d')
    mat = numpy.array(buffer_).view(full.matrix)
    if timing:
        t1 = time.clock()
        print "Time used in rspvec.read: %g" % (t1-t0)
    return mat

def tomat(N, ifc, tmpdir='/tmp'):
    """Vector to matrix"""
    import os
    from util import full
    norbt = ifc.norbt
    new = full.matrix((norbt, norbt))
    lwop = len(N)/2
    ij = 0
    for i, j in jwop(ifc):
        new[i, j] = N[ij]
        new[j, i] = N[ij+lwop]
        ij += 1
    return new

def jwop(ifc):
    """Generate  orbital excitation list"""
    for i in range(ifc.nisht):
        for j in range(ifc.nisht, ifc.norbt):
            yield (i, j)
    for i in range(ifc.nisht, ifc.nocct):
        for j in range(ifc.nocct, ifc.norbt):
            yield (i, j)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print "Usage: %s ncdim property file" % sys.argv[0]
        sys.exit(1)
    prop = sys.argv[1]
    filename = sys.argv[2]
    rvec = read(prop, filename, timing=True)
    print rvec
