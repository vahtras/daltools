#!/usr/bin/env python
"""Module for reading response vector data from RSPVEC"""

THRESHOLD = 1e-5

import numpy as np

def read(property_label, freq=0.0, propfile="RSPVEC"):
    """Read response vector given property"""
    import time, numpy
    from util import full, unformatted
    rspvec = unformatted.FortranBinary(propfile)
    while rspvec.find(property_label) is not None:
        # rspvec.rec contains matching record
        # check that frequncy matches
        veclabs = rspvec.rec.read(16,'c')
        vfreq = rspvec.rec.read(1, 'd')[0]
        if abs(vfreq-freq) < THRESHOLD:
            rspvec.next()
            kzyvar = rspvec.reclen / 8
            buffer_ = rspvec.readbuf(kzyvar,'d')
            mat = numpy.array(buffer_).view(full.matrix)
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

def tovec(mat, ifc, tmpdir='/tmp'):
    """Vector to matrix"""
    import os
    from util import full
    lwop = ifc.nisht*ifc.nasht + ifc.nocct*(ifc.norbt-ifc.nocct)
    N = full.matrix((2*lwop,))
    ij = 0
    for i, j in  jwop(ifc):
        N[ij] = mat[i, j]
        N[ij+lwop] = mat[j, i]
        ij += 1
    return N

def jwop(ifc):
    """Generate  orbital excitation list"""
    for i in range(ifc.nisht):
        for j in range(ifc.nisht, ifc.nocct):
            yield (i, j)
    for i in range(ifc.nisht):
        for j in range(ifc.nocct, ifc.norbt):
            yield (i, j)
    for i in range(ifc.nisht,ifc.nocct):
        for j in range(ifc.nocct, ifc.norbt):
            yield (i, j)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print "Usage: %s ncdim property file" % sys.argv[0]
        sys.exit(1)
    prop = sys.argv[1]
    filename = sys.argv[2]
    rvec = read(prop, filename)
    print rvec
