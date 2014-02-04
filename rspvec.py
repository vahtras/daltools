#!/usr/bin/env python
"""Module for reading response vector data from RSPVEC"""

THRESHOLD = 1e-5

import numpy as np
import time
from util import full, unformatted

class RspVecError(Exception): pass

def read(*args, **kwargs):
    """Read response vector given property"""
    freqs = kwargs.get('freqs', (0.0,))
    propfile = kwargs.get('propfile', 'RSPVEC')

    rspvec = unformatted.FortranBinary(propfile)
    vecs = {}

    for rec in rspvec:
        for lab in args:
            if lab in rec:
                rec.read(16,'c')
                vfreq = rec.read(1, 'd')[0]
                if vfreq in freqs:
                    rspvec.next()
                    kzyvar = rspvec.reclen / 8
                    buffer_ = rspvec.readbuf(kzyvar,'d')
                    vecs[(lab,vfreq)] = np.array(buffer_).view(full.matrix)
    # check that all required vectors are saved
    # print 'vecs',vecs.keys()
    for l in args:
        for v in freqs:
            if (l,v) not in vecs:
                raise RspVecError(
        "Linear response vector N(%s,%g) not found on file %s" %
        (l, v, propfile)
        )
    return [[vecs[(l, v)] for l in args] for v in freqs]

def readall(property_label, propfile="RSPVEC"):
    """Read response all vectors given property"""
    import time, numpy
    from util import full, unformatted
    rspvec = unformatted.FortranBinary(propfile)
    found=[]
    while rspvec.find(property_label) is not None:
        # rspvec.rec contains matching record
        veclabs = rspvec.rec.read(16,'c')
        vfreq = rspvec.rec.read(1, 'd')[0]
        rspvec.next()
        kzyvar = rspvec.reclen / 8
        buffer_ = rspvec.readbuf(kzyvar,'d')
        mat = numpy.array(buffer_).view(full.matrix)
        found.append((mat, vfreq))
    return found

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
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('prop', help='Response property')
    parser.add_argument('filename', help='File of response vectors(RSPVEC)')
    parser.add_argument('--w', type=float, default=0., help='Frequency')
    args = parser.parse_args()
    rvec = read(args.prop, propfile=args.filename, freqs=(args.w,))
    print args.prop, args.w, args.filename
    print rvec[0][0]
