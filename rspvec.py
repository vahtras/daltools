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

def tovec(k, luindf="LUINDF"):
    """Matrix to vector """
    from util import full
    wop = jwop(luindf=luindf)
    lwop = len(wop)
    kzyvar = 2*lwop
    new = full.matrix((kzyvar,))
    for j in range(lwop):
        r = wop[j][0]-1
        s = wop[j][1]-1
        new[j] = k[r, s]
        new[j+lwop] = k[s, r]
    return new

def tomat(N, ifc, tmpdir='/tmp'):
    """Vector to matrix"""
    import os
    from util import full
    norbt = ifc.norbt
    new = full.matrix((norbt, norbt))
    wop = jwop(luindf=os.path.join(tmpdir, "LUINDF"))
    lwop = len(wop)
    for j in range(lwop):
        r = wop[j][0]-1
        s = wop[j][1]-1
        new[r, s] = N[j]
        new[s, r] = N[j+lwop]
    return new

def jwop(luindf="LUINDF"):
    """Get orbital excitation list from LUINDF"""
    from util import unformatted
    luindf = unformatted.FortranBinary(luindf)
    import one
    table = luindf.find("EXOPSYM1")
    rec = luindf.next()
    nwopt, = rec.read(1, 'i')
    kzy = 2*nwopt
    rec = luindf.next()
    ints = rec.read(2*nwopt,'i')
    wop = []
    for i in range(0, kzy, 2):
        wop.append((ints[i], ints[i+1]))
    return wop
      

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print "Usage: %s ncdim property file" % sys.argv[0]
        sys.exit(1)
    prop = sys.argv[1]
    filename = sys.argv[2]
    rvec = read(prop, filename, timing=True)
    print rvec
