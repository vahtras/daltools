#!/usr/bin/env python
"""Linear response module"""
import os
from util.timing import timing
import prop, sirifc, dens, rspvec, one
def LR(A, B, w=0.0, tmpdir='/tmp'):
    """Calculate the linear response function <<A;B>> from response vector
    on RSPVEC for B, as <[kB,A]>"""
    ifcfile = os.path.join(tmpdir, 'SIRIFC')
    propfile = os.path.join(tmpdir, 'AOPROPER')
    ifc = sirifc.sirifc(ifcfile)
    a, = prop.read(A, filename=propfile, unpack=True)
    dkb =  Dk(B, w, ifc, tmpdir)
    return a&dkb

def Dk(label, freq=0.0, ifc=None, tmpdir='/tmp'):
    """Calculate the density transformed with response vector kb
    on RSPVEC for B, as d S kb - kb S d"""
    #
    # Read interface data from SIRIFC if not provided
    #
    if ifc is None: ifc = sirifc.sirifc(name=os.path.join(tmpdir, 'SIRIFC'))
    cmo = ifc.cmo
    #
    # Get densities in AO basis
    #
    dc, do = dens.ifc(os.path.join(tmpdir, "SIRIFC"), ifc)
    d = dc+do
    #
    # Get response vector (no symmetry)
    #
    kzywop = 2*ifc.nwopt
    NB = rspvec.read(
        label, freqs=(freq,), 
        propfile=os.path.join(tmpdir, "RSPVEC")
        )[0][0]
    #
    # Vector to matrix
    #
    kB = rspvec.tomat(NB, ifc, tmpdir=tmpdir).T
    cmo = ifc.cmo.unblock()
    t_kb = timing("kB")
    kb = cmo*kB*cmo.T
    t_kb.stop()
    S = one.read(filename=os.path.join(tmpdir, 'AOONEINT')).unpack().unblock()
    #(D S kb - kb S D)
    t_dkb = timing("dkb")
    if 1:
        dS = d*S
        dkb = dS*kb-kb*dS.T
    else:
        dkb =  d*S*kb-kb*S*d
    t_dkb.stop()
    return dkb

if __name__ == "__main__":
    import sys
    from optparse import OptionParser
    op = OptionParser()
    op.add_option(
          '-t','--tmpdir',
          dest='tmpdir', default='/tmp',
          help='scratch directory [/tmp]'
          )
    try:
        opt, arg = op.parse_args(sys.argv[1:])
        Aop, Bop = arg
        tmpdir = opt.tmpdir
    except (IndexError, ValueError):
        print "Usage: %s A B" % sys.argv[0]
        sys.exit(1)

    print LR(Aop, Bop, tmpdir)
