#!/usr/bin/env python
def AB(A,B):
   """Calculate the linear response function <<A;B>> from response vector
   on RSPVEC for B, as <[kB,A]>"""
   import prop,sirifc,dens,rspvec,one,oli
   from timing import timing
   t_AB=timing("lr.AB")
   t_ifc=timing("sirifc.sirifc")
   ifc=sirifc.sirifc()
   t_ifc.stop()
   t_read=timing("prop.read")
   a=prop.read(ifc.nbast,A).unpack()
   t_read.stop()
   t_Dk=timing("lr.Dk")
   dkb= Dk(B,ifc)
   t_Dk.stop()
   if 1:
      print t_ifc
      print t_read
      print t_Dk
      print t_AB
   return a&dkb
def Dk(label,ifc=None,tmpdir='/tmp'):
   """Calculate the density transformed with response vector kb
   on RSPVEC for B, as d S kb - kb S d"""
   import os,prop,sirifc,dens,rspvec,one,oli
   from util.timing import timing 
   #
   # Read interface data from SIRIFC if not provided
   #
   if ifc is None: ifc=sirifc.sirifc(name=os.path.join(tmpdir,'SIRIFC'))
   cmo=ifc.cmo
   #
   # Get densities in AO basis
   #
   dc,do=dens.ifc(os.path.join(tmpdir,"SIRIFC"),ifc)
   d=dc+do
   #
   # Get response vector (no symmetry)
   #
   kzywop=2*ifc.nwopt
   NB=rspvec.read(label,propfile=os.path.join(tmpdir,"RSPVEC"))
   #
   # Vector to matrix
   #
   kB=oli.tomat(NB,ifc,tmpdir=tmpdir).T
   cmo=ifc.cmo.unblock()
   t_kb=timing("kB")
   kb=cmo*kB*cmo.T
   t_kb.stop()
   S=one.read(filename=os.path.join(tmpdir,'AOONEINT')).unpack().unblock()
   #(D S kb - kb S D)
   t_dkb=timing("dkb")
   if 1:
      dS=d*S
      dkb=dS*kb-kb*dS.T
   else:
      dkb= d*S*kb-kb*S*d
   t_dkb.stop()
   return dkb

if __name__ == "__main__":
   import sys,getopt
   verbose=False
   try:
      opt,arg=getopt.getopt(sys.argv[1:],'v',['verbose'])
      for o,v in opt:
	 if o in ('-v','--verbose'):
	    verbose=True
      A=arg[0]; B=arg[1]
   except IndexError:
      print "Usage: %s A B" % sys.argv[0]
      sys.exit(1)
   print AB(A,B)
