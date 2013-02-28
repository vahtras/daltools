#!/usr/bin/env python
import struct
def read(property,propfile="RSPVEC",timing=False):
   import string,time,numpy
   from util import full,unformatted
   if timing:
      t0=time.clock()
   rspvec=unformatted.FortranBinary(propfile)
   rspvec.find(property)
   rspvec.readrec()
   kzyvar = rspvec.reclen / 8
   buffer_=rspvec.readbuf(kzyvar,'d')
   mat=numpy.array(buffer_).view(full.matrix)
   if timing:
      t1=time.clock()
      print "Time used in rspvec.read: %g"%(t1-t0)
   return mat

def tovec(k,luindf="LUINDF"):
   from util import full
   wop=jwop(luindf=luindf)
   lwop=len(wop)
   kzyvar=2*lwop
   new=full.matrix((kzyvar,))
   for j in range(lwop):
      r=wop[j][0]-1
      s=wop[j][1]-1
      new[j]=k[r,s]
      new[j+lwop]=k[s,r]
   return new

def tomat(N,ifc,tmpdir='/tmp'):
   import os
   from util import full
   norbt=ifc.norbt
   new=full.matrix((norbt,norbt))
   wop=jwop(luindf=os.path.join(tmpdir,"LUINDF"))
   lwop=len(wop)
   #print "wop.shape",N.shape,lwop
   #print "wop",wop
   for j in range(lwop):
      #print "j wop",j,wop[j]
      r=wop[j][0]-1
      s=wop[j][1]-1
      #print "r s",r,s
      new[r,s]=N[j]
      new[s,r]=N[j+lwop]
   return new

def jwop(luindf="LUINDF"):
   from util import unformatted
   luindf=unformatted.FortranBinary(luindf)
   import one
   table=luindf.find("EXOPSYM1")
   #print "table",table
   luindf.readrec()
   ints=luindf.readbuf(8,'i')
   nwopt=ints[0]
   nwoph=ints[1]
   jwopsy=ints[2]
   nklwop=ints[3]
   #print "sirset",nwopt,nwoph,jwopsy,nklwop
   kzy=2*nwopt
   #print "kzy",kzy
   luindf.readrec()
   ints=luindf.readbuf(kzy,'i')
   wop=[]
   for i in range(0,kzy,2):
      wop.append((ints[i],ints[i+1]))
   return wop
      
   return ints

if __name__ == "__main__":
   import sys,string
   if len(sys.argv) != 4:
      print "Usage: %s ncdim property file" % sys.argv[0]
      sys.exit(1)
   kzyvar=string.atoi(sys.argv[1])
   prop=sys.argv[2]
   file=sys.argv[3]
   rvec=read(kzyvar,prop,file,timing=True)
   print rvec
